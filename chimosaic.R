# =============================================================================
# Chimosaic — micro-substitution chimera / PCR artifact detection
# Taxonomy-free companion to DADA2's removeBimeraDenovo, for full-length 16S.
#
# THE SIGNAL. Two abundant sequences A and B co-amplify, and a rare sequence C
# is entirely A except for one or a few SNPs whose bases are carried
# (near-exclusively) by abundant B. removeBimeraDenovo is structurally blind to
# this: it relies on abundance-ratio heuristics and clean breakpoint detection,
# and a 1-5 SNP micro-substitution has neither.
#
# HOW IT WORKS.
#   1. No taxonomy in any decision. "Most abundant" replaces "taxonomy".
#      Nothing is auto-REAL — abundance is a SCORE, not a short-circuit.
#   2. Per-SNP donor provenance attribution: of everyone carrying that minority
#      base at that position, how exclusively does it trace to one abundant
#      donor? Exclusive is suspicious; mixed is probably a real/common base.
#   3. Fixed reference coordinate frame via cmalign against a bacterial 16S
#      covariance model (Rfam RF00177). Consensus/match columns = coordinate
#      frame; insert-state columns = indels (dropped from SNP comparison).
#      This gives coordinate stability and O(N) scaling — no global de novo MSA.
#
# CLASSES. LIKELY_ARTIFACT (crossover chimera evidence) / POINT_ERROR
# (single-parent polymerase substitution — mechanistically distinct, labelled
# separately, never silently deleted) / UNCLEAR / REAL.
#
# INPUTS.  A FASTA of ASV sequences and a CSV count table (rows = samples,
#          cols = ASVs) whose column names match the FASTA headers exactly.
# OUTPUTS. A per-ASV results table, a per-sample table, and an audit file of
#          sequences removed before scoring.
#
# REQUIRES. R (dplyr, stringr) and Infernal's cmalign on PATH, plus a
#           covariance model (RF00177.cm for bacterial 16S).
# =============================================================================

CHIMOSAIC_VERSION <- "0.1.0"

library(dplyr)
library(stringr)
# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

# --- Rarity / parent candidacy ---
RARE_THRESHOLD   <- 0.01   # rel-abundance cutoff for "rare" in a sample (<1%)
IQR_MULTIPLIER   <- 1.5    # multiplier for IQR-based per-ASV outlier detection
MAX_SNPS_CHIMERA <- 5      # >this many SNPs from nearest parent -> not a simple chimera

# --- Classification thresholds ---
LOW_ARTIFACT_THRESHOLD  <- 0.30   # artifact_score below this -> REAL
HIGH_ARTIFACT_THRESHOLD <- 0.65   # artifact_score above this -> LIKELY_ARTIFACT

# --- POINT_ERROR: a separate, mechanistic class for single-parent PCR errors ---
# NOT the same as a (bimolecular) chimera. A chimera crosses over between two
# co-amplifying parents and typically shows a RUN of donor bases (>1 SNP); its
# minority-base donor is physically co-present. A polymerase POINT substitution
# is 1 SNP off ONE hyper-abundant template, needs no second co-present parent,
# and — being regenerated at a low per-cycle rate — never blooms into a real
# fraction of any community. We label these separately (not auto-deleted) so the
# call is auditable and does not get lumped with true chimeras.
#
# A candidate is POINT_ERROR when ALL hold (and it is not already REAL):
#   - n_snps <= POINT_ERROR_MAX_SNPS               (a point substitution, not a crossover)
#   - nearest_parent reads >= K x focal reads      (an abundant template to mis-copy)
#   - focal NEVER exceeds R rel-abundance in ANY sample  (never blooms -> not a real strain)
#
# The R gate is the strain-variation safety valve (esp. for even communities like
# soil): anything 1 SNP off a parent that DOES bloom past R in some sample escapes
# this label and stays on the normal scored path. Toggle the whole rule off to
# recover 3-class behavior.
#
# TODO (sample-aware R, deferred): R is a flat fraction. A better gate would
# scale with the sample's dominant-ASV share — in a host-dominated community
# where one taxon holds 80% of reads, a genuine low-abundance strain tops out
# far lower than it would in an even community, so a flat R is comparatively
# aggressive there. Lower R for even communities (e.g. soil). Single flat knob
# for now; the hook is here so this is a small change later, not a rewrite.
POINT_ERROR_ENABLED        <- TRUE   # FALSE -> no POINT_ERROR class (3-class behavior)
POINT_ERROR_MAX_SNPS       <- 1      # max SNPs from parent to count as a point substitution
POINT_ERROR_MIN_PARENT_RATIO <- 20   # K: parent total reads must be >= this x focal total reads
POINT_ERROR_MAX_RELAB      <- 0.01   # R: if focal ever reaches this within-sample rel-abundance
                                     #    in ANY sample, it is NOT called a point error (may be a
                                     #    genuine low-abundance strain). 0.01 = 1% of a sample.

# --- Scoring weights (chimera-evidence terms; positive weights sum to 1) ---
# Provenance (exclusivity) and donor-abundance are kept as SEPARATE scores, as
# in v1. Folding abundance into the per-SNP quality as a multiplier deflated
# exclusive single-SNP calls (an abundant-enough donor is required, but a modest
# donor should not zero out a base that traces exclusively to it). W_ATTRIBUTION
# + W_DONOR_ABUNDANCE reproduce v1's 0.30/0.15 split (== the old lumped 0.45).
W_SNP               <- 0.25   # snp_count_score
W_ATTRIBUTION       <- 0.30   # mean per-SNP provenance (exclusivity of the minority base)
W_DONOR_ABUNDANCE   <- 0.15   # mean per-SNP donor abundance (is the donor a credible parent)
W_COOCCURRENCE      <- 0.20   # top donor co-occurs with the focal ASV
W_SITE_CONSERVATION <- 0.10   # parent base is near-universal at the SNP columns

# --- Abundance as a SCORE (not a short-circuit) ---
# The focal ASV's own abundance pulls its artifact_score DOWN: a sequence that
# holds a real share of the reads is unlikely to be a mere PCR artifact. This is
# how "most abundant" earns REAL status through scoring instead of a hardcoded rule.
#
# The penalty is a SOFT, SATURATING function of the ASV's GLOBAL READ SHARE
# (its reads / all reads in the dataset), NOT a normalization against the single
# most-abundant ASV and NOT a hard "top 1%" cutoff. A share-based term is
# dataset-relative, so it ports to other studies unchanged; the old
# log(reads)/log(max_reads) form was pinned to whatever the biggest ASV happened
# to be (here fd31de's 3.4M reads), which deflated exclusive single-SNP calls.
#
#   focal_ab_score = share / (share + ABUNDANCE_HALFSHARE)   # Hill, n = 1
#
# -> share >> HALFSHARE : score -> 1  (abundant; strong evidence against artifact)
# -> share == HALFSHARE : score  = 0.5 (half protection)
# -> share <<  HALFSHARE : score -> 0  (rare; no abundance protection)
W_ABUNDANCE_PENALTY <- 0.50   # 0 = ignore focal abundance; 1 = full subtraction
ABUNDANCE_HALFSHARE <- 1e-3   # global read share (0.1% of all reads) at half protection

# --- Attribution internals ---
# Donors are ranked by total reads; provenance = fraction of the minority-base
# read mass that comes from the single top donor ASV.
DONOR_MIN_READS <- 1          # a carrier must have >= this many total reads to count

# --- Overlength pre-filter (runs BEFORE cmalign) -----------------------------
# A real full-length 16S gene is ~1,500 nt. Sequences much longer than that are
# concatemers / read-through PCR products — i.e. artifacts by construction — and
# they are pathological for cmalign: the truncated-alignment DP matrix scales
# ~O(length^2), so a single 2,800 bp read needs several GB and minutes, and a
# couple in parallel OOM-kill the process (exit 137). We therefore flag any ASV
# longer than OVERLENGTH_BP, EXCLUDE it from cmalign, and classify it directly
# as LIKELY_ARTIFACT on length alone (taxonomy-free). Set OVERLENGTH_BP <- Inf
# to disable (then raise --mxsize / lower threads to survive the long ones).
OVERLENGTH_BP        <- 1600    # > this many non-gap nt -> overlength artifact
CALL_OVERLENGTH_ARTIFACT <- TRUE

# --- Input-integrity guard ---------------------------------------------------
# Mismatched FASTA headers vs. count-table column names is the most likely way
# to get a confidently wrong answer, because it silently removes the abundant
# ASVs that parent-finding depends on. Both checks are on READ MASS, not ASV
# counts — losing 10 dominant sequences matters far more than losing 500 rare
# ones. Deliberate removals (overlength, alignment-failed) are exempt.
MAX_UNMATCHED_READ_FRACTION <- 0.01  # fatal if unexplained missing ASVs exceed this share of reads
UNMATCHED_TOP_N             <- 10    # fatal if any of the top-N ASVs by reads is missing
# The reverse direction (aligned sequences with no count column) cannot be
# weighed by read mass — the counts are exactly what is missing — so it is
# gated on count instead. Raise only for a deliberate superset FASTA.
MAX_UNCOUNTED_ALN_FRACTION  <- 0.01  # fatal if this share of aligned seqs lack counts

# --- Optional abundance pre-filter (DEFAULT OFF) -----------------------------
# WARNING: this method exists to preserve rare-but-real diversity and adjudicate
# rare sequences on evidence, not to gate them by abundance. Leave these at 0
# unless you specifically need to shrink the candidate set for scalability on a
# very large dataset. Applied as an explicit, logged step after counts load.
MIN_TOTAL_READS <- 0            # drop ASVs with fewer total reads across samples
MIN_SAMPLES     <- 0            # drop ASVs detected in fewer samples

# --- cmalign coordinate frame ---
RUN_CMALIGN       <- TRUE                      # FALSE -> use provided MSA (fallback)
# Resolved from the environment so the script is portable. Override with:
#   export CMALIGN_BIN=/path/to/cmalign
#   export CM_FILE=/path/to/RF00177.cm
# or edit the unset= defaults below.
CMALIGN_BIN       <- Sys.getenv("CMALIGN_BIN", unset = unname(Sys.which("cmalign")))
CM_FILE           <- Sys.getenv("CM_FILE",     unset = "RF00177.cm")
CMALIGN_STK       <- "asvs_cmalign.stk"        # base name for per-batch output + cache
UNALIGNED_FASTA   <- "asvs_unaligned.fasta"    # degapped, overlength-filtered input

# MEMORY MODEL (learned the hard way): cmalign peak RAM has TWO independent
# drivers, and both must be bounded or it OOM-kills mid-run:
#   1. Output accumulation — it buffers the whole growing alignment until the end.
#      --matchonly caps output to the fixed 1,533 match columns (we drop inserts
#      anyway), which controls this.
#   2. Concurrent DP matrices — with --cpu N, up to N threads each allocate a DP
#      matrix (up to --mxsize MB) at once. A cluster of hard-to-band sequences on
#      16 threads x ~2 GB = ~32 GB in seconds. THIS is what kept killing us.
# We bound BOTH by aligning in BATCHES (caps accumulation + gives progress) and
# keeping threads modest (caps concurrent DP). Peak ~= CMALIGN_THREADS x mxsize.
CMALIGN_BATCH     <- 500                        # sequences per cmalign invocation
CMALIGN_THREADS   <- 2                          # concurrent DP matrices per batch
# mxsize (MB) caps a SINGLE DP matrix. Most seqs band to ~85 MB; a few divergent
# ones need several GB (one needed 4.3 GB). A seq needing more than the cap makes
# cmalign ABORT cleanly (naming the seq) -> we quarantine & flag it, not fatal.
#
# CRITICAL INVARIANT: CMALIGN_THREADS x mxsize must fit UNDER --mem, so cmalign's
# recoverable abort is the binding limit — NOT the OS OOM-killer (which is
# unrecoverable and bypasses the quarantine). Here 2 x 12000 = 24 GB, safe in a
# 32 GB job. If you raise mxsize or threads, raise --mem to match (see mem check
# below). 12000 still aligns anything needing <=12 GB (well past the 4.3 GB seen).
CMALIGN_EXTRA     <- "--noprob --dnaout --matchonly --mxsize 12000"
CMALIGN_CACHE_RDS <- paste0(CMALIGN_STK, ".matchcols.rds")  # combined frame cache
# Parsed once here so the DP cap can be reported in logs, reason strings, and
# the params output — alignment failures are only interpretable alongside it.
.MXSIZE_MB <- suppressWarnings(as.numeric(sub(".*--mxsize ([0-9]+).*", "\\1", CMALIGN_EXTRA)))

# --- Working directory -------------------------------------------------------
# Every path below — inputs, outputs, temp files, and the cmalign cache — is
# relative to this directory. Two supported ways to run:
#
#   1. cd into your data directory and run the script from wherever it lives:
#        cd /path/to/my/data && Rscript /path/to/chimera_detection_public.R
#   2. Set the working directory explicitly, and run from anywhere:
#        export CHIMOSAIC_WORK_DIR=/path/to/my/data
#      (or hardcode the unset= default below)
#
# NOTE: the cmalign coordinate frame is cached in this directory. Running from a
# different directory than a previous run will not find the cache and will
# re-align from scratch — the expensive step. Keep runs on the same dataset in
# the same directory to reuse it.
WORK_DIR <- Sys.getenv("CHIMOSAIC_WORK_DIR", unset = "")
if (nzchar(WORK_DIR)) {
  if (!dir.exists(WORK_DIR)) stop("CHIMOSAIC_WORK_DIR does not exist: ", WORK_DIR)
  setwd(WORK_DIR)
}

# --- File paths (relative to the working directory above) --------------------
# ASV sequences (aligned or unaligned — gaps are stripped before cmalign) and a
# count table whose COLUMN NAMES must match the FASTA headers exactly.
FASTA_FILE    <- "asv_sequences.fasta"
COUNTS_FILE   <- "asv_read_counts.csv"   # rows = samples, cols = ASVs
OUT_MAIN      <- "chimosaic_results.tsv"
OUT_PERSAMPLE <- "chimosaic_results_per_sample.tsv"

# Also write an .xlsx alongside the TSVs (requires the writexl package).
WRITE_XLSX    <- FALSE

# =============================================================================
# PREFLIGHT: fail fast, and fail informatively
# cmalign runs in batches, so an unset binary or a missing CM would otherwise
# surface as an opaque failure on batch 1 of N, minutes in.
# =============================================================================
cat(sprintf("chimosaic v%s\n", CHIMOSAIC_VERSION))
cat(sprintf("  Working dir: %s%s\n", getwd(),
            if (nzchar(WORK_DIR)) " (set via CHIMOSAIC_WORK_DIR)" else ""))

if (!file.exists(FASTA_FILE)) {
  stop("FASTA_FILE not found: '", FASTA_FILE, "' (working dir: ", getwd(), ")")
}
if (!file.exists(COUNTS_FILE)) {
  stop("COUNTS_FILE not found: '", COUNTS_FILE, "' (working dir: ", getwd(), ")")
}
if (RUN_CMALIGN) {
  # If a cached coordinate frame is present, cmalign is never invoked — so
  # requiring the binary and the covariance model here would block a perfectly
  # valid re-run (e.g. a scoring-only change, or a different machine that has
  # the cache but not Infernal). Check for the cache first.
  .have_cache <- file.exists(CMALIGN_CACHE_RDS)
  if (.have_cache) {
    cat(sprintf("  Cached coordinate frame found (%s) — cmalign will not run.\n",
                CMALIGN_CACHE_RDS))
  } else {
    if (!nzchar(CMALIGN_BIN)) {
      stop("cmalign not found on PATH and CMALIGN_BIN is unset.\n",
           "  Install Infernal (conda install -c bioconda infernal), or set:\n",
           "    export CMALIGN_BIN=/path/to/cmalign\n",
           "  (Not needed if a cached frame exists at ", CMALIGN_CACHE_RDS, ")")
    }
    if (!file.exists(CMALIGN_BIN)) {
      stop("CMALIGN_BIN does not exist: '", CMALIGN_BIN, "'")
    }
    if (!file.exists(CM_FILE)) {
      stop("Covariance model not found: '", CM_FILE, "'\n",
           "  CM_FILE is resolved relative to the working directory (", getwd(), ")\n",
           "  unless given as an absolute path. Set it explicitly, e.g.:\n",
           "    export CM_FILE=/absolute/path/to/RF00177.cm\n",
           "  Bacterial 16S model: Rfam RF00177 (https://rfam.org/family/RF00177).\n",
           "  (Not needed if a cached frame exists at ", CMALIGN_CACHE_RDS, ")")
    }
    cat(sprintf("  cmalign: %s\n  CM:      %s\n", CMALIGN_BIN, CM_FILE))
  }
}
if (WRITE_XLSX && !requireNamespace("writexl", quietly = TRUE)) {
  stop("WRITE_XLSX is TRUE but the 'writexl' package is not installed. ",
       "install.packages('writexl'), or set WRITE_XLSX <- FALSE.")
}

# =============================================================================
# STEP 0: BUILD THE FIXED COORDINATE FRAME WITH cmalign
# =============================================================================
# We align the (degapped) ASV sequences to a bacterial 16S covariance model.
# cmalign returns a Stockholm alignment in which MATCH (consensus) columns are
# the fixed reference coordinate frame and INSERT columns are indels relative
# to the model. We keep only the match columns, so every ASV is described in
# the SAME ~1,533-column 16S coordinate system regardless of its own indels.
# =============================================================================

parse_fasta <- function(filepath) {
  lines   <- readLines(filepath)
  headers <- which(startsWith(lines, ">"))
  seqs <- vector("character", length(headers))
  names(seqs) <- sub("^>", "", lines[headers])
  for (i in seq_along(headers)) {
    start <- headers[i] + 1
    end   <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    seqs[i] <- paste(lines[start:end], collapse = "")
  }
  seqs
}

# Parse a Stockholm/Pfam alignment (one line per sequence, as cmalign emits with
# --outformat Pfam). Returns the alignment restricted to MATCH columns, upper-
# cased, with all gaps as "-". Match columns are those where #=GC RF != '.'.
parse_stockholm_match_cols <- function(filepath) {
  lines <- readLines(filepath)
  seq_names <- character(0)
  seq_data  <- character(0)
  rf_line   <- NULL
  for (ln in lines) {
    if (startsWith(ln, "#=GC RF")) {
      rf_line <- sub("^#=GC\\s+RF\\s+", "", ln)
      next
    }
    if (startsWith(ln, "#") || startsWith(ln, "//") || !nzchar(trimws(ln))) next
    # sequence line: "name   ALIGNMENT"
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(parts) < 2) next
    nm  <- parts[1]
    seq <- paste(parts[-1], collapse = "")
    idx <- match(nm, seq_names)
    if (is.na(idx)) {                    # new sequence
      seq_names <- c(seq_names, nm)
      seq_data  <- c(seq_data, seq)
    } else {                             # wrapped block: append
      seq_data[idx] <- paste0(seq_data[idx], seq)
    }
  }
  if (is.null(rf_line)) {
    stop("No '#=GC RF' consensus annotation found in ", filepath,
         " — cannot identify match columns.")
  }
  rf_chars   <- strsplit(rf_line, "")[[1]]
  match_cols <- which(rf_chars != "." & rf_chars != "-" & rf_chars != "~")

  mat <- do.call(rbind, strsplit(seq_data, ""))
  rownames(mat) <- seq_names
  mat <- mat[, match_cols, drop = FALSE]      # fixed coordinate frame
  mat <- toupper(mat)
  mat[mat == "." | mat == "~"] <- "-"
  colnames(mat) <- as.character(seq_along(match_cols))  # consensus coordinates
  mat
}

# Returns a list: aln_mat (match-column matrix for NON-overlength seqs),
# asv_length (named non-gap length of every input seq), overlength_ids.
build_coordinate_frame <- function() {
  raw_seqs   <- parse_fasta(FASTA_FILE)
  unaligned  <- gsub("[-.~]", "", toupper(raw_seqs))
  asv_length <- nchar(unaligned)

  # --- LENGTH PRE-FILTER (runs first, right after removeBimera) --------------
  # A 16S gene is ~1,500 nt. Anything much longer is physically impossible as a
  # single 16S (concatemer / read-through) — not a diversity judgement call, so
  # we remove it up front and REPORT the count before proceeding. These commonly
  # slip past removeBimera (which models a single two-parent join on exact seqs).
  overlength_ids <- names(asv_length)[asv_length > OVERLENGTH_BP]
  keep_ids       <- names(asv_length)[asv_length <= OVERLENGTH_BP]
  cat(sprintf("  [Length pre-filter] %d of %d sequences > %d nt flagged as probable errors (median length %d nt).\n",
              length(overlength_ids), length(asv_length), OVERLENGTH_BP,
              as.integer(median(asv_length))))
  if (length(overlength_ids) > 0) {
    cat(sprintf("    Longest: %d nt. These are set aside before alignment; see overlength_removed.csv.\n",
                max(asv_length)))
  }

  if (!RUN_CMALIGN) {
    cat("  RUN_CMALIGN = FALSE -> using provided MSA as coordinate frame.\n")
    m <- do.call(rbind, strsplit(toupper(raw_seqs[keep_ids]), ""))
    rownames(m) <- keep_ids
    m[m == "." | m == "~"] <- "-"
    return(list(aln_mat = m, asv_length = asv_length, overlength_ids = overlength_ids,
                aln_failed_ids = character(0)))
  }

  # Reuse the combined coordinate frame if we've already built it.
  if (file.exists(CMALIGN_CACHE_RDS)) {
    cat("  Reusing cached coordinate frame: ", CMALIGN_CACHE_RDS, "\n")
    m <- readRDS(CMALIGN_CACHE_RDS)
    kept <- rownames(m)
    return(list(aln_mat = m, asv_length = asv_length, overlength_ids = overlength_ids,
                aln_failed_ids = base::setdiff(keep_ids, kept)))
  }

  # Align in batches. Every batch is aligned to the SAME CM, so its match columns
  # are the same fixed 1,533 model positions regardless of batch membership —
  # after we keep only match columns, the per-batch matrices stack (rbind) into
  # one coordinate frame. Batching caps both memory drivers (see MEMORY MODEL).
  kseq    <- unaligned[keep_ids]
  n       <- length(kseq)
  n_batch <- ceiling(n / CMALIGN_BATCH)
  cat(sprintf("  Aligning %d seqs in %d batch(es) of <=%d (cmalign --cpu %d, %s)...\n",
              n, n_batch, CMALIGN_BATCH, CMALIGN_THREADS, CMALIGN_EXTRA))
  log_file <- paste0(CMALIGN_STK, ".log")
  if (file.exists(log_file)) file.remove(log_file)
  bat_fa  <- paste0(UNALIGNED_FASTA, ".batch")
  bat_stk <- paste0(CMALIGN_STK, ".batch")
  mats <- vector("list", n_batch)
  # Sequences cmalign could not align within the --mxsize DP-matrix cap. NOTE:
  # this is a RESOURCE limit, not a verdict from the model — a larger --mxsize
  # would align some of them. On 16S-only input a sequence that will not band
  # against a 16S model is nonetheless suspect, so we flag these as artifacts by
  # default, but the reason string records the true trigger and the cap used so
  # the call stays auditable and reproducible.
  aln_failed_ids <- character(0)

  for (b in seq_len(n_batch)) {
    lo <- (b - 1) * CMALIGN_BATCH + 1
    hi <- min(b * CMALIGN_BATCH, n)
    batch_ids <- names(kseq)[lo:hi]

    # Align this batch; if cmalign aborts because a specific sequence's DP matrix
    # exceeds --mxsize (it names the sequence), quarantine that sequence and retry.
    # On 16S-only input these are strong artifact candidates, so we flag them
    # rather than let one sequence kill the run — but see the note above: the
    # trigger is a resource cap, so the cap is recorded with every such call.
    repeat {
      if (length(batch_ids) == 0) { mats[[b]] <- NULL; break }
      bseq <- kseq[batch_ids]
      fa <- character(2 * length(bseq))
      fa[c(TRUE, FALSE)] <- paste0(">", names(bseq))
      fa[c(FALSE, TRUE)] <- bseq
      writeLines(fa, bat_fa)

      cmd <- sprintf("%s --cpu %d %s --outformat Pfam -o %s %s %s",
                     shQuote(CMALIGN_BIN), CMALIGN_THREADS, CMALIGN_EXTRA,
                     shQuote(bat_stk), shQuote(CM_FILE), shQuote(bat_fa))
      out    <- suppressWarnings(system(cmd, intern = TRUE))
      status <- attr(out, "status"); if (is.null(status)) status <- 0L
      cat(out, file = log_file, sep = "\n", append = TRUE)

      ok <- status == 0 && file.exists(bat_stk) && file.size(bat_stk) > 0
      if (ok) {
        mats[[b]] <- parse_stockholm_match_cols(bat_stk)
        cat(sprintf("    batch %d/%d: %d seqs aligned (%d match cols)%s\n",
                    b, n_batch, nrow(mats[[b]]), ncol(mats[[b]]),
                    if (length(batch_ids) < (hi - lo + 1))
                      sprintf("; %d quarantined", (hi - lo + 1) - length(batch_ids)) else ""))
        break
      }

      # Try to recover: parse the offending sequence name(s) from cmalign's error.
      bad <- unique(sub(".*sequence[: ]+(\\S+).*", "\\1",
                        grep("Problem during alignment of sequence", out, value = TRUE)))
      bad <- base::intersect(bad, batch_ids)
      is_mem_abort <- any(grepl("DP mxes need|mxsize", out))
      if (length(bad) == 0 || !is_mem_abort) {
        err <- grep("Error:|Problem", out, value = TRUE)
        stop("cmalign failed on batch ", b, "/", n_batch, " (status ", status,
             if (status == 137) " = SIGKILL / out-of-memory: lower CMALIGN_THREADS, --mxsize, or CMALIGN_BATCH" else "",
             "). ",
             if (length(err)) paste0("cmalign says: ", paste(err, collapse = " | "), ". ") else "",
             "Full log: ", log_file, "\n",
             "  Checked: CMALIGN_BIN='", CMALIGN_BIN, "', CM_FILE='", CM_FILE, "'.")
      }
      aln_failed_ids <- c(aln_failed_ids, bad)
      batch_ids <- base::setdiff(batch_ids, bad)
      cat(sprintf("    batch %d/%d: quarantined %s (DP matrix > --mxsize %s MB); retrying %d seqs...\n",
                  b, n_batch, paste(bad, collapse = ", "),
                  .MXSIZE_MB, length(batch_ids)))
    }
  }
  if (file.exists(bat_fa))  file.remove(bat_fa)
  if (file.exists(bat_stk)) file.remove(bat_stk)
  if (length(aln_failed_ids) > 0) {
    cat(sprintf("  [Alignment] %d sequence(s) failed to align within --mxsize %s MB; flagged as artifacts.\n",
                length(aln_failed_ids), .MXSIZE_MB))
    cat("    (Resource limit, not a model verdict — raising --mxsize may align some. See removed_sequences.csv.)\n")
  }

  mats <- mats[!vapply(mats, is.null, logical(1))]
  # All batches must share the same fixed coordinate frame.
  ncols <- vapply(mats, ncol, integer(1))
  if (length(unique(ncols)) != 1) {
    stop("Batches produced differing match-column counts (", paste(unique(ncols), collapse = ", "),
         ") — cannot combine into one coordinate frame.")
  }
  m <- do.call(rbind, mats)
  saveRDS(m, CMALIGN_CACHE_RDS)
  cat(sprintf("  Combined coordinate frame: %d seqs x %d cols -> cached in %s\n",
              nrow(m), ncol(m), CMALIGN_CACHE_RDS))
  list(aln_mat = m, asv_length = asv_length, overlength_ids = overlength_ids,
       aln_failed_ids = aln_failed_ids)
}

# --- Memory-safety guard -----------------------------------------------------
# cmalign's real peak is MORE than CMALIGN_THREADS x mxsize: mxsize caps only the
# main DP matrix, but each thread ALSO holds HMM banding matrices, and the batch's
# growing alignment + base overhead add on top. Empirically (a 2 x 12000 = 24 GB
# nominal config OOM-killed a 32 GB job) the true multiplier is ~1.5x. We budget
# with CMALIGN_MEM_FACTOR so cmalign's recoverable mxsize-abort — not the OS
# OOM-killer, which bypasses the quarantine — is always the binding limit.
CMALIGN_MEM_FACTOR <- 1.5
if (RUN_CMALIGN) {
  .mxsize_mb <- .MXSIZE_MB
  .slurm_mem <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")))  # MB, if under SLURM
  .peak <- CMALIGN_MEM_FACTOR * CMALIGN_THREADS * .mxsize_mb
  if (!is.na(.mxsize_mb) && !is.na(.slurm_mem) && .slurm_mem > 0) {
    if (.peak > 0.85 * .slurm_mem) {
      stop(sprintf(paste0(
        "Unsafe cmalign config: est peak %d MB (%.1f x %d threads x %d MB mxsize) exceeds 85%% of job memory (%d MB).\n",
        "  This risks an OS OOM-kill that bypasses the quarantine backstop.\n",
        "  Fix: lower CMALIGN_THREADS or --mxsize, or raise --mem to >= %d MB."),
        as.integer(.peak), CMALIGN_MEM_FACTOR, CMALIGN_THREADS, as.integer(.mxsize_mb),
        as.integer(.slurm_mem), as.integer(ceiling(.peak / 0.85))))
    }
    cat(sprintf("  [Mem check] cmalign peak est %.0f GB (%.1f x %d threads x %.0f GB) vs job mem %.0f GB — OK\n",
                .peak/1024, CMALIGN_MEM_FACTOR, CMALIGN_THREADS, .mxsize_mb/1024, .slurm_mem/1024))
  } else if (!is.na(.mxsize_mb)) {
    cat(sprintf("  [Mem check] cmalign peak est ~%.0f GB (%.1f x %d threads x %.0f GB); ensure --mem exceeds this.\n",
                .peak/1024, CMALIGN_MEM_FACTOR, CMALIGN_THREADS, .mxsize_mb/1024))
  }
}

cat("Building fixed coordinate frame (cmalign / RF00177)...\n")
.cf             <- build_coordinate_frame()
aln_mat         <- .cf$aln_mat
asv_length      <- .cf$asv_length        # non-gap length of every input seq
overlength_ids  <- .cf$overlength_ids    # excluded from alignment; flagged as artifacts
aln_failed_ids <- .cf$aln_failed_ids   # exceeded the DP cap; flagged as artifacts
if (is.null(aln_failed_ids)) aln_failed_ids <- character(0)
aln_len         <- ncol(aln_mat)
cat(sprintf("  %d aligned sequences x %d consensus columns\n", nrow(aln_mat), aln_len))

# =============================================================================
# STEP 1: LOAD READ COUNTS
# =============================================================================
cat("Loading read counts...\n")
counts_raw  <- read.csv(COUNTS_FILE, check.names = FALSE, row.names = 1)
counts_mat  <- as.matrix(counts_raw)          # rows = samples, cols = ASVs
storage.mode(counts_mat) <- "double"

# Sequences removed BEFORE the chimera-evidence step because they're structurally
# implausible or unalignable — overlength (concatemers) and sequences that would
# not band within the DP cap. Both are flagged, removed from the scored set, and recorded together in
# ONE consolidated audit file (removed_sequences.csv) so anyone can hand-check the
# full set of exclusions in one place. Per-set stat vectors are also kept for the
# output table below.
overlength_in_counts  <- base::intersect(overlength_ids,  colnames(counts_mat))
aln_failed_in_counts <- base::intersect(aln_failed_ids, colnames(counts_mat))

col_stats <- function(ids) list(
  total = colSums(counts_mat[, ids, drop = FALSE]),
  max   = apply(counts_mat[, ids, drop = FALSE], 2, max),
  nsamp = colSums(counts_mat[, ids, drop = FALSE] > 0))

removed_rows <- list()
if (length(overlength_in_counts) > 0) {
  s <- col_stats(overlength_in_counts)
  ol_total <- s$total; ol_max <- s$max; ol_nsamp <- s$nsamp
  removed_rows[["overlength"]] <- data.frame(
    asv_id      = overlength_in_counts,
    filter      = "overlength",
    length_bp   = as.integer(asv_length[overlength_in_counts]),
    total_reads = as.integer(ol_total[overlength_in_counts]),
    n_samples   = as.integer(ol_nsamp[overlength_in_counts]),
    reason      = sprintf("length %d nt > %d (probable concatemer/read-through)",
                          as.integer(asv_length[overlength_in_counts]), OVERLENGTH_BP),
    stringsAsFactors = FALSE, row.names = NULL)
}
if (length(aln_failed_in_counts) > 0) {
  s <- col_stats(aln_failed_in_counts)
  un_total <- s$total; un_max <- s$max; un_nsamp <- s$nsamp
  removed_rows[["alignment_failed"]] <- data.frame(
    asv_id      = aln_failed_in_counts,
    filter      = "alignment_failed",
    length_bp   = as.integer(asv_length[aln_failed_in_counts]),
    total_reads = as.integer(un_total[aln_failed_in_counts]),
    n_samples   = as.integer(un_nsamp[aln_failed_in_counts]),
    reason      = sprintf("cmalign DP matrix exceeded --mxsize %s MB (resource limit, not a model verdict); would not band against the 16S CM",
                          .MXSIZE_MB),
    stringsAsFactors = FALSE, row.names = NULL)
}
if (length(removed_rows) > 0) {
  removed_df <- bind_rows(removed_rows) %>% arrange(filter, desc(length_bp))
  write.csv(removed_df, "removed_sequences.csv", row.names = FALSE)
  cat(sprintf("  [Pre-filter] %d sequences removed & flagged for hand-check (%d overlength, %d alignment-failed) -> removed_sequences.csv\n",
              nrow(removed_df), length(overlength_in_counts), length(aln_failed_in_counts)))
}

# Keep only ASVs present in BOTH the alignment and the count table.
#
# GUARD (read-mass, not ASV count). A mismatch between FASTA headers and
# count-table column names surfaces here as a silently smaller intersection, not
# an error: the run completes and reports confident calls on whatever survived.
#
# Counting ASVs is the wrong test. The damaging case is losing a SMALL NUMBER OF
# ABUNDANT sequences — parent candidacy depends entirely on them, so dropping a
# handful of dominant ASVs guts the method while the total ASV loss still looks
# unremarkable. (Pairing a compartment-subset count table with a full FASTA does
# exactly this: ~89% of ASVs survive, but the dominant parents do not.) We
# therefore gate on the READ MASS that goes missing and on whether any top-N ASV
# is among the casualties.
#
# Sequences we removed on purpose (overlength, alignment-failed) are expected to
# be absent from the alignment and are excluded from the accounting.
.counts_totals    <- colSums(counts_mat)
.grand_raw        <- sum(.counts_totals)
.expected_missing <- base::union(overlength_in_counts, aln_failed_in_counts)
.unexplained      <- base::setdiff(
                       base::setdiff(colnames(counts_mat), rownames(aln_mat)),
                       .expected_missing)
.lost_mass        <- if (length(.unexplained) > 0)
                       sum(.counts_totals[.unexplained]) / .grand_raw else 0
.top_ids          <- names(sort(.counts_totals, decreasing = TRUE))[
                       seq_len(min(UNMATCHED_TOP_N, length(.counts_totals)))]
.top_missing      <- base::intersect(.top_ids, .unexplained)

common_asvs <- base::intersect(colnames(counts_mat), rownames(aln_mat))

.hdr_hint <- paste0(
  "  Count table columns: ", ncol(counts_mat), " (e.g. ",
  paste(utils::head(colnames(counts_mat), 2), collapse = ", "), ")\n",
  "  Alignment rows:      ", nrow(aln_mat), " (e.g. ",
  paste(utils::head(rownames(aln_mat), 2), collapse = ", "), ")\n",
  "  FASTA headers must match count-table column names exactly.")

if (length(common_asvs) == 0) {
  stop("No ASVs shared between the count table and the alignment.\n", .hdr_hint)
}
if (length(.top_missing) > 0) {
  stop(sprintf(paste0(
    "%d of the %d most abundant ASVs in the count table are missing from the alignment.\n",
    "  Missing: %s\n",
    "  These are exactly the sequences parent-finding depends on; results would be\n",
    "  wrong rather than merely incomplete, so this is fatal.\n",
    "  Usual cause: the count table and FASTA describe different ASV sets (e.g. a\n",
    "  subset count table against a full FASTA), or headers were altered.\n%s"),
    length(.top_missing), UNMATCHED_TOP_N,
    paste(utils::head(.top_missing, 5), collapse = ", "), .hdr_hint))
}
if (.lost_mass > MAX_UNMATCHED_READ_FRACTION) {
  stop(sprintf(paste0(
    "ASVs missing from the alignment hold %.2f%% of all reads (limit %.2f%%).\n",
    "  %d ASV(s) unaccounted for; largest: %s\n",
    "  Raise MAX_UNMATCHED_READ_FRACTION only if you know why they are absent.\n%s"),
    100 * .lost_mass, 100 * MAX_UNMATCHED_READ_FRACTION,
    length(.unexplained),
    paste(utils::head(.unexplained[order(-.counts_totals[.unexplained])], 5),
          collapse = ", "),
    .hdr_hint))
}
if (length(.unexplained) > 0) {
  cat(sprintf("  [Note] %d ASV(s) in the count table are absent from the alignment (%.3f%% of reads) — dropped.\n",
              length(.unexplained), 100 * .lost_mass))
}
.aln_only <- base::setdiff(rownames(aln_mat), colnames(counts_mat))
# The OTHER direction, and the more dangerous one: sequences we aligned that have
# no column in the count table. This is the failure that pairing a subset count
# table with a full FASTA actually produces — and it is undetectable by read mass,
# because the missing ASVs have no counts to weigh. Their abundance is precisely
# what we cannot see. So the threshold here is on COUNT, and it is deliberately
# tight: if the FASTA and the count table were built from the same run they should
# describe the same ASV set, and any sizeable divergence means one of them is the
# wrong file. Raise the tolerance only if you know the FASTA is a deliberate
# superset (e.g. a shared reference FASTA across several count tables).
.aln_only_frac <- length(.aln_only) / max(1, nrow(aln_mat))
if (.aln_only_frac > MAX_UNCOUNTED_ALN_FRACTION) {
  stop(sprintf(paste0(
    "%d of %d aligned sequences (%.1f%%) have no column in the count table (limit %.1f%%).\n",
    "  Examples: %s\n",
    "  Their read counts are unknown, so this cannot be checked by abundance — but if\n",
    "  any of them are dominant ASVs, parent-finding silently loses the very sequences\n",
    "  the method depends on and every downstream call is wrong.\n",
    "  Usual cause: a subset count table (e.g. one compartment) paired with a full FASTA.\n",
    "  If the FASTA is an intentional superset, raise MAX_UNCOUNTED_ALN_FRACTION.\n%s"),
    length(.aln_only), nrow(aln_mat), 100 * .aln_only_frac,
    100 * MAX_UNCOUNTED_ALN_FRACTION,
    paste(utils::head(.aln_only, 5), collapse = ", "), .hdr_hint))
}
if (length(.aln_only) > 0) {
  cat(sprintf("  [Note] %d aligned sequence(s) have no column in the count table — dropped.\n",
              length(.aln_only)))
}

counts_mat  <- counts_mat[, common_asvs, drop = FALSE]
aln_mat     <- aln_mat[common_asvs, , drop = FALSE]
all_asvs    <- colnames(counts_mat)
all_samples <- rownames(counts_mat)
cat(sprintf("  %d samples x %d ASVs (intersection of counts and alignment)\n",
            nrow(counts_mat), ncol(counts_mat)))

sample_totals     <- rowSums(counts_mat)
asv_totals        <- colSums(counts_mat)
asv_max_reads     <- apply(counts_mat, 2, max)
asv_n_samples_det <- colSums(counts_mat > 0)
max_total_reads   <- max(asv_totals, na.rm = TRUE)
grand_total_reads <- sum(asv_totals, na.rm = TRUE)   # all reads, all ASVs, all samples
# Peak within-sample relative abundance per ASV (best sample). Sample-depth-robust
# "does it ever bloom" signal used by the POINT_ERROR rule. Empty samples -> 0.
asv_max_relab     <- apply(counts_mat / pmax(sample_totals, 1), 2, max)

# --- Optional abundance pre-filter (default off; see WARNING at config) -------
if (MIN_TOTAL_READS > 0 || MIN_SAMPLES > 0) {
  keep <- asv_totals >= MIN_TOTAL_READS & asv_n_samples_det >= MIN_SAMPLES
  cat(sprintf("  [Abundance pre-filter] MIN_TOTAL_READS=%d, MIN_SAMPLES=%d -> dropping %d of %d ASVs.\n",
              MIN_TOTAL_READS, MIN_SAMPLES, sum(!keep), length(keep)))
  cat("    NOTE: abundance filtering removes rare candidates this method is meant to adjudicate.\n")
  counts_mat <- counts_mat[, keep, drop = FALSE]
  aln_mat    <- aln_mat[colnames(counts_mat), , drop = FALSE]
  all_asvs   <- colnames(counts_mat)
  sample_totals     <- rowSums(counts_mat)
  asv_totals        <- colSums(counts_mat)
  asv_max_reads     <- apply(counts_mat, 2, max)
  asv_n_samples_det <- colSums(counts_mat > 0)
  max_total_reads   <- max(asv_totals, na.rm = TRUE)
  grand_total_reads <- sum(asv_totals, na.rm = TRUE)
  asv_max_relab     <- apply(counts_mat / pmax(sample_totals, 1), 2, max)
}

# =============================================================================
# STEP 2: RARE vs. ABUNDANT (per sample) — no taxonomy, no priors
# =============================================================================
cat("Classifying rare vs. abundant ASVs per sample...\n")
rel_abund_mat  <- sweep(counts_mat, 1, sample_totals, "/")
rare_fixed_mat <- rel_abund_mat < RARE_THRESHOLD & counts_mat > 0

rare_iqr_mat <- matrix(FALSE, nrow(counts_mat), ncol(counts_mat),
                       dimnames = dimnames(counts_mat))
for (asv in all_asvs) {
  nz <- counts_mat[counts_mat[, asv] > 0, asv]
  if (length(nz) < 4) next
  lv  <- log10(nz + 0.5)
  q1  <- quantile(lv, 0.25); q3 <- quantile(lv, 0.75)
  fence <- q1 - IQR_MULTIPLIER * (q3 - q1)
  rare_iqr_mat[names(nz)[lv < fence], asv] <- TRUE
}

is_rare_fixed_any     <- apply(rare_fixed_mat, 2, any)
is_rare_iqr_any       <- apply(rare_iqr_mat,   2, any)
is_rare_consensus     <- is_rare_fixed_any | is_rare_iqr_any
is_abundant_somewhere <- apply(!rare_fixed_mat & counts_mat > 0, 2, any)  # parent candidacy

cat(sprintf("  Rare (either threshold) in >=1 sample: %d | Parent candidates: %d\n",
            sum(is_rare_consensus), sum(is_abundant_somewhere)))

# =============================================================================
# STEP 3: PER-SAMPLE PARENT SEARCH (nearest abundant co-occurring ASV)
# =============================================================================
count_snps <- function(a, b) {
  valid <- a != "-" & b != "-"
  diffs <- a[valid] != b[valid]
  list(n_snps = sum(diffs), positions = which(valid)[diffs])
}

cat("Running per-sample parent search...\n")
rare_asv_ids          <- all_asvs[is_rare_consensus]
parent_candidates_ids <- all_asvs[is_abundant_somewhere]
parent_mat <- aln_mat[parent_candidates_ids, , drop = FALSE]

results_list <- vector("list", length(rare_asv_ids))
names(results_list) <- rare_asv_ids
pb_step <- max(1, floor(length(rare_asv_ids) / 20))

for (ri in seq_along(rare_asv_ids)) {
  if (ri %% pb_step == 0) cat(sprintf("  ... %d/%d\n", ri, length(rare_asv_ids)))
  rare_id    <- rare_asv_ids[ri]
  rare_chars <- aln_mat[rare_id, ]

  rare_in_samples <- rownames(counts_mat)[
    (rare_fixed_mat[, rare_id] | rare_iqr_mat[, rare_id]) & counts_mat[, rare_id] > 0]
  if (length(rare_in_samples) == 0) next

  sample_results <- vector("list", length(rare_in_samples))
  names(sample_results) <- rare_in_samples
  for (samp in rare_in_samples) {
    cands <- parent_candidates_ids[
      counts_mat[samp, parent_candidates_ids] > 0 & parent_candidates_ids != rare_id]
    if (length(cands) == 0) {
      sample_results[[samp]] <- list(parent = NA, n_snps = NA, snp_cols = integer(0)); next
    }
    best_parent <- NA_character_; best_snps <- Inf; best_cols <- integer(0)
    for (cand in cands) {
      res <- count_snps(rare_chars, parent_mat[cand, ])
      if (res$n_snps < best_snps) {
        best_snps <- res$n_snps; best_parent <- cand; best_cols <- res$positions
      }
    }
    sample_results[[samp]] <- list(parent = best_parent, n_snps = best_snps,
                                   snp_cols = best_cols)
  }

  parents_found <- sapply(sample_results, function(x) x$parent)
  snps_found    <- sapply(sample_results, function(x) x$n_snps)
  if (all(is.na(parents_found))) {
    results_list[[ri]] <- list(
      rare_id = rare_id, nearest_parent = NA_character_, n_snps_vs_parent = NA_integer_,
      snp_aln_cols = integer(0), n_samples_rare = length(rare_in_samples),
      n_samples_parent_found = 0L,
      parent_search_note = "No abundant co-occurring parent found in any sample")
    next
  }
  parent_freq <- sort(table(parents_found[!is.na(parents_found)]), decreasing = TRUE)
  top_parent  <- names(parent_freq)[1]
  snps_top    <- snps_found[!is.na(parents_found) & parents_found == top_parent]
  first_samp  <- names(parents_found)[!is.na(parents_found) & parents_found == top_parent][1]
  results_list[[ri]] <- list(
    rare_id = rare_id, nearest_parent = top_parent,
    n_snps_vs_parent = as.integer(round(median(snps_top))),
    snp_aln_cols = sample_results[[first_samp]]$snp_cols,
    n_samples_rare = length(rare_in_samples),
    n_samples_parent_found = sum(!is.na(parents_found)),
    parent_search_note = sprintf("Parent '%s' found in %d/%d samples with median %d SNPs",
                                 top_parent, sum(!is.na(parents_found)),
                                 length(rare_in_samples),
                                 as.integer(round(median(snps_top)))))
}
cat("  Parent search complete.\n")

# =============================================================================
# STEP 4: DONOR ATTRIBUTION
#   per-SNP quality = donor_abundance x provenance  (reported; the two terms are
#   also scored SEPARATELY in Step 5 — folding them into one multiplier deflated
#   exclusive single-SNP calls, because a modest-but-exclusive donor should not
#   be zeroed out by the abundance factor)
#     donor_abundance = log10(top_donor_reads+1) / log10(max_reads+1)
#     provenance      = top_donor_reads / (total reads of ALL carriers of the
#                       minority base at that column)  [1 = exclusive donor]
# =============================================================================
cat("Attributing SNPs to donors (donor_abundance x provenance)...\n")

log_max <- log10(max_total_reads + 1)

for (ri in seq_along(results_list)) {
  r <- results_list[[ri]]
  if (is.null(r) || is.na(r$nearest_parent) || length(r$snp_aln_cols) == 0) {
    if (!is.null(r)) {
      results_list[[ri]]$n_snps_attributed    <- 0L
      results_list[[ri]]$fraction_attributed  <- NA_real_
      results_list[[ri]]$mean_snp_quality     <- 0
      results_list[[ri]]$mean_provenance      <- 0
      results_list[[ri]]$mean_donor_abundance <- 0
      results_list[[ri]]$top_donor            <- NA_character_
      results_list[[ri]]$top_donor_total_reads<- NA_real_
      results_list[[ri]]$snp_positions_cons   <- NA_character_
      results_list[[ri]]$snp_changes          <- NA_character_
      results_list[[ri]]$per_snp_quality      <- NA_character_
    }
    next
  }

  rare_id   <- r$rare_id
  snp_cols  <- r$snp_aln_cols
  rare_chars<- aln_mat[rare_id, ]
  par_chars <- aln_mat[r$nearest_parent, ]

  snp_changes_vec <- paste0(par_chars[snp_cols], "→", rare_chars[snp_cols])
  per_snp_q  <- numeric(length(snp_cols))
  provenance <- numeric(length(snp_cols))
  donor_ab   <- numeric(length(snp_cols))
  snp_top_donor <- character(length(snp_cols))
  attributed <- logical(length(snp_cols))

  for (si in seq_along(snp_cols)) {
    col <- snp_cols[si]; X <- rare_chars[col]
    carriers <- rownames(aln_mat)[aln_mat[, col] == X & rownames(aln_mat) != rare_id]
    carriers <- carriers[carriers %in% names(asv_totals)]
    carriers <- carriers[asv_totals[carriers] >= DONOR_MIN_READS]
    if (length(carriers) == 0) next
    attributed[si] <- TRUE
    creads   <- asv_totals[carriers]
    top_id   <- carriers[which.max(creads)]
    top_reads<- creads[which.max(creads)]
    provenance[si]   <- as.numeric(top_reads / sum(creads))       # exclusivity
    donor_ab[si]     <- as.numeric(log10(top_reads + 1) / log_max) # abundance
    per_snp_q[si]    <- provenance[si] * donor_ab[si]
    snp_top_donor[si]<- top_id
  }

  n_attr <- sum(attributed)
  # Report the ASV that is the top donor for the most SNP columns (tiebreak: reads)
  if (n_attr > 0) {
    tbl  <- table(snp_top_donor[attributed])
    tied <- names(tbl)[tbl == max(tbl)]
    top_donor_id <- tied[which.max(asv_totals[tied])]
  } else {
    top_donor_id <- NA_character_
  }

  results_list[[ri]]$n_snps_attributed     <- n_attr
  results_list[[ri]]$fraction_attributed   <- n_attr / length(snp_cols)
  results_list[[ri]]$mean_snp_quality      <- if (length(snp_cols) > 0) mean(per_snp_q) else 0
  # Separated components (scored independently, v1-style):
  results_list[[ri]]$mean_provenance       <- if (length(snp_cols) > 0) mean(provenance) else 0
  results_list[[ri]]$mean_donor_abundance  <- if (length(snp_cols) > 0) mean(donor_ab)  else 0
  results_list[[ri]]$top_donor             <- top_donor_id
  results_list[[ri]]$top_donor_total_reads <- if (!is.na(top_donor_id)) asv_totals[top_donor_id] else NA_real_
  results_list[[ri]]$snp_positions_cons    <- paste(snp_cols, collapse = ";")
  results_list[[ri]]$snp_changes           <- paste(snp_changes_vec, collapse = ";")
  results_list[[ri]]$per_snp_quality       <- paste(sprintf("%.3f", per_snp_q), collapse = ";")
}
cat("  Donor attribution complete.\n")

# =============================================================================
# STEP 5: SCORING — abundance is a penalty score, not a short-circuit
# =============================================================================
cat("Computing artifact scores...\n")

for (ri in seq_along(results_list)) {
  r <- results_list[[ri]]
  if (is.null(r)) next
  n_snps <- r$n_snps_vs_parent

  # Score 1: SNP-count (1 SNP -> 1.0, MAX_SNPS -> 0.0)
  snp_score <- if (is.na(n_snps) || n_snps > MAX_SNPS_CHIMERA) 0 else
    max(0, 1 - (n_snps - 1) / (MAX_SNPS_CHIMERA - 1))

  # Score 2: mean per-SNP provenance (exclusivity of the minority base)
  att_score <- if (is.na(n_snps) || n_snps > MAX_SNPS_CHIMERA) 0 else r$mean_provenance

  # Score 2b: mean per-SNP donor abundance (separate term, v1-style)
  dab_score <- if (is.na(n_snps) || n_snps > MAX_SNPS_CHIMERA) 0 else r$mean_donor_abundance

  # Score 3: co-occurrence of top donor with the focal ASV
  rare_id <- r$rare_id; top_donor <- r$top_donor
  co_score <- 0
  if (!is.na(top_donor) && top_donor %in% colnames(counts_mat)) {
    rp <- counts_mat[, rare_id] > 0; dp <- counts_mat[, top_donor] > 0
    if (sum(rp) > 0) co_score <- sum(rp & dp) / sum(rp)
  }

  # Score 4: site conservation of the parent base at the SNP columns
  site_score <- 0
  if (length(r$snp_aln_cols) > 0 && !is.na(r$nearest_parent)) {
    par_chars <- aln_mat[r$nearest_parent, ]
    cs <- numeric(length(r$snp_aln_cols))
    for (si in seq_along(r$snp_aln_cols)) {
      col <- r$snp_aln_cols[si]; pb <- par_chars[col]
      cb  <- aln_mat[, col]; valid <- cb != "-"
      cs[si] <- if (sum(valid) > 0) sum(cb[valid] == pb) / sum(valid) else 0
    }
    site_score <- mean(cs)
  }

  # Focal-abundance penalty score (abundance treated as evidence AGAINST artifact).
  # Soft, saturating function of GLOBAL READ SHARE — dataset-relative, no top-N
  # cutoff, and not pinned to the single most-abundant ASV. See config notes.
  focal_share    <- as.numeric(asv_totals[rare_id]) / grand_total_reads
  focal_ab_score <- focal_share / (focal_share + ABUNDANCE_HALFSHARE)

  evidence <- W_SNP * snp_score + W_ATTRIBUTION * att_score +
              W_DONOR_ABUNDANCE * dab_score +
              W_COOCCURRENCE * co_score + W_SITE_CONSERVATION * site_score
  composite <- evidence - W_ABUNDANCE_PENALTY * focal_ab_score
  composite <- min(1, max(0, composite))

  # >MAX_SNPS cannot be a simple chimera -> zero everything (routes to REAL)
  if (!is.na(n_snps) && n_snps > MAX_SNPS_CHIMERA) {
    snp_score <- 0; att_score <- 0; dab_score <- 0; co_score <- 0; site_score <- 0; composite <- 0
  }

  results_list[[ri]]$snp_count_score         <- snp_score
  results_list[[ri]]$attribution_score       <- att_score
  results_list[[ri]]$donor_abundance_score   <- dab_score
  results_list[[ri]]$cooccurrence_score      <- co_score
  results_list[[ri]]$site_conservation_score <- site_score
  results_list[[ri]]$focal_read_share        <- as.numeric(focal_share)
  results_list[[ri]]$focal_abundance_score   <- as.numeric(focal_ab_score)
  # POINT_ERROR inputs (see classify_asv): focal peak within-sample rel-abundance,
  # focal total reads, and the abundant backbone parent's total reads.
  results_list[[ri]]$focal_max_relab         <- as.numeric(asv_max_relab[rare_id])
  results_list[[ri]]$focal_total_reads       <- as.numeric(asv_totals[rare_id])
  results_list[[ri]]$parent_total_reads      <- if (!is.na(r$nearest_parent) &&
                                                     r$nearest_parent %in% names(asv_totals))
                                                  as.numeric(asv_totals[r$nearest_parent]) else NA_real_
  results_list[[ri]]$artifact_score          <- composite
}
cat("  Scoring complete.\n")

# =============================================================================
# STEP 6: CLASSIFICATION — nothing auto-REAL; everything is scored
# =============================================================================
cat("Classifying ASVs...\n")

classify_asv <- function(r) {
  n_snps    <- r$n_snps_vs_parent
  art       <- r$artifact_score
  frac      <- r$fraction_attributed
  has_parent<- !is.na(r$nearest_parent)

  if (!has_parent || is.na(art)) {
    return(list(class = "UNCLEAR", reason = sprintf(
      "No abundant co-occurring parent in any of %d rare sample(s); cannot evaluate origin",
      r$n_samples_rare)))
  }
  if (!is.na(n_snps) && n_snps > MAX_SNPS_CHIMERA && art < LOW_ARTIFACT_THRESHOLD) {
    return(list(class = "REAL", reason = sprintf(
      "%d SNPs from nearest parent (>%d); artifact score %.2f — likely genuine variant",
      n_snps, MAX_SNPS_CHIMERA, art)))
  }
  if (art < LOW_ARTIFACT_THRESHOLD) {
    return(list(class = "REAL", reason = sprintf(
      "Low artifact score %.2f (< %.2f); abundance/evidence favor genuine",
      art, LOW_ARTIFACT_THRESHOLD)))
  }
  # POINT_ERROR: single-parent polymerase substitution, NOT a chimera. 1 SNP off
  # an abundant backbone template that never blooms into a real fraction anywhere.
  # Placed AFTER the REAL branches so a genuinely abundant sequence still wins.
  if (POINT_ERROR_ENABLED && !is.na(n_snps) && n_snps <= POINT_ERROR_MAX_SNPS) {
    prat  <- if (!is.na(r$parent_total_reads) && !is.na(r$focal_total_reads) &&
                 r$focal_total_reads > 0) r$parent_total_reads / r$focal_total_reads else NA_real_
    relab <- r$focal_max_relab
    if (!is.na(prat) && prat >= POINT_ERROR_MIN_PARENT_RATIO &&
        !is.na(relab) && relab < POINT_ERROR_MAX_RELAB) {
      return(list(class = "POINT_ERROR", reason = sprintf(
        "%d SNP from abundant parent '%s' (%.0fx focal reads); peak %.3f%% of any sample (< %.2f%%) — single-parent PCR/seq error, never blooms",
        n_snps, r$nearest_parent, prat, relab * 100, POINT_ERROR_MAX_RELAB * 100)))
    }
  }
  if (!is.na(n_snps) && n_snps <= MAX_SNPS_CHIMERA &&
      art >= HIGH_ARTIFACT_THRESHOLD && !is.na(frac) && frac > 0) {
    return(list(class = "LIKELY_ARTIFACT", reason = sprintf(
      "%d SNP(s) from parent '%s'; %.0f%% SNPs attributed; mean SNP quality %.2f; artifact score %.2f",
      n_snps, r$nearest_parent, frac * 100, r$mean_snp_quality, art)))
  }
  list(class = "UNCLEAR", reason = sprintf(
    "Artifact score %.2f in ambiguous band [%.2f, %.2f); %s SNPs; %.0f%% attributed",
    art, LOW_ARTIFACT_THRESHOLD, HIGH_ARTIFACT_THRESHOLD,
    if (is.na(n_snps)) "unknown" else as.character(n_snps),
    if (is.na(frac)) 0 else frac * 100))
}

for (ri in seq_along(results_list)) {
  r <- results_list[[ri]]; if (is.null(r)) next
  cl <- classify_asv(r)
  results_list[[ri]]$classification        <- cl$class
  results_list[[ri]]$classification_reason <- cl$reason
}

# =============================================================================
# STEP 7: BUILD MAIN TABLE
# =============================================================================
cat("Building output tables...\n")
safe_get <- function(lst, key, default = NA) {
  v <- lst[[key]]; if (is.null(v)) default else v
}

main_rows <- lapply(results_list, function(r) {
  if (is.null(r)) return(NULL)
  id <- r$rare_id
  data.frame(
    asv_id                  = id,
    classification          = safe_get(r, "classification", "UNCLEAR"),
    classification_reason   = safe_get(r, "classification_reason", NA_character_),
    artifact_score          = safe_get(r, "artifact_score", NA_real_),
    n_snps_vs_parent        = safe_get(r, "n_snps_vs_parent", NA_integer_),
    nearest_parent          = safe_get(r, "nearest_parent", NA_character_),
    n_snps_attributed       = safe_get(r, "n_snps_attributed", NA_integer_),
    fraction_attributed     = safe_get(r, "fraction_attributed", NA_real_),
    mean_snp_quality        = safe_get(r, "mean_snp_quality", NA_real_),
    top_donor               = safe_get(r, "top_donor", NA_character_),
    top_donor_total_reads   = safe_get(r, "top_donor_total_reads", NA_real_),
    snp_positions_cons      = safe_get(r, "snp_positions_cons", NA_character_),
    snp_changes             = safe_get(r, "snp_changes", NA_character_),
    per_snp_quality         = safe_get(r, "per_snp_quality", NA_character_),
    snp_count_score         = safe_get(r, "snp_count_score", NA_real_),
    attribution_score       = safe_get(r, "attribution_score", NA_real_),
    donor_abundance_score   = safe_get(r, "donor_abundance_score", NA_real_),
    cooccurrence_score      = safe_get(r, "cooccurrence_score", NA_real_),
    site_conservation_score = safe_get(r, "site_conservation_score", NA_real_),
    focal_read_share        = safe_get(r, "focal_read_share", NA_real_),
    focal_abundance_score   = safe_get(r, "focal_abundance_score", NA_real_),
    peak_within_sample_relab= safe_get(r, "focal_max_relab", NA_real_),
    n_samples_rare          = safe_get(r, "n_samples_rare", NA_integer_),
    n_samples_parent_found  = safe_get(r, "n_samples_parent_found", NA_integer_),
    parent_search_note      = safe_get(r, "parent_search_note", NA_character_),
    stringsAsFactors = FALSE)
})
main_df <- bind_rows(main_rows[!sapply(main_rows, is.null)])

# Overlength ASVs (excluded from cmalign): classified by length alone. A 16S
# gene is ~1,500 nt, so a much longer sequence is a concatemer/read-through
# artifact. This is a taxonomy-free, alignment-free artifact call.
if (length(overlength_in_counts) > 0) {
  ol_class  <- if (CALL_OVERLENGTH_ARTIFACT) "LIKELY_ARTIFACT" else "UNCLEAR"
  main_df <- bind_rows(main_df, data.frame(
    asv_id                = overlength_in_counts,
    classification        = ol_class,
    classification_reason = sprintf(
      "Overlength: %d nt (> %d); ~2x a 16S gene, not alignable to the model — probable concatemer/chimera",
      as.integer(asv_length[overlength_in_counts]), OVERLENGTH_BP),
    artifact_score        = 1,
    n_snps_vs_parent      = NA_integer_,
    stringsAsFactors      = FALSE))
}

# Alignment-failed ASVs: cmalign could not band these within the --mxsize DP
# cap. On 16S-only input, a sequence that will not align to a 16S model is a
# strong artifact candidate, so they are flagged on the same footing as
# overlength. The trigger is a resource limit rather than a model verdict, so
# the cap is recorded in the reason string and the params file — a run with a
# different --mxsize may resolve some of these differently.
if (length(aln_failed_in_counts) > 0) {
  main_df <- bind_rows(main_df, data.frame(
    asv_id                = aln_failed_in_counts,
    classification        = if (CALL_OVERLENGTH_ARTIFACT) "LIKELY_ARTIFACT" else "UNCLEAR",
    classification_reason = sprintf(
      "Alignment failed: cmalign DP matrix exceeded --mxsize %s MB; would not band against the 16S CM — probable non-target/chimeric",
      .MXSIZE_MB),
    artifact_score        = 1,
    n_snps_vs_parent      = NA_integer_,
    stringsAsFactors      = FALSE))
}

# ASVs never flagged rare: scored implicitly as REAL (abundant everywhere they
# occur). They still carry NO taxonomy-based shortcut — they are REAL because
# nothing rarer needs explaining.
unevaluated <- base::setdiff(all_asvs, main_df$asv_id)
if (length(unevaluated) > 0) {
  main_df <- bind_rows(main_df, data.frame(
    asv_id                = unevaluated,
    classification        = "REAL",
    classification_reason = "Not rare in any sample; no chimera evidence to evaluate",
    artifact_score        = 0,
    n_snps_vs_parent      = NA_integer_,
    stringsAsFactors      = FALSE))
}

# Length + structural-exclusion flags for every row
main_df$length_bp      <- as.integer(asv_length[main_df$asv_id])
main_df$is_overlength  <- main_df$asv_id %in% overlength_ids
main_df$is_alignment_failed <- main_df$asv_id %in% aln_failed_ids

# Shared per-ASV stats
main_df$total_reads          <- as.numeric(asv_totals[main_df$asv_id])
main_df$max_reads_any_sample <- as.numeric(asv_max_reads[main_df$asv_id])
main_df$n_samples_detected   <- as.integer(asv_n_samples_det[main_df$asv_id])
main_df$is_rare_fixed        <- as.logical(is_rare_fixed_any[main_df$asv_id])
main_df$is_rare_outlier      <- as.logical(is_rare_iqr_any[main_df$asv_id])
main_df$is_abundant_somewhere<- as.logical(is_abundant_somewhere[main_df$asv_id])

# Overlength rows aren't in the aligned/subset stat vectors — fill from the
# counts recorded before subsetting, and mark them non-rare / non-parent.
if (length(overlength_in_counts) > 0) {
  oi <- match(overlength_in_counts, main_df$asv_id)
  main_df$total_reads[oi]           <- as.numeric(ol_total[overlength_in_counts])
  main_df$max_reads_any_sample[oi]  <- as.numeric(ol_max[overlength_in_counts])
  main_df$n_samples_detected[oi]    <- as.integer(ol_nsamp[overlength_in_counts])
  main_df$is_rare_fixed[oi]         <- FALSE
  main_df$is_rare_outlier[oi]       <- FALSE
  main_df$is_abundant_somewhere[oi] <- FALSE
}
if (length(aln_failed_in_counts) > 0) {
  ui <- match(aln_failed_in_counts, main_df$asv_id)
  main_df$total_reads[ui]           <- as.numeric(un_total[aln_failed_in_counts])
  main_df$max_reads_any_sample[ui]  <- as.numeric(un_max[aln_failed_in_counts])
  main_df$n_samples_detected[ui]    <- as.integer(un_nsamp[aln_failed_in_counts])
  main_df$is_rare_fixed[ui]         <- FALSE
  main_df$is_rare_outlier[ui]       <- FALSE
  main_df$is_abundant_somewhere[ui] <- FALSE
}

main_df <- main_df %>%
  mutate(sort_order = case_when(classification == "LIKELY_ARTIFACT" ~ 1,
                                classification == "POINT_ERROR"     ~ 2,
                                classification == "UNCLEAR"         ~ 3,
                                TRUE                                ~ 4)) %>%
  arrange(sort_order, desc(artifact_score)) %>% select(-sort_order)

# =============================================================================
# STEP 8: PER-SAMPLE TABLE
# =============================================================================
persample_rows <- list()
for (samp in all_samples) {
  ids <- colnames(counts_mat)[counts_mat[samp, ] > 0]
  if (length(ids) == 0) next
  persample_rows[[samp]] <- data.frame(
    sample_id          = samp,
    asv_id             = ids,
    reads              = counts_mat[samp, ids],
    relative_abundance = counts_mat[samp, ids] / sample_totals[samp],
    is_rare_fixed      = rare_fixed_mat[samp, ids],
    is_rare_outlier    = rare_iqr_mat[samp, ids],
    stringsAsFactors   = FALSE)
}
persample_df <- bind_rows(persample_rows) %>%
  left_join(main_df %>% select(asv_id, classification, artifact_score,
                               nearest_parent, n_snps_vs_parent, top_donor),
            by = "asv_id") %>%
  arrange(sample_id, desc(reads))

# =============================================================================
# STEP 9: WRITE OUTPUT
# =============================================================================
cat("Writing output files...\n")

# Every parameter that can change a call is recorded, so a results file is
# self-describing and a run is reproducible from the table alone.
params_df <- data.frame(
  Parameter = c("chimosaic_version","run_date",
                "work_dir","FASTA_FILE","COUNTS_FILE",
                "RARE_THRESHOLD","IQR_MULTIPLIER","MAX_SNPS_CHIMERA",
                "LOW_ARTIFACT_THRESHOLD","HIGH_ARTIFACT_THRESHOLD",
                "W_SNP","W_ATTRIBUTION","W_DONOR_ABUNDANCE","W_COOCCURRENCE","W_SITE_CONSERVATION",
                "W_ABUNDANCE_PENALTY","ABUNDANCE_HALFSHARE","DONOR_MIN_READS",
                "POINT_ERROR_ENABLED","POINT_ERROR_MAX_SNPS","POINT_ERROR_MIN_PARENT_RATIO",
                "POINT_ERROR_MAX_RELAB",
                "OVERLENGTH_BP","CALL_OVERLENGTH_ARTIFACT",
                "MIN_TOTAL_READS","MIN_SAMPLES",
                "MAX_UNMATCHED_READ_FRACTION","UNMATCHED_TOP_N","MAX_UNCOUNTED_ALN_FRACTION",
                "RUN_CMALIGN","CM_FILE","CMALIGN_EXTRA","CMALIGN_BATCH","CMALIGN_THREADS",
                "n_samples","n_asvs_scored"),
  Value = as.character(c(CHIMOSAIC_VERSION, format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            getwd(), FASTA_FILE, COUNTS_FILE,
            RARE_THRESHOLD, IQR_MULTIPLIER, MAX_SNPS_CHIMERA,
            LOW_ARTIFACT_THRESHOLD, HIGH_ARTIFACT_THRESHOLD,
            W_SNP, W_ATTRIBUTION, W_DONOR_ABUNDANCE, W_COOCCURRENCE, W_SITE_CONSERVATION,
            W_ABUNDANCE_PENALTY, ABUNDANCE_HALFSHARE, DONOR_MIN_READS,
            POINT_ERROR_ENABLED, POINT_ERROR_MAX_SNPS, POINT_ERROR_MIN_PARENT_RATIO,
            POINT_ERROR_MAX_RELAB,
            OVERLENGTH_BP, CALL_OVERLENGTH_ARTIFACT,
            MIN_TOTAL_READS, MIN_SAMPLES,
            MAX_UNMATCHED_READ_FRACTION, UNMATCHED_TOP_N, MAX_UNCOUNTED_ALN_FRACTION,
            RUN_CMALIGN, CM_FILE, CMALIGN_EXTRA, CMALIGN_BATCH, CMALIGN_THREADS,
            nrow(counts_mat), nrow(main_df))),
  stringsAsFactors = FALSE)

# TSV is the default so output drops straight into a pipeline; xlsx is opt-in.
write.table(main_df,     OUT_MAIN,      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(persample_df, OUT_PERSAMPLE, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
OUT_PARAMS <- sub("\\.tsv$", "", OUT_MAIN)
OUT_PARAMS <- paste0(OUT_PARAMS, "_params.tsv")
write.table(params_df, OUT_PARAMS, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

if (WRITE_XLSX) {
  xlsx_path <- paste0(sub("\\.tsv$", "", OUT_MAIN), ".xlsx")
  writexl::write_xlsx(list("Results"    = main_df,
                           "Per-sample" = persample_df,
                           "Params"     = params_df),
                      path = xlsx_path)
  cat(sprintf("  Also wrote %s\n", xlsx_path))
}

cat(sprintf("\n=== DONE (chimosaic v%s) ===\n", CHIMOSAIC_VERSION))
cat(sprintf("Main results:     %s\n", OUT_MAIN))
cat(sprintf("Per-sample table: %s\n", OUT_PERSAMPLE))
cat(sprintf("Run parameters:   %s\n", OUT_PARAMS))
if (length(removed_rows) > 0) cat("Removed (audit):  removed_sequences.csv\n")
cat("\nClassification summary:\n"); print(table(main_df$classification))
cat("\nTop LIKELY_ARTIFACT ASVs:\n")
main_df %>% filter(classification == "LIKELY_ARTIFACT") %>%
  select(asv_id, artifact_score, n_snps_vs_parent, mean_snp_quality,
         nearest_parent, top_donor) %>% head(20) %>% print()
