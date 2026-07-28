# Chimosaic
Micro-substitution chimera / PCR artifact detection for PacBio full-length 16S rRNA amplicon sequencing; to be run on sequence tables after DADA2 or other denoisers. 

A taxonomy-free tool for detecting micro-substitution chimeras in full-length 16S data, designed to run alongside DADA2's removeBimeraDenovo. 

> **Status: v0.1.0 — early release.** The method is described in Bedwell et al 2026, in review (BioRxiv DOI: https://doi.org/10.64898/2025.12.02.691918). Formal benchmarking against simulated and mock-community data is in progress.


---

## The problem

In low-complexity, host-dominated communities (within plants, insect guts, and similar systems where one or two taxa hold most of the reads, and particularly with PacBio long-read amplicons that add length and read depth), a specific artifact appears.

Two abundant sequences **A** and **B** co-amplify. A rare sequence **C** shows up that is *entirely* A except for one to five SNPs, and the bases at those SNPs are carried, near-exclusively, by abundant B. 

This falls outside the bimera model. Tools such as removeBimeraDenovo target two-parent chimeras with a detectable breakpoint, but a 1–5 SNP micro-substitution offers neither signal: there is no breakpoint to find, and the sequence sits too close to its parent for an abundance-ratio test to separate it from a genuine rare variant.

Chimosaic asks, for each differing position: **who else carries this exact base, and how exclusively does it trace back to one abundant donor?** A base only ever carried by an abundant donor is suspicious. 

Chimosaic is meant to run alongside current chimera and sequence error removal tools such as removeBimeraDenovo. Currently it is tooled to run at the end of DADA2 or similar pipelines on polished sequence tables, but this can be toggled to suit the user. 


## How it works

1. **No taxonomy in any decision.** "Most abundant" replaces "taxonomy" throughout. Nothing is automatically called real; abundance is scored for all. 
2. **Per-SNP donor provenance attribution.** For each SNP, of all sequences carrying that minority base, what fraction of their reads belongs to the single most abundant carrier?
3. **A fixed reference coordinate frame** built with `cmalign` against a bacterial 16S covariance model (Rfam RF00177). Consensus/match columns define the coordinate system; insert-state columns are treated as indels and dropped from SNP comparison. This gives coordinate stability across samples and O(N) scaling and does not require intensive de-novo multiple sequence alignment. 

## Requirements

- **R** with `dplyr` and `stringr` (`writexl` only if you want `.xlsx` output)
- **[Infernal](http://eddylab.org/infernal/)** — `cmalign` must be on your `PATH`
- **A covariance model** — [Rfam RF00177](https://rfam.org/family/RF00177) for bacterial 16S

```bash
conda install -c bioconda infernal
```

## Installation

It's a single R script. Clone or download it.

```bash
git clone https://github.com/bedwellbacteria/Chimosaic.git
```

## Quick start

Put your FASTA and count table in one directory, then either `cd` there or point `CHIMOSAIC_WORK_DIR` at it:

```bash
export CM_FILE=/absolute/path/to/RF00177.cm
cd /path/to/my/data
Rscript /path/to/chimosaic.R
```

or

```bash
export CHIMOSAIC_WORK_DIR=/path/to/my/data
export CM_FILE=/absolute/path/to/RF00177.cm
Rscript /path/to/chimosaic.R
```

All other settings live in the config block at the top of the script. Edit `FASTA_FILE` and `COUNTS_FILE` to match your filenames.

### Inputs

| File | Format |
|---|---|
| `asvs.fasta` | ASV sequences, aligned or unaligned (gaps are stripped before `cmalign`) |
| `asv_read_counts.csv` | Rows = samples, columns = ASVs, first column = sample IDs |

**The count table's column names must match the FASTA headers exactly.** This is the most common way to get a confidently wrong answer, because a mismatch silently drops the abundant sequences that parent-finding depends on. The script checks for this and refuses to run if the overlap looks wrong — see [Input integrity](#input-integrity).

### Outputs

| File | Contents |
|---|---|
| `chimosaic_results.tsv` | One row per ASV: classification, score, parent, donor, per-SNP detail |
| `chimosaic_results_per_sample.tsv` | Per-sample abundances joined to classifications |
| `chimosaic_results_params.tsv` | Every parameter used, so a results file is self-describing |
| `removed_sequences.csv` | Audit trail for sequences excluded before scoring |

## Classifications

| Class | Meaning |
|---|---|
| `LIKELY_ARTIFACT` | Evidence of crossover: few SNPs from an abundant parent, with those SNPs attributable to an abundant donor |
| `POINT_ERROR` | A single substitution from a much more abundant parent that never blooms in any sample — a polymerase/sequencing error, mechanistically distinct from a chimera |
| `UNCLEAR` | Ambiguous score band, or no abundant co-occurring parent to compare against |
| `REAL` | No artifact evidence |

`POINT_ERROR` exists as its own class deliberately. These are single-parent errors, not template-switching products, and collapsing them into a chimera call would misattribute the mechanism. They are labelled, never silently deleted — you decide what to do with them.

## Key parameters

All in the config block at the top of the script.

| Parameter | Default | What it controls |
|---|---|---|
| `RARE_THRESHOLD` | `0.01` | Within-sample relative abundance below which an ASV is a candidate |
| `MAX_SNPS_CHIMERA` | `5` | Maximum SNPs from a parent still considered a micro-substitution |
| `HIGH_ARTIFACT_THRESHOLD` | `0.65` | Score at or above which an ASV is called `LIKELY_ARTIFACT` |
| `LOW_ARTIFACT_THRESHOLD` | `0.30` | Score below which an ASV is called `REAL` |
| `POINT_ERROR_MAX_RELAB` | `0.01` | Peak within-sample abundance a point error may never exceed |
| `POINT_ERROR_MIN_PARENT_RATIO` | `20` | How many times more abundant the parent must be |
| `OVERLENGTH_BP` | `1600` | Length above which a sequence is called a concatemer (tied to the 16S CM) |

The scoring weights (`W_SNP`, `W_ATTRIBUTION`, `W_DONOR_ABUNDANCE`, `W_COOCCURRENCE`, `W_SITE_CONSERVATION`, `W_ABUNDANCE_PENALTY`) are also exposed. They are currently set by reasoning about the mechanism rather than learned from labelled data — treat them as a starting point, and record any changes as a sensitivity analysis rather than a tuning result.

## Input integrity

Two guards run before scoring, both deliberately fatal:

- **Abundant ASVs missing from the alignment.** Gated on read mass, plus a check on whether any of the top-N ASVs by reads is absent.
- **Aligned sequences with no count column.** Gated on count.
  
Both thresholds are configurable (`MAX_UNMATCHED_READ_FRACTION`, `UNMATCHED_TOP_N`, `MAX_UNCOUNTED_ALN_FRACTION`) but the defaults are tight on purpose.

## Performance notes

The cmalign step requires significant memory and should be run on a cluster or serer. Once the coordinate frame is built, the remaining steps run comfortably with local compute.

The coordinate frame is cached to `asvs_cmalign.stk.matchcols.rds` in the working directory. Later runs reuse it and skip `cmalign` entirely — including on machines without Infernal installed. Keep runs on the same dataset in the same directory to benefit.

`cmalign` buffers output across all sequences before writing, so memory scales with dataset size and multiplies with `--cpu`. The script processes in batches (`CMALIGN_BATCH`) and estimates peak memory before starting; under SLURM it will refuse to run a configuration likely to be OOM-killed. If you hit memory trouble, lower `CMALIGN_THREADS` or `--mxsize` before lowering the batch size.

Single-SNP chimeras are classified as POINT_ERROR. The point-error rule is evaluated first and POINT_ERROR_MAX_SNPS defaults to 1, so a genuine 1-SNP template-switch product that meets the abundance criteria lands there. The classes are not cleanly separable at one SNP.

Point errors introducing a novel base have the possibility to escape to REAL. If no other ASV carries the new base, nothing can be attributed, the score falls below LOW_ARTIFACT_THRESHOLD, and the sequence is called real before the point-error rule is reached. 


## Citation

If you use Chimosaic, please cite the preprint:

> *Bedwell et al. Concurrent ecological and evolutionary processes contribute to mutualism breakdown between legumes and rhizobia. bioRxiv (2025). doi:10.64898/2025.12.02.691918*

The specific version of this script used in the preprint is available via Zenodo: *10.5281/zenodo.17410935*
This current version has been updated to be taxonomy free, and produces comparable results on the same data. 

## License

University of Illinois/NCSA Open Source License. See the license for more information. 

Copyright (c) 2026 The Board of Trustees of the University of Illinois.

