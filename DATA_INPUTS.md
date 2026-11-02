# Raw data inputs

Every input `extractionScripts/` reads, where it comes from, and whether it is
ours to deposit. Compiled 2026-10-28 by walking the submission wrappers and
their sample CSVs, then checking each path on disk.

Sizes are what the directory holds today on `/projects/b1042/GoyalLab/egrody`.

## 1. Our sequencing — deposit to SRA

| Accession | Path under `rawData/` | Figure | Libraries | Size |
|---|---|---|---|---|
| EGS016 | `EGS016/fastq_concat` | 1 | 24 Ca\*Re amplicon (env, SSenv, pol, D1 × D0, W0, W2, D195, Invitro, Pacute) | 233G |
| EGS018 | `EGS018/Sequencing/fastq_concat` | 1 | 6 Ca\*Re *tat* amplicon, one per sample | 11G |
| EGS026 | `EGS026/fastq_concat/tenx` | 3 | 3 10x GEX (SI-GA-A3/C3/D3 = HIV-A/C/D) | 129G |
| EGS026 | `EGS026/fastq_concat/dualIndex` | 3 | 15 Ca\*Re amplicon (D1, D4, pol, env, nef × 3 clades) | 137G |
| EGS026 | `EGS026/fastq_P4/tenx` | 3 | 10x GEX, P4 reseq | 123G |
| EGS007 | `EGS007/Sequencing/longRead` | 2 | long-read cDNA + CIe amplicons, Invitro only | 13G |

**≈646G across 6 directories.** `EGS018/Sequencing/fastq_concat_tat`, which the
Figure 1 CSV names, is a directory of symlinks into `fastq_concat` — deposit the
real files, not the link farm. `EGS026/fastq_concat` and `fastq_P4` hold more
than the subdirectories listed above (396G and 370G in total); only the
subdirectories the wrappers read are needed.

## 2. Our sequencing — FASTQs no longer on disk

These six produced data the pipeline still depends on, but their FASTQ
directories are empty or gone: zero `.fastq.gz`, no broken symlinks, the first
five emptied 2026-02-01 (EGL002 on 2026-06-01). For all but LD002 the derived
cellranger output survives.

| Accession | Figure | What it was | Derived output that survives |
|---|---|---|---|
| EGS002 | 1 | 10x GEX: W0, W2, Invitro | `extractedData/Figure1/singleCell/counts/` |
| EGS004 | 1 | 10x GEX: LP29_D0, LP29_D195 | same |
| EGS009 | 1 | 10x GEX: P_acute_combined | same |
| EGL001 | 4 | 10x GEX: Normal_SHAM24_R2, Esrp2KO_sham24_R1 (mouse liver) | `extractedData/Figure4/liver/singleCell/counts/` |
| EGL002 | 4 | Ca\*Re: 12 host-target liver libraries | `extractedData/Figure4/liver/CaRe/counts/nointron/` |
| LD002 | 4 | bulk RNA-seq: 15 single-cell-derived clones, 2 lanes each (`20240925_LD02/Emmie/`) | `extractedData/LD002/nf-core/` |

**This blocks both goals.** SRA takes reads, so nothing here can be deposited
until the FASTQs are recovered from backup, the sequencing core, or whoever
holds the liver data. And the pipeline cannot be rerun from the top for these
arms — `runAll_TenX.sh` says so in its own header:

> There is no stepOne here: the TenX cellranger runs are the original EGS002 /
> EGS004 / EGS009 counts, not something this pipeline produces.

Figures 1 and 4 can still be rerun **from stepTwo**, which consumes the
surviving BAMs.

LD002's FASTQ directory `20240925_LD02/` no longer exists, but its nf-core
output was restored to `extractedData/LD002/nf-core/` on 30 October 2026 — 259G
including the per-clone BAMs, of which the two salmon matrices
(`star_salmon/salmon.merged.gene_{tpm,counts}.tsv`) are all Figure 4 Panel J
needs. `memoryseq_variability.sh` regenerates
`potential_heritable_genes_filtered.csv` from them and reproduces the published
252-row file exactly, so the only thing the missing reads block is the SRA
deposit.

## 3. Published data — cite the accession, deposit nothing

All three download themselves; no upload obligation.

| Dataset | Accession | Runs | Downloader |
|---|---|---|---|
| Collora et al. 2022 | PRJNA777601 / GSE187515 | 56 | `publicData/download_collora_sra.sh` |
| Giron et al. 2024 | PRJNA1070656 / SRP486539 / GSE254483 | 8 | `publicData/download_ena.sh` |
| HIV-Seq | GSE266455 | 35 | `publicData/download_sra.sh`, `download_GSE266455.sh` |

Giron also needs `publicData/Gironetal2024/feature_ref/feature_reference.csv`,
the hashing feature-barcode reference. It is 512 bytes and ours — include it in
the repository rather than the data deposit.

## 4. References — cite, do not deposit

| Reference | Used by | Status |
|---|---|---|
| `genomes/Mmul_10_mac239` (cellranger, Mmul_10 + SIVmac239) | Fig 1, Fig 2 polyA | present, 35G |
| `genomes/GRCh38_HXB2` (cellranger, GRCh38 + HIV HXB2) | Fig 3, all public data | present, 35G |
| `genomes/mac239` (viral only) | Fig 2 polyA | present, 178K |
| `genomes/SE_events_Mmul_10_targets.gff3` | Fig 4 BRIE2 | present, 2.6K |
| `refdata-gex-GRCh38-2024-A` (10x, plain GRCh38) | Fig 4 LD002 | **absent** — public 10x download |

The references are not deposited; the methods describe how each was built. That
is the normal arrangement -- nobody deposits a cellranger reference -- and no
new mkref scripts are needed for the release. Two things the prose has to carry
for that to be enough:

- **Pin the mac239 contig exactly.** Every coordinate in
  `pipelineScripts/bamsort/coordinates/` is valid only against the 9602 bp
  mac239 contig as built under cellranger 7.2.0, and the Ca\*Re and TenX arms
  are comparable only because both align to that same contig. Give the
  accession, the length, and how it was appended to Mmul_10, or the published
  coordinate tables cannot be used.
- **Pin the source builds**: Mmul_10 and GRCh38 release versions, and the HXB2
  accession.

The LD002 reference is the stock 10x GRCh38 build, not one of ours, and it is
the only reference a wrapper names that is not under `genomes/`. It is a public
download, so nothing is lost; `memoryseq_nfcore_rnaseq.sh` fails with the URL
rather than guessing, and its header gives `GRCh38_HXB2/` as the on-disk
substitute. The substitution is safe for an uninfected human clone panel — no
read maps to the HXB2 contig — but TPM is computed over whichever annotation is
supplied, so whichever build is used has to be the one the methods name.

`submissionScripts/PublicData/HIV/mkref_GRCh38_HXB2.sh` stays in the repository.
It is already written and a working script reproduces a reference more exactly
than prose can; `Mmul_10_mac239` is described in the methods instead.

## 5. Imaging — not SRA

| Accession | Figure | Contents | Size |
|---|---|---|---|
| EGS023 | 2 | smFISH TIFFs, 3 dates + `masks/` | 188G |
| EGS027 | 2 | smFISH TIFFs, 2 dates + `masks/` | 281G |

469G of image data. SRA does not accept it; it needs an image repository (BioStudies /
BioImage Archive) or a generalist archive such as Zenodo, and the segmentation
masks should travel with the TIFFs since `bigfish_manual.sh` reads both.

## 6. Intermediates with no regeneration path

Derived data the pipeline consumes but does not produce, beyond the stepOne gap
in §2:

- `extractedData/Figure2/longRead/nf-core/cDNA/bam/dedup/cDNA.dedup.sorted.bam`
  — the nf-core long-read pipeline that made it is not in this repository.
- `extractedData/Figure4/singleCell/counts/concat/` (GFP_invitro, 20G) — no
  stepOne wrapper.
- `extractedData/LD002/nf-core/` is the restored output of an nf-core/rnaseq run
  whose FASTQs are gone (§2). `memoryseq_nfcore_rnaseq.sh` is in the repository
  and would reproduce it given the reads, and everything downstream of it runs
  from the salmon matrices, so this is a gap in the deposit rather than in the
  pipeline.
- The liver cellranger `outs/` hold only `filtered_` and `raw_feature_bc_matrix`.
  There is no `possorted_genome_bam.bam`, so the liver GEX arm cannot be
  re-bamsorted; `liver_stepTwo_bamsort.sh` reads the Ca\*Re counts under
  `Figure4/liver/CaRe/counts/nointron/` instead, and those are present.

## Portability note

`submissionScripts/Figure1/stepOne_mkcounts.csv` hardcodes absolute
`/projects/b1042/GoyalLab/egrody/...` paths in all 30 rows and the wrapper uses
them verbatim, so `DATA_ROOT` does not redirect it — unlike
`stepOne_mkcounts_KLRB1.csv`, whose `Fastq_Path` is relative to `DATA_ROOT` and
is joined at `stepOne_mkcounts_KLRB1.sh:81`. Make the Figure 1 CSV relative
before release, or the repository is not clonable-and-runnable as the README
claims.

## 7. Processed data -- GEO deposit, and stepTwo readiness

Checked 2026-10-28 by resolving every path the stepTwo wrappers build from
their manifests and testing it on disk. A `raw_feature_bc_matrix` counts as
present only if it actually contains a `matrix.mtx*`; several exist as empty
directories.

| Arm | Libraries | BAM | Matrix |
|---|---:|---:|---:|
| Figure 1 Ca\*Re (Mmul_10_mac239) | 30 | 30 | 30 |
| Figure 1 TenX | 6 | 6 | **0** |
| Figure 3 Ca\*Re (GRCh38_HXB2) | 15 | 15 | 15 |
| Figure 3 TenX | 3 | 3 | 3 |
| Figure 4 Ca\*Re (SIV in vitro) | 3 | 3 | 3 |
| Figure 4 TenX (GFP_invitro) | 1 | 1 | 1 |
| Figure 4 liver Ca\*Re | 12 | **0** | **0** |
| Figure 4 liver TenX | 2 | **0** | **0** |
| **Total (ours)** | **72** | **58** | **52** |

The matrices are ~0.5 GB in total; the BAMs are ~90 GB. GEO takes processed
data, so the deposit is the matrices -- `barcodes`, `features` and `matrix.mtx`
per library. The BAMs are not normally deposited: the reads go to SRA (section
1) and the BAM is reproducible from them against the reference the methods
describe.

`extractedData/Figure3/CaRe/counts` is a symlink to
`extractedData/EGS026/SALVE/counts` (20 GB). Follow it when packaging.

Public-data counts (Collora 18, HIVseq 24, Giron 4) are re-analyses of
already-deposited data. Cite the accessions in section 3; deposit nothing.

### Can the pipeline start at stepTwo?

**Yes for every arm except Figure 4 liver.** stepTwo needs only the cellranger
BAM -- the stepTwo workers never read `RAW_CELLIDS_*`, which is a stepThree
input (`bamsort_process.R`). 103 of the 115 BAMs the manifests name are
present; the 12 that are missing are all liver Ca\*Re:

    extractedData/Figure4/liver/CaRe/counts/nointron/GRCm39_EGL002_{WT,KO}_{Kras,Nf2,Slkv3,Slkv4,Lsm14b,Dgkd}/outs/possorted_genome_bam.bam

Their run directories exist but `outs/` holds only the two matrix folders, and
those folders are empty. The liver TenX runs are in the same state. With
EGL001/EGL002 FASTQs also gone (section 2), **the liver arm cannot be rerun
from any step**, and has nothing to deposit.

stepThree then succeeds for Figures 1, 3 and 4 even though many
`raw_feature_bc_matrix` folders are empty, because `bamsort_process.R` prefers
a cached `<sample>_raw_cellIDs.csv` and those caches exist:

| Cache | Samples |
|---|---|
| `Figure1/singleCell/raw_cellIDs/` | D0, D195, D83, Invitro, Pacute, W0, W2 |
| `Figure3/singleCell/raw_cellIDs/` | HIV-A, HIV-C, HIV-D (and A_, C_, D_) |
| `Figure4/singleCell/raw_cellIDs/` | GFP_invitro, GEX_GFP_invitro |
| `Figure4/liver/singleCell/raw_cellIDs/` | **empty** |

Those caches are load-bearing. Delete them and stepThree fails for the six
Figure 1 TenX libraries and the Figure 4 in vitro library, whose matrices are
gone -- so they belong in the deposit, or at least in the backup, alongside the
matrices.

## Summary for the deposit

| Category | Volume | Action |
|---|---|---|
| Our sequencing, in hand | ≈646G, 6 directories | upload to SRA |
| Our sequencing, FASTQs lost | 6 accessions | recover first, or declare unavailable |
| Our imaging | 469G | image repository, not SRA |
| Published data | 99 runs | cite accessions |
| Our processed matrices | 52 of 72 libraries, ~0.5 GB | deposit to GEO |
| Libraries with no matrix | 20 (6 Figure 1 TenX, 14 liver) | nothing to deposit |
| References | 70G | not deposited; methods must pin the mac239 contig |
