# Python environments

Two environments are enough to run everything: one for smFISH and one for the
rest. The published analysis ran in four, and those four are recorded here as
well, because what was tested and what is recommended are not the same thing.

## Recommended

| environment | python | manifest | covers |
|---|---|---|---|
| `care` | 3.13.11 | [`care.txt`](care.txt) | read extraction, splice, polyA, cNMF, BRIE2 |
| `bigfish_env` | 3.6.15 | [`bigfish.txt`](bigfish.txt) | smFISH spot detection and segmentation |

```sh
conda create -n care python=3.13.11 && conda activate care
pip install -r environments/care.txt
```

The workers activate three different names -- `cellranger`, `jupyter` and
`brie2` -- so either create `care.txt` under each of those three names, or
create it once and change the `source activate` line in the workers that
reference it. The table below says which worker uses which name.

smFISH stays separate and cannot be folded in. big-fish 0.6.2 calls numpy
aliases (`np.int`, `np.float`, `np.bool`) removed in numpy 1.24, and imports
`skimage.morphology.selem`, removed in scikit-image 0.19; it needs the last
stack in which those still exist. Only the spot-detection and segmentation
scripts need it. `hcr_positivity_overlay.py` and `hcr_positivity_panel.py` do
not import big-fish and run in either environment.

One pin in `care.txt` is carried over rather than observed: gffutils, which
only `polya/polya_utils.py` imports, ran as 0.12 under python 3.9 and is listed
as 0.14 here, the first release supporting python 3.13.

## What the published analysis ran on

| environment | python | manifest | activated in |
|---|---|---|---|
| `cellranger` | 3.9.18 | [`bamsort.txt`](bamsort.txt) | the bamsort and cellranger workers, and `Figure2/{bamsort_splice,bamsort_isoform_longread,run_polya}.sh` |
| `jupyter` | 3.13.11 | [`cnmf.txt`](cnmf.txt) | `Figure3/cNMF_*.sh` |
| `brie2` | 3.11.14 | [`brie2.txt`](brie2.txt) | `isoform/brie2_worker.sh` |
| `bigfish_env` | 3.6.15 | [`bigfish.txt`](bigfish.txt) | `Figure2/bigfish_*.sh` and `hcr_positivity_plots.sh` |

To build one, substitute the environment name and the python version named at
the top of its manifest:

```sh
conda create -n cellranger python=3.9.18 && conda activate cellranger
pip install -r environments/bamsort.txt
```

Merging the first three changes the versions the extraction arm sees, most
consequentially pysam 0.19.0 to 0.23.3 and numpy 1.24.4 to 2.3.5. Every import
the pipeline makes has been checked under the merged environment; the published
numbers were produced under the four.

## Nextflow

The bulk RNA-seq arm uses none of the above. `nfcore_rnaseq_worker.sh`
activates a separate environment, named `nextflow` by default and overridable
with `NEXTFLOW_ENV`, holding only Nextflow 25.10.4 and OpenJDK 21.0.10;
nf-core supplies the rest through Singularity containers. The Figure 2
long-read arm uses the same environment, through
`submissionScripts/Figure2/longread_scnanoseq.sh`.
