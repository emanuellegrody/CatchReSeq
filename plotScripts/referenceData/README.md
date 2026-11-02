# referenceData

Inputs SuppFigs.Rmd reads directly, kept here so the panels run from a clone.

## Amino-acid diversity (SuppFig1 Panel E, SuppFig3 Panel F)

`GeneticDiversity_AminoAcids_SIVmac239.txt` and
`GeneticDiversity_AminoAcids_HIV1cladeB.txt` are per-position amino-acid
diversity tables, with the alignments they were computed from
(`Mac239Genome_AminoAcidAlignment.fasta`, `cladeB_Genome_AminoAcidAlignment.fasta`).

They are outputs of the HIVToolkit MATLAB pipeline, which is third-party and not
part of this repository: the SIVmac239 table came from
`extractionScripts/HIVToolkit/` and the clade B table from that toolkit's
`EGS_026/` run. Regenerating them needs the toolkit; the four files here are
what the figures consume.
