# PEL-seq

This repository contains the scripts used to generate the main and supplementary figure analyses for the PEL-seq study.

PEL-seq (Probe-Enriched LINE-1 Sequencing) is a targeted Nanopore-based approach for sensitive, locus-resolved profiling of full-length LINE-1 transcripts. The analyses in this repository examine locus-specific LINE-1 transcription across cell types, its epigenetic regulation, remodeling during acquisition of osimertinib resistance, and CRISPRa-mediated LINE-1 reactivation.

## Scripts

| Script | Description |

| `Fig1_S1_PEL_seq_performance.R` | PEL-seq performance, read characteristics, LINE-1 detection, subfamily expression, and contribution of highly expressed loci. |

| `Fig2_S2_cell_type_specificity.R` | Cell-type-specific LINE-1 transcription and associations between LINE-1 loci and neighboring gene expression. |

| `Fig3_S3_epigenetic_regulation.R` | Integration of PEL-seq with H3K27ac, H3K9me3, ATAC-seq, KAS-seq, and DNA methylation data. |

| `Fig4_S4_S5_osimertinib_resistance.R` | Locus-specific LINE-1 transcription and epigenetic remodeling during acquisition of osimertinib resistance. |

| `Fig5_S6_CRISPRa_reactivation.R` | CRISPRa-mediated LINE-1 reactivation, neighboring gene expression, interferon-response signatures, and osimertinib sensitivity. |

## Data availability

Raw sequencing data generated in this study are deposited in the public repositories described in the manuscript.

Large sequencing files, including FASTQ, BAM, and bigWig files, are not included in this GitHub repository. The scripts use processed count matrices, genomic annotations, differential-expression results, and signal matrices derived from PEL-seq, RNA-seq, CUT&Tag, KAS-seq, ATAC-seq, and WGBS data.



