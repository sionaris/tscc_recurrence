<h1 align="center">Transcriptomic signatures of recurrence in oral squamous cell carcinoma</h1>

<p align="center">
  <img alt="R" src="https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white">
  <img alt="Bioconductor" src="https://img.shields.io/badge/Bioconductor-DESeq2%20%C2%B7%20IHW%20%C2%B7%20clusterProfiler-1a81c2">
  <img alt="nf-core/rnaseq" src="https://img.shields.io/badge/nf--core-rnaseq%203.14-24B064">
  <img alt="Data" src="https://img.shields.io/badge/data-bulk%20RNA--seq%20(single--end)-lightgrey">
</p>

---

## Scope

This repository contains the analysis code and outputs for a bulk RNA-seq study of **oral squamous cell carcinoma (OSCC)**. It asks which genes and pathways distinguish tumours that went on to **recur** from those that did not, and whether that signal depends on **tumour stage** (pathological T stage).

| | |
|---|---|
| **Cohort** | 40 samples (32 tumours and 8 adjacent normal tissues) from 34 patients. Each sample was sequenced on two lanes (read 1 only), and the lane counts were summed per sample |
| **Analysed set** | 32 tumours: early stage (T1–T2) and late stage (T3–T4), with and without recurrence |
| **Model** | `DESeq2` Wald test, `~ Stage * Recurrence`, with IHW p-value weighting and `ashr` LFC shrinkage |
| **Downstream** | DEG overlap statistics, ORA, GSEA and active-subnetwork enrichment (`pathfindR`) |

<table>
  <tr><th></th><th>Early (T1–T2)</th><th>Late (T3–T4)</th></tr>
  <tr><td><b>Recurrence</b></td><td align="center">7</td><td align="center">8</td></tr>
  <tr><td><b>No recurrence</b></td><td align="center">9</td><td align="center">8</td></tr>
</table>

### Contrasts

| Code | Comparison |
|---|---|
| `TR_vs_TNR` | Recurrent vs non-recurrent tumours (stage-balanced) |
| `ETR_vs_ETNR` | Early-stage: recurrent vs non-recurrent |
| `LTR_vs_LTNR` | Late-stage: recurrent vs non-recurrent |
| `LTR_vs_ETR` | Recurrent tumours: late vs early stage |

---

## Workflow

```mermaid
flowchart LR
    A["FASTQ<br/>(40 samples × 2 lanes)"] --> B["nf-core/rnaseq<br/>STAR + Salmon, GRCh38"]
    B --> C["1 · Preprocessing &<br/>exploratory analysis"]
    K["clinchar.csv"] --> D["1b · Clinical<br/>characteristics (Table 1)"]
    C --> E["2 · Differential<br/>expression (DESeq2)"]
    E --> F["2b · DEG overlaps<br/>across contrasts"]
    E --> G["3 · Pathway analysis<br/>ORA · GSEA · pathfindR"]
    C & E & G --> H["HTML report"]
```

| Step | Script | What it does |
|---|---|---|
| 0 | [`nextflow/nf_job.sh`](nextflow/nf_job.sh) | SLURM array job running `nf-core/rnaseq` per sample (Ensembl GRCh38.112) |
| 1 | [`1. Preprocessing_and_exploratory_analysis.R`](R_Scripts/1.%20Preprocessing_and_exploratory_analysis.R) | `tximport` of Salmon quantifications, filtering, summing the two lanes of each sample, normalisation, VST/rlog, distance heatmaps, PCA / GLM-PCA / MDS |
| 1b | [`1.b. clinchar_stats.R`](R_Scripts/1.b.%20clinchar_stats.R) | Group comparisons of clinicopathological variables (Fisher, Kruskal–Wallis, Wilcoxon; BH-adjusted), exported as Table 1 |
| 2 | [`2. DGEA.R`](R_Scripts/2.%20DGEA.R) | `DESeq2` interaction model on tumours, IHW diagnostics, `ashr` shrinkage, MA / volcano plots, DEG heatmaps |
| 2b | [`2.b. DGEA_overlaps.R`](R_Scripts/2.b.%20DGEA_overlaps.R) | Pairwise DEG overlap between contrasts: Fisher tests, Jaccard / signed Jaccard, sign concordance, recovery; multi-panel figure |
| 3 | [`3. Pathway_analysis.R`](R_Scripts/3.%20Pathway_analysis.R) | ORA and GSEA (`clusterProfiler`, `ReactomePA`) plus `pathfindR` active-subnetwork analysis across KEGG, Reactome, BioCarta and GO |
| — | [`fast_pathfindR_hclust.R`](R_Scripts/fast_pathfindR_hclust.R) | Faster drop-in kappa-matrix / hierarchical clustering for `pathfindR` term clustering (via `fastcluster`) |

---

## Repository structure

```text
oscc_recurrence/
├── nextflow/                     # HPC job script for nf-core/rnaseq
├── R_Scripts/                    # Analysis scripts, numbered in run order
├── Resources_and_environments/   # Inputs and saved R objects
│   ├── samplesheet.csv           #   nf-core sample sheet
│   ├── Pheno.xlsx                #   sample phenotype data
│   ├── summed_counts.xlsx        #   per-sample count matrix (lanes summed)
│   ├── final_dds.rds             #   final DESeq2 object
│   └── *_sessionInfo.txt         #   package versions for each step
├── Results/
│   ├── Exploratory analysis/     # QC, density, distance and dim-reduction plots
│   ├── DGEA/
│   │   ├── DGEA_raw.xlsx         #   DE results, unshrunk LFCs (one sheet per contrast)
│   │   ├── DGEA_ashr.xlsx        #   DE results, ashr-shrunk LFCs
│   │   ├── IHW_diagnostics/
│   │   ├── MA_plots/  Volcano_plots/  Post_DGEA_heatmaps/
│   │   └── Overlaps/             #   cross-contrast overlap matrices and figure
│   └── Pathways/
│       ├── ashr_<contrast>/      #   ORA/, GSEA/, pathfindR/<database>/
│       └── archived/             #   earlier stage-only contrasts
├── Reports/                      # R Markdown report (+ rendered HTML) and its figures
├── Tables/                       # Table 1: clinicopathological characteristics (.csv, .docx)
├── clinchar.csv                  # Clinical characteristics per sample
└── oscc_seq_rec.Rproj            # RStudio project (run scripts from the repo root)
```

---

## Reproducing the analysis

1. **Alignment and quantification.** Run [`nextflow/nf_job.sh`](nextflow/nf_job.sh) on an HPC cluster with the sample sheet in `Resources_and_environments/`. Outputs go to `NF-processed_samples/batch_*`, which is git-ignored because of its size. FASTQ files can be found in SRA: [PRJNA###](#LINK). Summarized counts available in `Resources_and_environments/summed_counts.xlsx`.
2. **R analysis.** Open `oscc_seq_rec.Rproj` and run the scripts in `R_Scripts/` in numerical order. All paths are relative to the repository root.

Package versions for each stage are recorded in `Resources_and_environments/*_sessionInfo.txt`.

> **Note:** Raw sequencing data and intermediate `.RData` environments are not tracked in this repository.
