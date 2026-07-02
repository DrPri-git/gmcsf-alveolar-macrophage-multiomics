# GM-CSF vs M-CSF Human Macrophage Transcriptome Reanalysis

**Part of a self-directed multi-omics learning project**, undertaken to build hands-on transcriptomic and multi-omics analysis skills relevant to alveolar macrophage immunobiology, innate immunity, and host-pathogen response — the core themes of my PhD research (Immunology, Newcastle University).

## Background

GM-CSF and M-CSF drive distinct macrophage differentiation programs. GM-CSF-conditioned macrophages (GM-MØ) display proinflammatory, immunogenic, and antimicrobial properties, and this differentiation pathway is essential for alveolar macrophage maturation and lung innate host defense (Shibata et al., 2001; Trapnell & Whitsett, 2002). In contrast, M-CSF-conditioned macrophages (M-MØ) are associated with an anti-inflammatory, tissue-remodeling phenotype. GM-CSF additionally links innate and adaptive immunity in the lung via PU.1-dependent regulation of Fc-receptor-mediated phagocytosis (Berclaz et al., 2002). Nieto et al. (2018) profiled the transcriptomes of human GM-CSF- and M-CSF-conditioned monocyte-derived macrophages and identified an activin A–PPARγ axis shaping the GM-MØ gene signature; their microarray data are publicly deposited as **GEO GSE68061**.

## Research Question

Which transcriptional programs distinguish GM-CSF-conditioned from M-CSF-conditioned human macrophages, and do they reproduce the macrophage polarization markers reported in the literature?

## Aim

To independently reanalyze GSE68061 and characterize the differential gene expression signature distinguishing GM-CSF- and M-CSF-conditioned human monocyte-derived macrophages.

## Objectives

1. Retrieve and quality-check the GSE68061 dataset and sample metadata.
2. Fit a linear model controlling for monocyte subset (CD14+CD16⁻ vs CD16+) to isolate the GM-CSF vs M-CSF treatment effect.
3. Identify and annotate differentially expressed genes.
4. Visualize and interpret the result against the published macrophage polarization literature.

## Data Source

| Field | Detail |
|---|---|
| Accession | GEO GSE68061 |
| Platform | Agilent GPL6480 (4x44K one-color microarray) |
| Organism | *Homo sapiens* |
| Design | 2 (CD14+CD16⁻ / CD16+ monocyte-derived macrophages) × 2 (GM-CSF 1000 U/mL / M-CSF 10 ng/mL) × 3 donor replicates = 12 samples |
| Source publication | Nieto et al., 2018, *Front Immunol* |

## Methods

Data were retrieved in R (v4.5.2) using `GEOquery::getGEO()`. Sample metadata (`treatment:ch1`, `cell type:ch1`) were extracted from `pData()` and verified for balanced group sizes (3 replicates per condition; confirmed by cross-tabulation). A linear model was fit using `limma::lmFit()` with the design `~ treatment + celltype`, controlling for monocyte subset so that the treatment effect estimate is not confounded by CD14/CD16 origin. Treatment was releveled with M-CSF as reference so that positive log-fold-change values indicate higher expression in GM-CSF-conditioned macrophages. Moderated t-statistics were computed with `eBayes()`, and results extracted with `topTable()`. Probes were annotated to gene symbols using the GPL6480 platform feature data (`fData()`); probes without a gene symbol were removed, and where multiple probes mapped to the same gene, the probe with the lowest p-value was retained. This yielded 19,595 uniquely annotated genes from an initial 41,076 probes. Differential expression was visualized using `EnhancedVolcano` (log2 fold-change vs. −log10 adjusted p-value) and `pheatmap` (row-scaled expression of the top 30 genes by significance, with unsupervised hierarchical clustering of samples and genes).

## Results

Differential expression analysis identified genes significantly associated with GM-CSF vs. M-CSF conditioning after multiple-testing correction (Benjamini-Hochberg adjusted p < 0.05). Among the most significant genes, **MAFB**, **THBD**, **CD93**, **MPEG1**, and **RNASE1** were more highly expressed in M-CSF-conditioned macrophages — consistent with their established roles as M-CSF-driven, anti-inflammatory/tissue-resident macrophage markers. **CXCL10**, typically associated with proinflammatory macrophage activation, was unexpectedly higher in M-CSF-conditioned cells in this dataset rather than GM-CSF, which runs counter to its canonical M1-associated role; this is noted as an area for further investigation rather than forced into the expected polarization framework. Unsupervised hierarchical clustering of the top 30 differentially expressed genes separated samples primarily by treatment condition (GM-CSF vs. M-CSF) rather than by monocyte subset of origin (Figure 2), indicating that GM-CSF/M-CSF conditioning — not CD14/CD16 monocyte subset — is the dominant driver of transcriptional variation in this dataset.

**Figure 1.** Volcano plot of differential gene expression, GM-CSF- vs. M-CSF-conditioned human macrophages (GSE68061 reanalysis). Dashed lines indicate significance (adj. p < 0.05) and fold-change (|log2FC| > 1) thresholds. Selected polarization-relevant genes are labeled.
`volcano_GSE68061.png`

**Figure 2.** Heatmap of the top 30 differentially expressed genes (row-scaled expression), with unsupervised hierarchical clustering of samples and genes. Column annotations indicate treatment (GM-CSF/M-CSF) and monocyte subset (CD14+CD16⁻/CD16+).
`heatmap_GSE68061.png`

## Interpretation

This reanalysis independently recovers established GM-CSF/M-CSF macrophage polarization biology from public transcriptomic data, validating the analysis pipeline against the published literature (Nieto et al., 2018; Shibata et al., 2001; Berclaz et al., 2002). This dataset characterizes the "priming" state of human macrophages and forms the first component of a planned multi-omics comparison with a second dataset (GSE213547 / PXD035269; Zec et al., 2023) examining alveolar macrophage transcriptomic and proteomic responses during *Streptococcus pneumoniae* infection, with the aim of relating GM-CSF-dependent macrophage maturation to downstream innate immune responses relevant to neutrophil recruitment.

## Limitations

This is a self-directed reanalysis of publicly available data using standard, well-established methods (limma), undertaken for skill-building purposes. It does not represent novel wet-lab data generation. The choice to retain only the lowest-p-value probe per gene during deduplication is one reasonable convention among several possible approaches. Given the modest sample size (n=3/group), results should be interpreted as hypothesis-generating rather than definitive.

## Software

R 4.5.2; Bioconductor packages `GEOquery`, `limma`, `Biobase`; `EnhancedVolcano`; `pheatmap`.

## References

1. Shibata Y, Berclaz PY, Chroneos ZC, Yoshida M, Whitsett JA, Trapnell BC. GM-CSF regulates alveolar macrophage differentiation and innate immunity in the lung through PU.1. *Immunity*. 2001;15(4):557-567.
2. Berclaz PY, Shibata Y, Whitsett JA, Trapnell BC. GM-CSF, via PU.1, regulates alveolar macrophage FcγR-mediated phagocytosis and the IL-18/IFN-γ-mediated molecular connection between innate and adaptive immunity in the lung. *Blood*. 2002;100(12):4193-4200.
3. Trapnell BC, Whitsett JA. GM-CSF regulates pulmonary surfactant homeostasis and alveolar macrophage-mediated innate host defense. *Annu Rev Physiol*. 2002;64:775-802.
4. Nieto C, Bragado R, Municio C, Sierra-Filardi E, Alonso B, Escribese MM, Domínguez-Andrés J, Ardavín C, Castrillo A, Vega MA, Puig-Kröger A, Corbí AL. The Activin A-Peroxisome Proliferator-Activated Receptor Gamma Axis Contributes to the Transcriptome of GM-CSF-Conditioned Human Macrophages. *Front Immunol*. 2018;9:31. (Data: GEO GSE68061)
5. Zec K, Thiebes S, Bottek J, Siemes D, Spangenberg P, Trieu DV, Kirstein N, Subramaniam N, Christ R, Klein D, Jendrossek V, Loose M, Wagenlehner F, Jablonska J, Bracht T, Sitek B, Budeus B, Klein-Hitpass L, Theegarten D, Shevchuk O, Engel DR. Comparative transcriptomic and proteomic signature of lung alveolar macrophages reveals the integrin CD11b as a regulatory hub during pneumococcal pneumonia infection. *Front Immunol*. 2023;14:1227191. (Data: GEO GSE213547; PRIDE PXD035269)
6. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47.
