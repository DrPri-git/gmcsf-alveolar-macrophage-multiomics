# GM-CSF vs M-CSF Human Macrophage Transcriptome Reanalysis

**Self-directed multi-omics learning project**, undertaken to build hands-on transcriptomic and multi-omics analysis skills relevant to alveolar macrophage immunobiology, innate immunity, and host-pathogen response — core themes of my PhD research (Immunology, Newcastle University).

---

## Introduction

Macrophages are functionally shaped by the colony-stimulating factor they differentiate under. Granulocyte-macrophage colony-stimulating factor (GM-CSF) and macrophage colony-stimulating factor (M-CSF) drive distinct, well-characterized macrophage differentiation programs, producing cells with different transcriptional profiles, surface markers, and functional roles in host defense versus tissue homeostasis. Understanding this axis is foundational to my PhD work on innate and adaptive immune responses in the lung, particularly in the context of alveolar macrophage biology and the host response to respiratory pathogens such as *Streptococcus pneumoniae*.

This project independently reanalyzes a publicly available human macrophage transcriptomic dataset to (1) build practical skill in transcriptomic data handling and differential expression analysis in R, and (2) generate a reference signature of GM-CSF- versus M-CSF-conditioned macrophages that can later be related to alveolar macrophage responses during infection, as part of a broader planned multi-omics learning exercise.

## Background

GM-CSF is essential for the terminal differentiation of alveolar macrophages and for innate lung host defense. Shibata et al. (2001) demonstrated that GM-CSF-deficient mice develop alveolar macrophages with defective phagocytosis, pathogen killing, and pattern-recognition receptor expression, mediated through the transcription factor PU.1. Subsequent work from the same group showed that GM-CSF, acting via PU.1, also regulates macrophage Fc-receptor-mediated phagocytosis and forms a molecular bridge between innate and adaptive immunity in the lung through an IL-18/IFN-γ-dependent mechanism (Berclaz et al., 2002), and that GM-CSF is required more broadly for pulmonary surfactant homeostasis and alveolar macrophage-mediated host defense (Trapnell & Whitsett, 2002).

*In vitro*, macrophages generated under GM-CSF (GM-MØ) display proinflammatory, immunogenic, and antimicrobial properties, in contrast to M-CSF-conditioned macrophages (M-MØ), which display a more anti-inflammatory, tissue-remodeling phenotype. Nieto et al. (2018) profiled the transcriptomes of human GM-MØ and M-MØ and identified an activin A–PPARγ regulatory axis underlying the GM-MØ-specific gene signature, with PPARγ expression and activity differing markedly between the two macrophage subtypes. Their microarray dataset is publicly deposited in the Gene Expression Omnibus as **GSE68061**, and forms the basis of this reanalysis.

## Research Question

Which transcriptional programs distinguish GM-CSF-conditioned from M-CSF-conditioned human macrophages, and does independent reanalysis of public data reproduce the macrophage polarization markers reported in the literature?

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
| Source publication | Nieto et al., 2018, *Frontiers in Immunology* |

## Methods

Data were retrieved in R (v4.5.2) using `GEOquery::getGEO()`. Sample metadata (`treatment:ch1`, `cell type:ch1`) were extracted from `pData()` and cross-tabulated to confirm a balanced design (3 replicates per condition, 12 samples total). A linear model was fit using `limma::lmFit()` with the design formula `~ treatment + celltype`, controlling for monocyte subset of origin so that the estimated treatment effect is not confounded by CD14/CD16 status. The treatment factor was releveled with M-CSF as the reference level, so that positive log2 fold-change values indicate higher expression in GM-CSF-conditioned macrophages. Moderated t-statistics and p-values were computed using `eBayes()`, and results were extracted with `topTable()`, ranked by p-value, with Benjamini-Hochberg adjustment for multiple testing across probes.

Probes were annotated to gene symbols using the GPL6480 platform feature data (`fData()`). Probes lacking a gene symbol (controls or unannotated features) were removed. Where multiple probes mapped to the same gene symbol, the probe with the lowest p-value was retained as the representative measurement for that gene. This filtering and deduplication reduced the dataset from 41,076 probes to 19,595 uniquely annotated genes.

Results were visualized as a volcano plot (log2 fold-change vs. −log10 adjusted p-value) using `EnhancedVolcano`, and as a heatmap of row-scaled expression for the top 30 genes by significance using `pheatmap`, with unsupervised hierarchical clustering applied to both samples and genes, and column annotations indicating treatment and monocyte subset.

## Results

Differential expression analysis identified genes significantly associated with GM-CSF versus M-CSF conditioning after multiple-testing correction (adjusted p < 0.05). Among the most statistically significant genes, **MAFB**, **THBD**, **CD93**, **MPEG1**, and **RNASE1** were more highly expressed in M-CSF-conditioned macrophages, consistent with their established roles as markers of M-CSF-driven, anti-inflammatory and tissue-resident macrophage phenotypes. **CXCL10**, a chemokine typically associated with proinflammatory (M1-like) macrophage activation, was unexpectedly found to be higher in M-CSF-conditioned cells in this dataset rather than GM-CSF; this departs from its canonical association and is noted here as a finding warranting further investigation rather than being forced to fit the expected polarization framework.

Unsupervised hierarchical clustering of the top 30 differentially expressed genes separated samples primarily by treatment condition (GM-CSF vs. M-CSF), with monocyte subset of origin (CD14+CD16⁻ vs. CD16+) intermixed within each treatment cluster rather than driving the separation. This indicates that GM-CSF/M-CSF conditioning, not monocyte subset, is the dominant source of transcriptional variation captured in this dataset — supporting the validity of the modeling approach used.

**Figure 1.** Volcano plot of differential gene expression, GM-CSF- vs. M-CSF-conditioned human macrophages (GSE68061 reanalysis). Dashed lines indicate significance (adjusted p < 0.05) and fold-change (|log2FC| > 1) thresholds. Selected polarization-relevant genes are labeled.
`volcano_GSE68061.png`

**Figure 2.** Heatmap of the top 30 differentially expressed genes (row-scaled expression), with unsupervised hierarchical clustering of samples and genes. Column annotations indicate treatment (GM-CSF/M-CSF) and monocyte subset (CD14+CD16⁻/CD16+).
`heatmap_GSE68061.png`

## Summary

This reanalysis independently recovers established GM-CSF/M-CSF macrophage polarization biology from a public transcriptomic dataset, validating the analysis pipeline against previously published findings (Nieto et al., 2018; Shibata et al., 2001; Berclaz et al., 2002). The resulting gene signature characterizes the "priming" state of human macrophages under GM-CSF versus M-CSF conditioning and forms the first component of a planned multi-omics comparison with a second, infection-focused dataset (GSE213547 / PXD035269; Zec et al., 2023), examining alveolar macrophage transcriptomic and proteomic responses during *Streptococcus pneumoniae* infection. The intended next step is to relate GM-CSF-dependent macrophage priming to downstream innate immune responses relevant to neutrophil recruitment during infection — directly connecting this exercise to the innate immunity, alveolar macrophage, and neutrophil-recruitment themes of my PhD research.

## Limitations

This is a self-directed reanalysis of publicly available data using standard, well-established statistical methods (limma), undertaken for skill-building purposes, and does not represent novel wet-lab data generation. Retaining the lowest-p-value probe per gene during deduplication is one reasonable convention among several possible approaches, and this choice is stated explicitly for transparency. Given the modest sample size (n=3 per group), findings should be interpreted as hypothesis-generating rather than definitive, and any discrepancies from the published literature (e.g., the CXCL10 direction) are reported openly rather than reinterpreted to fit expectations.

## Software

R 4.5.2; Bioconductor packages `GEOquery`, `limma`, `Biobase`; `EnhancedVolcano`; `pheatmap`.

---

## References

1. Shibata Y, Berclaz PY, Chroneos ZC, Yoshida M, Whitsett JA, Trapnell BC. GM-CSF regulates alveolar macrophage differentiation and innate immunity in the lung through PU.1. *Immunity*. 2001;15(4):557-567.
https://doi.org/10.1016/S1074-7613(01)00218-7

2. Berclaz PY, Shibata Y, Whitsett JA, Trapnell BC. GM-CSF, via PU.1, regulates alveolar macrophage FcγR-mediated phagocytosis and the IL-18/IFN-γ-mediated molecular connection between innate and adaptive immunity in the lung. *Blood*. 2002;100(12):4193-4200.
https://doi.org/10.1182/blood-2002-04-1102

3. Trapnell BC, Whitsett JA. GM-CSF regulates pulmonary surfactant homeostasis and alveolar macrophage-mediated innate host defense. *Annu Rev Physiol*. 2002;64:775-802.
https://doi.org/10.1146/annurev.physiol.64.090601.113847

4. Nieto C, Bragado R, Municio C, Sierra-Filardi E, Alonso B, Escribese MM, Domínguez-Andrés J, Ardavín C, Castrillo A, Vega MA, Puig-Kröger A, Corbí AL. The Activin A-Peroxisome Proliferator-Activated Receptor Gamma Axis Contributes to the Transcriptome of GM-CSF-Conditioned Human Macrophages. *Front Immunol*. 2018;9:31.
https://doi.org/10.3389/fimmu.2018.00031
Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68061

5. Zec K, Thiebes S, Bottek J, Siemes D, Spangenberg P, Trieu DV, Kirstein N, Subramaniam N, Christ R, Klein D, Jendrossek V, Loose M, Wagenlehner F, Jablonska J, Bracht T, Sitek B, Budeus B, Klein-Hitpass L, Theegarten D, Shevchuk O, Engel DR. Comparative transcriptomic and proteomic signature of lung alveolar macrophages reveals the integrin CD11b as a regulatory hub during pneumococcal pneumonia infection. *Front Immunol*. 2023;14:1227191.
https://doi.org/10.3389/fimmu.2023.1227191
Transcriptomic dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213547
Proteomic dataset: https://www.ebi.ac.uk/pride/archive/projects/PXD035269

6. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47.
https://doi.org/10.1093/nar/gkv007

7. Rohart F, Gautier B, Singh A, Lê Cao KA. mixOmics: An R package for 'omics feature selection and multiple data integration. *PLOS Comput Biol*. 2017;13(11):e1005752.
https://doi.org/10.1371/journal.pcbi.1005752
