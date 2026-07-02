# GM-CSF Priming and Pneumococcal Infection Response: A Multi-Dataset Transcriptomic Reanalysis

**Self-directed multi-omics learning project**, undertaken to build hands-on transcriptomic analysis and cross-dataset comparison skills relevant to alveolar macrophage immunobiology, innate immunity, and host-pathogen response.

---

## Introduction

Alveolar macrophages are shaped by two distinct biological processes relevant to lung host defense: their developmental priming by colony-stimulating factors, and their acute functional response upon pathogen encounter. This project independently reanalyzes two publicly available datasets addressing each side of this axis, then asks whether genes associated with macrophage priming reappear in the acute infection response — a comparison that, to my knowledge, has not been performed between these two specific datasets.

## Background

GM-CSF drives terminal alveolar macrophage differentiation and is required for effective innate lung host defense (Shibata et al., 2001; Trapnell & Whitsett, 2002), and additionally bridges innate and adaptive immunity via PU.1-dependent regulation of Fc-receptor phagocytosis (Berclaz et al., 2002). Nieto et al. (2018) profiled human GM-CSF- and M-CSF-conditioned macrophage transcriptomes and identified an activin A–PPARγ axis underlying the GM-CSF-specific gene signature (GEO: GSE68061).

Separately, Zec et al. (2023) used comparative transcriptomic and proteomic profiling to show that murine alveolar macrophages intrinsically upregulate CD11b (Itgam) during *Streptococcus pneumoniae* lung infection, identifying CD11b as a regulatory hub influencing neutrophil recruitment, activation, and migration (GEO: GSE213547; PRIDE: PXD035269).

These two studies characterize different biological states of the same cell type — macrophage priming versus macrophage response to infection — but have not previously been analyzed together.

## Research Question

Do genes associated with GM-CSF-driven macrophage priming reappear among genes differentially expressed during alveolar macrophage response to *Streptococcus pneumoniae* infection?

## Aim

To independently reanalyze GSE68061 and GSE213547, and to test whether the resulting gene signatures show greater overlap than expected by chance.

## Objectives

1. Reanalyze GSE68061 (GM-CSF vs. M-CSF human macrophages) — completed (see Part 1).
2. Reanalyze GSE213547 (*S. pneumoniae*-infected vs. mock-infected murine alveolar macrophages).
3. Compare the two resulting differential expression gene lists using a formal statistical enrichment test.
4. Interpret any overlap — or lack thereof — in the context of the source literature.

---

## Part 1: GM-CSF vs. M-CSF Human Macrophages (GSE68061)

See separate README (`README_GSE68061.md`) for full methods, results, and figures for this dataset. In summary: differential expression analysis (limma, controlling for monocyte subset) identified 1,156 genes significant at adjusted p < 0.05, with established M-CSF/M2-associated markers (MAFB, THBD, CD93, MPEG1, RNASE1) among the top hits, consistent with Nieto et al. (2018).

## Part 2: *S. pneumoniae* Infection vs. Mock (GSE213547)

### Data Source

| Field | Detail |
|---|---|
| Accession | GEO GSE213547 |
| Platform | Affymetrix Clariom S Mouse (GPL23038) |
| Organism | *Mus musculus* |
| Design | 5 infected (*S. pneumoniae*) vs. 5 mock-infected lung macrophage samples |
| Source publication | Zec et al., 2023, *Frontiers in Immunology* |

### Methods

Data were retrieved via `GEOquery::getGEO()`. Metadata confirmed a simple two-group design (genotype labels were redundant with treatment and were not modeled separately). Expression values were already log2-transformed (range 3.6–13.8), consistent with the deposited GC-RMA/vsn-normalized data. A linear model (`~ treatment`) was fit using `limma`, with mock infection as the reference level. Because this Affymetrix Clariom platform's annotation lacks a clean gene-symbol field, gene symbols were extracted programmatically from the first entry of the `mrna_assignment` text field (pattern: symbol given in parentheses immediately preceding a comma); probes without a resolvable symbol (1,091 of 22,206, 4.9%) were excluded. Where multiple probes mapped to the same gene, the probe with the lowest p-value was retained (154 duplicates resolved), yielding 20,961 uniquely annotated genes.

### Results

Infection produced a substantially larger transcriptional response than the GM-CSF/M-CSF comparison (5,135 genes significant at adjusted p < 0.05, versus 1,156 in GSE68061), with markedly larger effect sizes (top t-statistics exceeding 36, versus a maximum of ~14 in GSE68061) — consistent with an acute bacterial infection producing a broader, more synchronized transcriptional shift than a difference in culture growth factor. Top upregulated genes included classic acute innate immune response genes: **Il1a**, **Tnfaip3**, **Nfkb2** (NF-κB pathway activation), **Irg1** (itaconate production), **Clec4e**/Mincle (pattern recognition), **Lcn2** (antimicrobial protein), and neutrophil-recruiting chemokines **Cxcl1**, **Ccl9**, and **Ccl4**. **Itgam** (CD11b) was significantly upregulated during infection (log2FC = 1.30, adjusted p = 8.8×10⁻⁴), directly reproducing the central finding of Zec et al. (2023) in an independent reanalysis. **Inhba** (Activin A) — the central molecule of the Nieto et al. (2018) GM-CSF study — also emerged among the most significant infection-response genes (adjusted p = 1.3×10⁻¹⁰), an initial point of connection between the two datasets explored formally below.

**Figure 3.** Volcano plot, *S. pneumoniae* infection vs. mock (GSE213547 reanalysis). `volcano_GSE213547.png`

**Figure 4.** Heatmap of top 30 differentially expressed genes; infected and mock samples formed two complete, non-overlapping clusters. `heatmap_GSE213547.png`

---

## Part 3: Cross-Dataset Comparison

### Methods

Gene symbols from both datasets were case-standardized (converted to uppercase) to allow direct symbol matching between human (GSE68061) and mouse (GSE213547) gene sets. This is a simplified approach relying on the convention that conserved orthologous genes generally share symbol names across species; it does not account for orthologs with divergent nomenclature and does not constitute a validated orthology mapping (e.g., via Ensembl biomaRt), which would be a more rigorous approach for a definitive analysis. Genes significant (adjusted p < 0.05) in both datasets were identified, restricted to the background of genes measured on both platforms (n = 14,478). Statistical enrichment of the overlap was tested using a one-sided hypergeometric test.

### Results

Of 14,478 genes common to both platforms, 1,156 were significant in the GM-CSF dataset and 5,135 in the infection dataset. The observed overlap was 362 genes, closely matching the 410 genes expected under random chance given the size of both lists, and the hypergeometric test found no statistical enrichment (p = 0.999). This indicates that, at the level of bulk gene-list overlap, GM-CSF priming and pneumococcal infection response do not share more transcriptional signal than expected by chance — a result likely driven in part by the very large size and breadth of the infection-response gene list, which alone limits the discriminating power of a simple overlap test.

Despite the absence of bulk statistical enrichment, several individual genes of direct mechanistic relevance to the GM-CSF/macrophage-priming literature were independently identified as significant in the infection dataset, including **INHBA** (Activin A, the central molecule of the Nieto et al. 2018 GM-CSF mechanism), **CSF1** and **CSF1R** (the M-CSF ligand-receptor pair), **PPARG** (the transcription factor central to the GM-CSF/activin A axis), **CEBPA/CEBPD** (macrophage-lineage transcription factors), pattern-recognition receptors **TLR2/TLR4/TLR5** and **NLRP3**, and **CD14**. These represent targeted, mechanistically motivated points of overlap rather than evidence of a broadly shared transcriptional program, and are presented as hypothesis-generating observations warranting focused follow-up rather than as a proven mechanistic link.

## Summary

This project independently reanalyzed two public transcriptomic datasets addressing complementary aspects of alveolar macrophage biology — GM-CSF-driven priming and acute response to *Streptococcus pneumoniae* infection — recovering the principal published findings of both source studies (Nieto et al., 2018; Zec et al., 2023). A formal statistical test of gene-list overlap between the two datasets found no significant bulk enrichment, a result reported transparently rather than reinterpreted; however, several specific genes central to the GM-CSF/activin A/PPARγ priming mechanism were independently found to be significant in the infection response, representing plausible, targeted hypotheses for future investigation into how macrophage priming state may shape the innate response to pneumococcal infection.

## Limitations

- Cross-species gene matching used direct symbol comparison rather than a validated orthology database; this may both miss true orthologs with divergent symbols and is not the gold-standard approach for definitive cross-species comparison.
- The two datasets differ in species (human vs. mouse), platform (Agilent vs. Affymetrix), and biological context (in vitro culture vs. in vivo infection), all of which limit direct comparability beyond the gene-symbol level.
- The infection dataset's very large significant gene list (5,135 genes) inherently limits the power of a simple overlap-based enrichment test to detect meaningful specificity; a gene-set-level approach (e.g., GSEA) applied to curated pathway sets may be more informative than whole-list overlap for future work.
- As throughout this project, no novel wet-lab data were generated; this is a reanalysis of existing public data for self-directed skill-building.

## Software

R 4.5.2; Bioconductor packages `GEOquery`, `limma`, `Biobase`; `EnhancedVolcano`; `pheatmap`; base R `phyper()` for hypergeometric testing.

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
