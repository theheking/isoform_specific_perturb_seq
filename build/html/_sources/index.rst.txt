Isoform-specific Perturb-seq
============================

.. attention::
   ***Under construction***

Isoform-specific single-cell Perturb-seq reveals that alternative promoters (AP) are not just redundant regulatory elements but drive distinct biological programs. This documentation covers the methods used to identify these promoters and quantify their functional impact on gene regulation and drug response.


Abstract and Key Findings
-------------------------
CRISPR-dCas9 technologies are typically designed to modulate gene expression at the gene level. However, many hits are promoter and isoform-specific. By leveraging the spatial specificity of CRISPRi—which typically influences transcription within ~1000 nucleotides of the guide binding site—we developed an isoform-specific Perturb-Seq screen[cite: 66, 68].

**Major Biological Insights:**
* **Widespread Functional Divergence:** Alternative promoters drive functionally distinct biological programs in **51.6%** of surveyed genes.
* **Limited Transcriptional Overlap:** Promoter-specific knockdowns (P1 vs P2) typically perturb dozens to hundreds of genes, but only a small minority (median ~10) overlap between the two promoters of the same gene[cite: 636].
* **Cell Cycle Regulation:** Alternative promoters frequently control divergent pathways in cell cycle regulation and proliferation[cite: 177, 180].

Functional Case Study: ESR1 and Tamoxifen Response
--------------------------------------------------
The Estrogen Receptor 1 (*ESR1*) serves as a primary example of how alternative promoters influence clinical outcomes and drug sensitivity in breast cancer.


* **Isoform Structure:** Coordinated splicing between the P1 promoter and an alternative last exon produces a protein isoform with differences in the **AF2 domain**, potentially modifying interactions with estrogen and selective estrogen receptor modulators (SERMs)[cite: 185, 186].
* **Clinical Significance:** High expression of the P1 promoter strongly correlates with decreased survival in **Luminal-A** breast cancer patients (HR = 1.9), while P2 expression shows no such association[cite: 189].
* **Drug Response:**
    * **P2 Knockdown:** Increases sensitivity to tamoxifen and significantly reduces cellular proliferation[cite: 184, 197].
    * **P1 Knockdown:** Leads to increased proliferation in the presence of tamoxifen, suggesting a role in drug resistance[cite: 198].

Analysis Pipeline Overview
--------------------------

The analysis is structured into three main phases:

1. **Promoter Identification**
   Integration of RNA-seq, ChIP-seq, and CAGE-seq to identify targetable distal alternative promoters missed by standard CRISPRi libraries[cite: 144, 145].

2. **Guide Design**
   Utilizing FlashFry for promoter-specific dual-guide design to ensure spatial targeting within the window of CRISPRi effectiveness.

3. **APU Perturb-seq Analysis**
   Quantifying functional divergence through:
   * **Transcriptomics:** Using **Whippet** for isoform-specific quantification post-UMI deduplication[cite: 126].
   * **Pathway Activity:** Using **Spectra** to identify coordinated gene expression programs (e.g., cell cycle phases) associated with specific promoters[cite: 174, 175].
   * **Chromosomal Dynamics:** Using **inferCNV** to identify increases in copy number variations among cell cycle-related genes specifically in P2 knockdown populations[cite: 183].

.. toctree::
   :maxdepth: 2
   :caption: Analysis Pipeline:

   promoter_identification
   guide_design
   analysis_scripts