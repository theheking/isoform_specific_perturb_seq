.. _analysis-scripts:

Analysis Scripts & Pipelines
============================

This section contains the core computational notebooks and scripts used for data processing, 
statistical analysis, and visualization.

.. toctree::
   :maxdepth: 1
   :caption: Analysis Scripts

   scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb
   scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides.ipynb
   scripts/3_analysis_scripts/12_negativecontrol.ipynb
   scripts/3_analysis_scripts/3C_whippet_pseudobulk
   scripts/3_analysis_scripts/4_estatistic_gene_guide.ipynb
   scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression.ipynb
   scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap.ipynb
   scripts/3_analysis_scripts/8_perturbseq_differentialexpression.ipynb
   scripts/3_analysis_scripts/7_cellphase_model.ipynb
   scripts/3_analysis_scripts/9A_cnv_score.ipynb
   scripts/3_analysis_scripts/9B_velocity_loom.ipynb
   scripts/3_analysis_scripts/11_spectra.ipynb
   scripts/3_analysis_scripts/10_survival_curve_isoform

Script Categories
-----------------

**Preprocessing & Guide Calling (1, 2A, 2B, 12)**
   * Cell Ranger processing and guide calling workflows
   * **Shell Pipeline:** :download:`1_cellranger.sh <scripts/3_analysis_scripts/1_cellranger.sh>`

**Transcript Quantification (3A-3C)**
   * Isoform-specific expression analysis using Whippet
   * **Shell Scripts:** :download:`3A_seperate_bam_file_ideal.sh <scripts/3_analysis_scripts/3A_seperate_bam_file_ideal.sh>` | :download:`3B_whippet_pseudobulk.sh <scripts/3_analysis_scripts/3B_whippet_pseudobulk.sh>`

**Transcriptional Divergence (4-6, 8)**
   * E-statistics, differential expression, and UMAP visualization

**Cell Phase & Pathway Modeling (7, 9A-9B, 11)**
   * Cell cycle analysis, CNV scoring, and pathway modeling

**Clinical Analysis (10)**
   * Survival curve and clinical correlation analysis