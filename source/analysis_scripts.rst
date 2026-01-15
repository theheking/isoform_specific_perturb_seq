Analysis Scripts & Pipelines
============================

This section contains the core computational notebooks and scripts used for data processing, statistical analysis, and visualization.

Preprocessing & Guide Calling
-----------------------------

* **Cellranger Pipeline:** :download:`1_cellranger.sh <./scripts/3_analysis_scripts/1_cellranger.sh>`
* **Guide Calling (Primary):** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/2A_cellranger_guidecalling>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb>`
* **Guide Calling (Ideal Guides):** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides.ipynb>`
* **Negative Control Analysis:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/12_negativecontrol>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/12_negativecontrol.ipynb>`

Splicing & RNA Velocity
-----------------------

* **BAM File Separation:** :download:`3A_seperate_bam_file_ideal.sh <./scripts/3_analysis_scripts/3A_seperate_bam_file_ideal.sh>`
* **Whippet Pseudobulk Pipeline:** :download:`3B_whippet_pseudobulk.sh <./scripts/3_analysis_scripts/3B_whippet_pseudobulk.sh>` | :doc:`View R Analysis <./scripts/3_analysis_scripts/3C_whippet_pseudobulk>`
* **RNA Velocity (Loom Generation):** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/9B_velocity_loom>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/9B_velocity_loom.ipynb>`

Transcriptional Divergence (E-Stats)
------------------------------------

* **E-Statistic Calculation:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/4_estatistic_gene_guide>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/4_estatistic_gene_guide.ipynb>`
* **Neighbouring Gene Expression:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression.ipynb>`
* **UMAP & KL-Divergence:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap.ipynb>`

Cell Phase & Pathway Modeling
-----------------------------

* **Cell Phase Model:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/7_cellphase_model>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/7_cellphase_model.ipynb>`
* **Spectra Pathway Analysis:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/11_spectra>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/11_spectra.ipynb>`
* **CNV Score Analysis:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/9A_cnv_score>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/9A_cnv_score.ipynb>`

Clinical & Differential Expression
----------------------------------

* **Differential Expression:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/8_perturbseq_differentialexpression>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/8_perturbseq_differentialexpression.ipynb>`
* **Survival Curve Analysis:** :doc:`View Rmd <./scripts/3_analysis_scripts/10_survival_curve>` | :download:`Download .Rmd <./scripts/3_analysis_scripts/10_survival_curve.Rmd>`

.. toctree::
   :hidden:
   :maxdepth: 1

   ../scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb
   ../scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides.ipynb
   ../scripts/3_analysis_scripts/3C_whippet_pseudobulk.Rmd
   ../scripts/3_analysis_scripts/4_estatistic_gene_guide.ipynb
   ../scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression.ipynb
   ../scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap.ipynb
   ../scripts/3_analysis_scripts/7_cellphase_model.ipynb
   ../scripts/3_analysis_scripts/8_perturbseq_differentialexpression.ipynb
   ../scripts/3_analysis_scripts/9A_cnv_score.ipynb
   ../scripts/3_analysis_scripts/9B_velocity_loom.ipynb
   ../scripts/3_analysis_scripts/10_survival_curve.Rmd
   ../scripts/3_analysis_scripts/11_spectra.ipynb
   ../scripts/3_analysis_scripts/12_negativecontrol.ipynb