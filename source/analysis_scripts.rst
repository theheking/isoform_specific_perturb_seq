Analysis Scripts & Pipelines
============================

This section contains the core computational notebooks and scripts.



Preprocessing & Guide Calling
-----------------------------
* **Cellranger Pipeline:** :download:`1_cellranger.sh <./scripts/3_analysis_scripts/1_cellranger.sh>`
* **Guide Calling (Primary):** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/2A_cellranger_guidecalling>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb>`
* **Guide Calling (Ideal Guides):** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides.ipynb>`

Transcriptional Divergence (E-Stats)
------------------------------------
* **E-Statistic Calculation:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/4_estatistic_gene_guide>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/4_estatistic_gene_guide.ipynb>`
* **Neighbouring Gene Expression:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression.ipynb>`
* **UMAP & KL-Divergence:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap.ipynb>`

Cell Phase & Pathway Modeling
-----------------------------
* **Cell Phase Model:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/7A_cellphase_model>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/7A_cellphase_model.ipynb>`
* **Spectra Pathway Analysis:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/11A_spectra>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/11A_spectra.ipynb>`
* **CNV Score Analysis:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/9A_cnv_score>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/9A_cnv_score.ipynb>`

Clinical & Differential Expression
----------------------------------
* **Differential Expression:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/8_perturbseq_differentialexpression>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/8_perturbseq_differentialexpression.ipynb>`
* **ESR1 Survival Curves:** :doc:`Visualize Notebook <./scripts/3_analysis_scripts/10D_esr1_survival_curve>` | :download:`Download .ipynb <./scripts/3_analysis_scripts/10D_esr1_survival_curve.ipynb>`

.. toctree::
   :hidden:
   :maxdepth: 1

   ../scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb
   ../scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides.ipynb
   ../scripts/3_analysis_scripts/4_estatistic_gene_guide.ipynb
   ../scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression.ipynb
   ../scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap.ipynb
   ../scripts/3_analysis_scripts/7A_cellphase_model.ipynb
   ../scripts/3_analysis_scripts/8_perturbseq_differentialexpression.ipynb
   ../scripts/3_analysis_scripts/9A_cnv_score.ipynb
   ../scripts/3_analysis_scripts/10D_esr1_survival_curve.ipynb
   ../scripts/3_analysis_scripts/11A_spectra.ipynb