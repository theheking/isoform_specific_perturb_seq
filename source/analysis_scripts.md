Analysis Scripts
================

This section contains the core processing pipeline for APU Perturb-seq.

Pre-processing & Quantification
-------------------------------
* **Cellranger:** :download:`1_cellranger.sh <../scripts/3_analysis_scripts/1_cellranger.sh>`
* **Guide Calling:** :download:`2A_cellranger_guidecalling.ipynb <../scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb>`
* **Guide Assignment Visuals:** :download:`2B_upset_plot_guide_assignment.r <../scripts/3_analysis_scripts/2B_upset_plot_guide_assignment.r>`



Alternative Splicing & Pseudobulk
---------------------------------
* **BAM Separation:** :download:`3A_seperate_bam_file_ideal.sh <../scripts/3_analysis_scripts/3A_seperate_bam_file_ideal.sh>`
* **Whippet Processing:** :download:`3B_whippet_pseudobulk.sh <../scripts/3_analysis_scripts/3B_whippet_pseudobulk.sh>`
* **Quality Control:** :download:`3C_whippet_qc.R <../scripts/3_analysis_scripts/3C_whippet_qc.R>`

Downstream Analysis
-------------------
* **Differential Expression:** :download:`8_perturbseq_differentialexpression.R <../scripts/3_analysis_scripts/8_perturbseq_differentialexpression.R>`
* **Cell Cycle Modeling:** :download:`7A_cellphase_model.ipynb <../scripts/3_analysis_scripts/7A_cellphase_model.ipynb>`
* **Spectra Visualisation:** :download:`11B_spectra_visualisation.R <../scripts/3_analysis_scripts/11B_spectra_visualisation.R>`