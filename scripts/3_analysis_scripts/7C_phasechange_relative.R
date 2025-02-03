

phase_assignment_gene <- read.csv("alt-prom-crispr-fiveprime/files/phase_assignment_gene.csv", header = TRUE, sep = ",", row.names = 1)
phase_assignment_gene_S <- phase_assignment_gene[(phase_assignment_gene$phase=="S")&(phase_assignment_gene$target_gene!="non-targeting"),]
phase_assignment_gene_S$relative_percentage <- phase_assignment_gene_S$percentage - phase_assignment_gene$percentage[(phase_assignment_gene$phase=="S")&(phase_assignment_gene$target_gene=="non-targeting")]
phase_assignment_gene_S <- phase_assignment_gene_S[order(phase_assignment_gene_S$relative_percentage),]
phase_assignment_gene_S$relative_percentage  <- as.numeric(as.character(phase_assignment_gene_S$relative_percentage))
phase_assignment_gene_S$target_gene_promoter <- paste0(phase_assignment_gene_S$target_gene, "_", phase_assignment_gene_S$promoter_type)
phase_assignment_gene_S

ggbarplot(phase_assignment_gene_S, x = "target_gene_promoter", y = "relative_percentage",
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts

)

#repeat G1
phase_assignment_gene_G1 <- phase_assignment_gene[(phase_assignment_gene$phase=="G1")&(phase_assignment_gene$target_gene!="non-targeting"),]
phase_assignment_gene_G1$relative_percentage <- phase_assignment_gene_G1$percentage - phase_assignment_gene$percentage[(phase_assignment_gene$phase=="G1")&(phase_assignment_gene$target_gene=="non-targeting")]
phase_assignment_gene_G1 <- phase_assignment_gene_G1[order(phase_assignment_gene_G1$relative_percentage),]
phase_assignment_gene_G1$relative_percentage  <- as.numeric(as.character(phase_assignment_gene_G1$relative_percentage))
phase_assignment_gene_G1$target_gene_promoter <- paste0(phase_assignment_gene_G1$target_gene, "_", phase_assignment_gene_G1$promoter_type)
phase_assignment_gene_G1

ggbarplot(phase_assignment_gene_G1, x = "target_gene_promoter", y = "relative_percentage",
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts

)

#repeat G2
phase_assignment_gene_G2 <- phase_assignment_gene[(phase_assignment_gene$phase=="G2")&(phase_assignment_gene$target_gene!="non-targeting"),]
phase_assignment_gene_G2$relative_percentage <- phase_assignment_gene_G2$percentage - phase_assignment_gene$percentage[(phase_assignment_gene$phase=="G2")&(phase_assignment_gene$target_gene=="non-targeting")]
phase_assignment_gene_G2 <- phase_assignment_gene_G2[order(phase_assignment_gene_G2$relative_percentage),]
phase_assignment_gene_G2$relative_percentage  <- as.numeric(as.character(phase_assignment_gene_G2$relative_percentage))
phase_assignment_gene_G2$target_gene_promoter <- paste0(phase_assignment_gene_G2$target_gene, "_", phase_assignment_gene_G2$promoter_type)
phase_assignment_gene_G2

ggbarplot(phase_assignment_gene_G2, x = "target_gene_promoter", y = "relative_percentage",
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts

)
