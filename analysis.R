library(DESeq2)
library(ggplot2)

file_path <- "counts.txt" #A modifier selon où est stocké counts.txt à la sortie
#du workflow

# Vérifie si le fichier existe dans le dossier spécifié
if (!file.exists(file_path)) {
  stop("Le fichier counts.txt est introuvable dans le dossier : ",
       print(file_path)
  )
}

# Chargement des données
data <- read.table(file_path, header = TRUE, skip = 1)

## DESeqDataSet --> dds

# 4 different input data type for DESeqDataSet creation 
# here we focus on count matrix of read counts (with featureCounts) input

# drop first columns in featureCounts file, keep only the counts_data_matrix and
#Geneid for orthology
data <- data[, -seq(1, 6)]


# definition of the experimental plan
sample <- c("SRR10379721.bam", "SRR10379722.bam", "SRR10379723.bam",
            "SRR10379724.bam", "SRR10379725.bam", "SRR10379726.bam")
conditions <- c("persister", "persister", "persister",
                "control", "control", "control")

# first if it is a control or test group -- replace with colnames of data
coldata <- matrix(conditions, dimnames = list(sample, 'condition'))
print(coldata)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ condition)


# performing differencial analysis 
dds <- DESeq(dds)

res <- results(dds)

# Affichage des résultats

### Solution 1, utiliser directement la fonction plotMA peu documentée et avec
#peu de paramètres

# plotMA(res,
#        alpha=0.05,
#        ylab=expression(log[2]*' fold change'),
#        colSig = "red", ylim=c(-4.25, 4.25)
#        )

### Solution 2, convertir le résultat de DESeq en dataframe,
# calculer la significativité, attribuer le label correspondant
# et afficher avec ggplot qui est bien plus documenté et maléable.

res_df <- as.data.frame(res)
alpha <- 0.05
res_df$significance <- ifelse(
  is.na(res_df$padj), "NA",
  ifelse(res_df$padj < alpha, "Significant", "Not Significant")
)

plt <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange,
                          color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(
    "Not Significant" = "black", 
    "Significant" = "red", 
    "NA" = "grey"
  )) +  
  labs(y = expression(log[2]*' fold change'),
       x = 'Mean of normalized counts') +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-4, 4, 2)) +
  theme_gray()

ggsave('Supplementary Figure 3.png',plt, scale=1)
