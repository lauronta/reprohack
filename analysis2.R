library(DESeq2)

help(DESeq2)

library("dplyr")
library(jsonlite)
library(ggplot2)

setwd("C:/Users/rapha/OneDrive/Bureau/IODAA/Cours Université Paris-Saclay/Hackathon")
getwd()


###############################################################################################"
# téléchargement gène impliqué dans les différentes fonction depuis KEGG-rest API
data <- fromJSON("https://rest.kegg.jp/get/br:sao00001/json")

# retourne dernier sous-élément d'une liste, ici el
f <- function(el) el[length(el)]

get_gene_name <- function(data,clef1,clef2,clef3){
  # on descend dans la hiérarchie
  col1 <- data$children$name
  print(col1)
  
  # on sélectionne selon le mask
  mask <- col1 == clef1
  data <- data$children$children[mask][[1]]
  
  #on descend encore dans la hiérachie 
  col2 <- data$name
  print(col2)
  
  # on sélectionne le masque deux
  mask <- col2 == clef2
  data <- data$children[mask][[1]]
  
  col3 <- data$name
  print(col3)
  
  
  mask <- col3 == clef3
  data <- data$children[mask][[1]]
  
  # split selon ;
  data <- data %>% mutate(name = sapply(strsplit(as.character(data$name), split=";"), "[",1))
  
  
  # split pour prendre la première colonne le code du gène et le dernier élément le nom du gène
  
  nom <- sapply(strsplit(as.character(data$name), split=" "), f)
  code_gene <- sapply(strsplit(as.character(data$name), split=" "), "[",1)
  
  data_kegg <- data.frame(code_gene, nom)
  
  return(data_kegg)
}


clef1 <- "09180 Brite Hierarchies"
clef2 <- "09182 Protein families: genetic information processing"
clef3 <- "03012 Translation factors [BR:sao03012]"
translation_factor_kegg <- get_gene_name(data,clef1,clef2,clef3)

clef1 <- "09120 Genetic Information Processing"
clef2 <- "09122 Translation"
clef3 <- "03010 Ribosome [PATH:sao03010]"
ribosome_kegg <- get_gene_name(data,clef1,clef2,clef3)

clef3 <- "00970 Aminoacyl-tRNA biosynthesis [PATH:sao00970]"
aa_synthetase_kegg <- get_gene_name(data,clef1,clef2,clef3)

#clef3 <- "03013 Nucleocytoplasmic transport"
#nucleo_transport_kegg <- get_gene_name(data,clef1,clef2,clef3)

#clef3 <- "03015 mRNA surveillance pathway"
#mrna_surveillance_kegg <- get_gene_name(data,clef1,clef2,clef3)

clef3 <- "03008 Ribosome biogenesis in eukaryotes"
ribosome_synthese_kegg <- get_gene_name(data,clef1,clef2,clef3)

print(translation_factor_kegg)             

# Construction de la liste des gènes impliquées dans la translation : ATTENTION UTILISATION DPLYR
trans_genes = rbind(ribosome_kegg, translation_factor_kegg,
                    aa_synthetase_kegg, ribosome_synthese_kegg)

print(trans_genes)
######################################################################################################

## DESeqDataSet --> dds

# 4 different input data type for DESeqDataSet creation 
# here we focus on count matrix of read counts (with featureCounts) input

data = read.table("counts.txt", header = TRUE, skip =1)

# Arrangement de la colonne Geneid pour supprimer le terme gene pour correspondance avec transgene
data$Geneid <- sapply(strsplit(as.character(data$Geneid), split="-"), "[", 2)

# Supprimer les lignes où 'Geneid' ne correspond pas à 'trans_genes'
data_filtrée <- data[data$Geneid %in% trans_genes$code_gene, ]

# Remplacement de la colonne "Geneid" de data par nom de "trans_genes", index respecté
data_filtrée$Geneid <- sapply(data_filtrée$Geneid, function(gene) {
  # Trouver le nom du gène correspondant à chaque code_gene
  matched_gene <- trans_genes$nom[trans_genes$code_gene == gene]
  if(length(matched_gene) > 0) {
    return(matched_gene)  # Remplacer par le nom du gène
  } else {
    return(NA)  # Si pas de correspondance, on met NA
  }
})

# drop first columns in featureCounts file, keep only the counts_data_matrix and Geneid for orthology
data = data_filtrée
data <- data[, -seq(2, 6)]


# definition trans_genes# definition of the experimental plan
sample <- c("SRR10379721.bam", "SRR10379722.bam", "SRR10379723.bam", "SRR10379724.bam", "SRR10379725.bam", "SRR10379726.bam")
coldata <- matrix(c("persister replicate 1", "persister replicate 2", "persister replicate 3", "control replicate 1", "control recplicate 2", "control replicate 3"), dimnames = list(sample, 'condition'))
print(coldata)


# first if it is a control or test group -- replace with colnames of data
coldata <- matrix( c("persister", "persister", "persister", "control", "control", "control"), 
                   dimnames = list(c("SRR10379721.bam", "SRR10379722.bam", "SRR10379723.bam", "SRR10379724.bam", "SRR10379725.bam", "SRR10379726.bam"), 'condition') )
print(coldata)

dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition, tidy = TRUE)
print(dds)

# correcting variance
# vst (n<30 ) or rlog (n>30) approach
# use 


# performing differencial analysis 
dds <- DESeq(dds)

# head(data)
res <- results(dds)
head(results(dds, tidy=TRUE))

plotMA(res, colSig = "red", ylim=c(-4.25, 4.25))

# Convertion les résultats en data.frame pour ggplot
res_df <- as.data.frame(res)
head(res_df)

# :og2 à la colonne baseMean directement dans le data.frame
res_df$log2_baseMean <- log2(res_df$baseMean)

# Remplacement de la colonne "Geneid" de data par nom de "trans_genes", index respecté
rownames(res_df) <- sapply(rownames(res_df), function(gene) {
  # Trouver le nom du gène correspondant à chaque code_gene
  matched_gene <- trans_genes$nom[trans_genes$code_gene == gene]
  if(length(matched_gene) > 0) {
    return(matched_gene)  # Remplacer par le nom du gène
  } else {
    return(NA)  # Si pas de correspondance, on met NA
  }
})

# Définition d'une liste de gènes à afficher
genes_to_annotate <- c("frr", "inf A", "inf B", "tsf", "inf C", "pth")

# Définition d'une liste d'annotation pour récupérer les noms de gènes à placer dans le graph
#rownames(res_df) %in% genes_to_annotate
#annotation_df <- res_df[rownames(res_df) %in% genes_to_annotate, ]
#print(annotation_df)

# Générer le MA plot
ggplot(res_df, aes(x = log2_baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(is.na(padj), "NA", 
             ifelse(padj < 0.05, "Significant", "Not Significant"))), 
             alpha = 0.9) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  theme_classic() +
  labs(title = "MA Plot", x = "log2(Base Mean)", y = "log2(Fold Change)") +
  theme(legend.title = element_blank()) +
  # xlim(c(0, max(res_df$log2_baseMean))) +  # Limites de l'axe X basées sur log2(baseMean)
  xlim(c(0, 20)) +
  ylim(c(-6, 5))  # Limites de l'axe Y

