library(DESeq2)
library("dplyr")
library(jsonlite)
library(ggplot2)
library(ggrepel)

setwd("~/Bureau/IODAA/iodaa/HACKATHON/reprohack/")


################################################################################
# download translation-related genes from KEGG-rest API
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
  
  name <- sapply(strsplit(as.character(data$name), split=" "), f)
  geneID <- sapply(strsplit(as.character(data$name), split=" "), "[",1)
  
  data_kegg <- data.frame(geneID, name,label=clef3)
  
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
# empty so commented 

#clef3 <- "03015 mRNA surveillance pathway"
# empty so commented 

clef3 <- "03008 Ribosome biogenesis in eukaryotes"
ribosome_synthese_kegg <- get_gene_name(data,clef1,clef2,clef3)

# Construction de la liste des gènes impliquées dans la translation : ATTENTION UTILISATION DPLYR
gene_names <- rbind(ribosome_kegg, translation_factor_kegg,
                    aa_synthetase_kegg, ribosome_synthese_kegg)

head(gene_names)

################################################################################

## DESeqDataSet --> dds

# 4 different input data type for DESeqDataSet creation 
# here we focus on count matrix of read counts (with featureCounts) input

data = read.table("counts.txt", header = TRUE, skip =1)

# Arrangement de la colonne Geneid  :  gene-SAOUHSC_00001 --> SAOUHSC_00001 
data$Geneid <- sapply(strsplit(as.character(data$Geneid), split="-"), "[", 2)

# Supprimer les lignes où 'Geneid' ne correspond pas à 'trans_genes'
data_filtrée <- data[data$Geneid %in% gene_names$geneID, ]

# drop comuns in counts file, keep only Geneid and the counts_data_matrix
data = data_filtrée
data <- data[, -seq(2, 6)]

#definition of the experimental plan
sample <- c("SRR10379721.bam", "SRR10379722.bam", "SRR10379723.bam", 
            "SRR10379724.bam", "SRR10379725.bam", "SRR10379726.bam")
conditions <- c("persister", "persister", "persister", "control", "control",
                "control")
coldata <- matrix(conditions, dimnames = list(sample, 'condition') )
print(coldata)

# create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, 
                              design = ~ condition, tidy=T)
head(dds)

#performing differencial analysis 
dds <- DESeq(dds)
res <- results(dds)
head(results(dds))

# Convert DDS-->data.frame for plotting
res_df <- as.data.frame(res)
head(res_df)

# compute log2BaseMean
res_df$log2_baseMean <- log2(res_df$baseMean)

# add information from KEGG orthology previously downloaded
res_df$geneID <- rownames(res_df)
res_df <- left_join(res_df,gene_names, by = c("geneID"="geneID"))


# Définition d'une liste de gènes à afficher
genes_to_annotate <- c("frr", "infA", "infB", "tsf", "infC", "PTH1")

# Définition d'une liste d'annotation pour récupérer les noms de gènes à placer dans le graph
#rownames(res_df) %in% genes_to_annotate
#annotation_df <- res_df[rownames(res_df) %in% genes_to_annotate, ]
#print(annotation_df)

################################################################################
# Plotting
ggplot(res_df) +
  theme_bw()+
  geom_point(aes(x = log2_baseMean, y = log2FoldChange, color = ifelse(is.na(padj), "NA", 
              ifelse(padj < 0.05, "Significant", "Not Significant"))), 
              size=2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey50")) +
  geom_point(data= filter(res_df, label=='00970 Aminoacyl-tRNA biosynthesis [PATH:sao00970]' ), 
             aes(x = log2_baseMean, y = log2FoldChange,shape = 'AA-tRNA synthetases'), stroke=1, colour = 'black')+
  scale_shape_manual(values = c("AA-tRNA synthetases" = 1))+
  geom_text_repel(data=filter(res_df, name %in% genes_to_annotate),aes(x=log2_baseMean,y=log2FoldChange,label=name),size=8,box.padding = 5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(x = expression(Log[2] *" Base Mean"), y = expression(Log[2] *" Fold Change")) +
  theme(legend.title = element_blank(), 
        panel.border = element_rect(linewidth = 1.5, fill = NA,colour = 'black'),
        legend.position = c(0.01, .01),
        legend.justification = c("left", "bottom"),
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = 'horizontal',
        legend.margin = margin(6, 6, 0, 6),
        legend.text = element_text(colour = 'black', size = 12),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(color='black',size = 12),
        axis.ticks = element_line(color='black', linewidth = 1.5),
        axis.ticks.length = unit(2, "mm")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  coord_cartesian(
    xlim = c(0, 20),
    ylim = c(-6, 5),
    expand = FALSE,
    clip = "off") +
  scale_x_continuous(breaks=seq(0, 20, 2))+
  scale_y_continuous(breaks=seq(-6, 5, 1))

