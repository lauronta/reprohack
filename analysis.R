library(DESeq2)

help(DESeq2)

setwd("C:/Users/rapha/OneDrive/Bureau/IODAA/Cours Université Paris-Saclay/Hackathon")
getwd()


## DESeqDataSet --> dds

# 4 different input data type for DESeqDataSet creation 
# here we focus on count matrix of read counts (with featureCounts) input

data = read.table("counts.txt", header = TRUE, skip =1)

# drop first columns in featureCounts file, keep only the counts_data_matrix and Geneid for orthology
data <- data[, -seq(1, 6)]


# definition of the experimental plan
sample <- c("SRR10379721.bam", "SRR10379722.bam", "SRR10379723.bam", "SRR10379724.bam", "SRR10379725.bam", "SRR10379726.bam")
coldata <- matrix(c("persister replicate 1", "persister replicate 2", "persister replicate 3", "control replicate 1", "control recplicate 2", "control replicate 3"), dimnames = list(sample, 'condition'))
print(coldata)


# first if it is a control or test group -- replace with colnames of data
coldata <- matrix( c("persister", "persister", "persister", "control", "control", "control"), 
                   dimnames = list(c("SRR10379721.bam", "SRR10379722.bam", "SRR10379723.bam", "SRR10379724.bam", "SRR10379725.bam", "SRR10379726.bam"), 'condition') )
print(coldata)

dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition)



# correcting variance
# vst (n<30 ) or rlog (n>30) approach
# use 


# performing differencial analysis 
dds <- DESeq(dds)

# head(data)
res <- results(dds)

plotMA(res, colSig = "red", ylim=c(-4.25, 4.25))

# Import de AureoWiki
aureowiki = read.table("OrthologueTable.tsv", header = TRUE, skip = 0, fill = TRUE)
aureowiki <- aureowiki[, -seq(3,4)] # suppression d'une colonne inutile
