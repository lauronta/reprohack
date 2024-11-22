library(DESeq2)

help(DESeq2)


## DESeqDataSet --> dds

# 4 different input data type for DESeqDataSet creation 
# here we focus on count matrix of read counts (with featureCounts) input

data = read.table("test_counts (2).txt", header = TRUE, skip =1)

# drop first columns in featureCounts file, keep only the counts_data_matrix
data <- data[, -seq(1, 6)]


# definition of the experimental plan
sample <- c("SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726")
coldata <- matrix(c('persister replicate 1', "persister replicate 2", "persister replicate 3", "control replicate 1", "control recplicate 2", "control replicate 3"), dimnames = list(sample, 'condition'))
print(coldata)


# first if it is a control or test group -- replace with colnames of data
coldata <- matrix(c("persister", "persister", "persister", "control", "control", "control"), dimnames = list(c("test_21.bam", "test_22.bam", "test_23.bam", "test_24.bam", "test_25.bam", "test_26.bam"), "condition"))
print(coldata)

dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition)



# correcting variance
# vst (n<30 ) or rlog (n>30) approach
# use 


# performing differencial analysis 
dds <- DESeq(dds)

# head(data)
res <- results(dds)

plotMA(res, colSig = "red", ylim=c(-6,5))
# 
