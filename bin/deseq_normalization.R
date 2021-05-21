#! /usr/bin/env Rscript

library("DESeq2")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")

# Parse CLI args
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Please input the location of the sample information file, and the gene_counts table", call.=FALSE)
}

# Parse sample informations
condition = read.csv(file = args[1], sep=",", header=T)$condition

# Parse gene counts
cts = read.csv(args[2], sep="\t")
coldata = data.frame(condition, row.names = colnames(cts))
coldata

# Create DESeq object and normalize
if(length(unique(condition)) == 1 ) {
	dds <- DESeqDataSetFromMatrix(countData = cts,
				      colData = coldata,
				      design = ~ 1)
	vsd <- vst(dds, nsub=100, blind=FALSE)
} else {
	dds <- DESeqDataSetFromMatrix(countData = cts,
				      colData = coldata,
				      design = ~ condition)
	vsd <- vst(dds, nsub=100, blind=FALSE)


	# Heatmap 1
	dds <- estimateSizeFactors(dds)
	select <- order(rowMeans(counts(dds,normalized=TRUE)),
			decreasing=TRUE)[1:40]
	df <- as.data.frame(colData(dds)[,c("condition")])
	rownames(df) = colnames(cts)

	# Heatmap 2
	sampleDists <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

	# Create pdf with plots
	pdf("normalization.pdf")
	plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name), vjust=2, check_overlap = TRUE, size = 3)
	pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
		 cluster_cols=FALSE, annotation_col=df)
	pheatmap(sampleDistMatrix,
		 clustering_distance_rows=sampleDists,
		 clustering_distance_cols=sampleDists,
		 col=colors)
	dev.off()
}

# Write normalized table
write.table(as.data.frame(assay(vsd)), file="counts_gene.normalized.vsd.csv", quote=FALSE, sep="\t")
