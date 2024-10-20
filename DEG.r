library(DESeq2)
library("rjson")

exptable <- read.table("E:/rna-seq/rawcounts.txt", header = TRUE, quote = '\t')
rownames(exptable) <-exptable$gene
exptable <- exptable[,-1]


json.data <- fromJSON(file="E:/rna-seq/sample.json")


# convert to DESeq DataSet

table.condition <- data.frame(name = colnames(exptable), condition=config$conditions)
rownames(table.condition) <- table.condition$name
dds <- DESeqDataSetFromMatrix(exptable, colData=table.condition, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ] # filter out the zero count genes

# output normalized counts
dds <- estimateSizeFactors(dds)
counts.normalized <- counts(dds, normalized=TRUE)
write.table(counts.normalized, "./normalized_counts.tsv", col.names=NA, row.names=TRUE, sep="\t", quote=FALSE)

# DESeq analyze
dds <- DESeq(dds)
deseq.analyze <- function(compare) {
    results(dds, contrast=c("condition", compare[2], compare[1]))
}
deseq.results <- lapply(config$comparisons, deseq.analyze)
deseq.save <- function(idx) {
    compare.name <- paste(config$comparisons[[idx]], collapse="-")
    out.filename <- paste(compare.name, ".tsv", sep="")
    write.table(deseq.results[idx], out.filename,
                                    col.names=NA, row.names=TRUE, sep="\t",
                                    quote=FALSE)
}

