# In R on Biomix
library(GenomicFeatures)
gtf <- "data/reference/mmus.gtf.gz"
txdb <- makeTxDbFromGFF(gtf, format="gtf")
k <- keys(txdb, keytype="TXNAME")
map <- select(txdb, keys=k, keytype="TXNAME", columns="GENEID")
tx2gene <- unique(data.frame(transcript=map$TXNAME, gene=map$GENEID))
write.csv(tx2gene, "data/reference/tx2gene.csv", row.names=FALSE)
