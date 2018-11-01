library("DESeq2");
setwd("~/RNASeqWS/transcriptomics");
d = read.table("rsem.gene.summary.count.txt", header=TRUE, row.names=1);

colData = as.data.frame((colnames(d)));
colData[1,2] = "control";
colData[2,2] = "expr";
colData[3,2] = "expr";
colData[4,2] = "control";
colData[5,2] = "expr";
colData[6,2] = "control";
colnames(colData) = c("Sample","group");
groups = factor(colData[,2]);

sumd = apply(X=d,MARGIN=1,FUN=sum);

filtd = subset(d, sumd > 6);
dds = DESeqDataSetFromMatrix(countData= as.matrix(filtd), colData=colData, design = ~ group);
dds <- DESeq(dds);
res <- results(dds);
mcols(res, use.names=TRUE)


######################################
## Visualizing a little
######################################
library("RColorBrewer");
library("gplots");
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100);

plotMA(dds,ylim=c(-2,2),main="DESeq2");
plotDispEsts(dds);

heatmap.2(counts(dds,normalized=TRUE), col = hmcol,
          Rowv = TRUE, Colv = TRUE, scale="row",
          dendrogram="both", trace="none", margin=c(10,6))

write.table(results(dds), file="DESeq.results.txt", quote=FALSE, sep="\t");      

# Lets now compare the normalized vs the raw counts
m =as.matrix(d);
nm = counts( dds, normalized=TRUE );
colnames(m) = colnames(d);
colnames(nm) =  colnames(d);

#Lets create a small simple routine to plot our data:
plot.counts <- function(um, nm, col) {
  min = min (c(um[col, ],nm[col, ]));
  max = max (c(um[col, ],nm[col, ]));
  
  plot(um[col, ], xaxt = "n", ylab="Counts", main = col, xlab="",pch=18, ylim=c(min, max), col="red");
  points(nm[col, ], ,pch=18,  col="black");
  axis(labels=colnames(m), side=1, at=c(1:6),cex.axis=0.75, las = 2);
}



plot.counts(m, nm, "Bcat2")
plot.counts(m, nm, "AK208554")
plot.counts(m, nm,"0610005C13Rik");
plot.counts(m, nm, "Nat15");
plot.counts(m,nm,"Crebbp");
plot.counts(m,nm,"Fgf21")
plot.counts(m,nm,"Coro7")

