#Tip: To run any line in RStudio you can use Ctrl+Enter or Command+Enter.
#To use any package in R, we need to load its library first.To run DESeq2 
#you just need to type the following.
library("DESeq2")
#We are going to use other libraries and lets load them too.
library("ggplot2")
library("RColorBrewer")
library("gplots")

#1. Reading the data.
# We already made gene quantifications and merged expected counts into 
# a single table. Here we prepared another table using whole reads 
# (not only reduced reads) to give you an idea about the complete 
# picture of DESeq analysis.

file <- "~/RNASeqWS/transcriptomics/data/deseq_dataset.tsv"
rsem<-read.table(file)
head(rsem)


#2. Read the data with header. We set the header=TRUE. In this way, 
# the first row will be in the header of the data structure.
rsem <- read.table(file,sep="\t", header=TRUE)


#3. Read the data with row.names. We have the gene names in the 
# first column in this table. When we set row.names=1 It will remove 
# the first column from the data and put them to the row name section.

rsem <- read.table(file,sep="\t", header=TRUE, row.names=1,
                  stringsAsFactors = TRUE)

#4. Creating the data structure for DESeq Analysis, we need to select 
# the columns that we are going to use in the analysis. Here we are 
# going to use experiment and control data that will be 6 coulmns. 
# Here in this step make sure to include the columns you want to use in 
# your analysis.

data <- data.frame(rsem[, c("exper_rep1","exper_rep2","control_rep3",
                            "control_rep1","control_rep2","exper_rep3")])

#5. DESeq analysis uses a Negative Binomial distribution to model 
# fragment COUNTS not TPMs. We used RSEM to estimate counts per transcript 
# or gene. RSEM uses a likelihood model, which results in non-integer 
# fragment coverage. We need to convert those values to integer. 

cols <- c(1:6);
data[,cols] <- apply(data[,cols], 2,
                    function(x) as.numeric(as.integer(x)))

#6. Just to see the distribution of the reads we can use hist function in 
# log10 scale.

hist(log10(data$exper_rep1), breaks=100)

#7. Scatter plots are a key approach to assess differences between conditions. 
# We plot the average expression value (counts or TPMs for example) of 
# two conditions.First calculate the average reads per gene. To find the 
# average we sum experiments per gene and divide it to the number of replicas

# In this case we will divide the sum to 3. We will calculate the average for 
# control libraries too. After that, we are ready to merge  the average values 
# from experiment and control to create a two column  data structure using the 
# cbind function that merges tables.

avgall<-cbind(rowSums(data[c("exper_rep1","exper_rep2","exper_rep3")])/3, 
              rowSums(data[c("control_rep1","control_rep2","control_rep3")])/3)

# We can also change the column names using colnames function.

colnames(avgall)<-c("Treat", "Control")

# Make a simple scatter plot.

plot(avgall)

#8. Hmm!!! The values are ranging from 0 to 800k. So, let's use log2.

plot(log2(avgall))

#9. Let's change the x and y labels and plot again.

log2vals <- log2(avgall)
colnames(log2vals)<-c("log2(Treat)", "log2(Control)")
plot(log2vals)

#10.There are several packages in R that can be used to 
# generate scatter plots. The same plot can be generated using 
# the ggplot package.

gdat<-data.frame(avgall)
ggplot() +
  geom_point(data=gdat, aes_string(x="Treat", y="Control"),
             colour="black", alpha=6/10, size=3) +
  scale_x_log10() +scale_y_log10()

######################################
## LETS START DESeq ANALYSIS
######################################
#
#11. The goal of Differential gene expression analysis is to find 
# genes or transcripts whose difference in expression, when accounting 
# for the variance within condition, is higher than expected by chance. 

# The first step is to indicate the condition that each column (experiment) 
# in the table represent. 

# Here we define the correspondence between columns and conditions. 
# Make sure the order of the columns matches to your table.


conds <- factor( c("Control","Control", "Control",
                   "Treat", "Treat","Treat") )

# DESeq function requires a special data sturcture called data frame.

colData <- as.data.frame((colnames(data)))
colData <- cbind(colData, conds)
colnames(colData) <- c("Libs","group")
groups <- factor(colData[,2])

#12. In Eukaryotes only a subset of all genes are expressed in 
# a given cell. Expression is therefore a bimodal distribution, 
# with non-expressed genes having counts that result from experimental 
# and biological noise. It is important to Filter out the genes 
# that are not expressed before doing differential gene expression. 

# You can decide which cutoff separates expressed vs non-expressed 
# genes by looking your histogram we created.

# In our case a total sum of 10 counts separates well expressed 
# from non-expressed genes


sumd <- rowSums(data)
hist(log10(sumd), breaks=100)
abline(v=1)

# To select the sum of the rows > 10, we can use subset function

filtd <- subset(data, sumd > 10)

#13. Create DESeq data set using prepared colData table for condition 
# definition and filtered data.
dds <- DESeqDataSetFromMatrix(countData=as.matrix(filtd),
                             colData=colData, design = ~ group)

#14. Run DESeq analysis
dr <- DESeq(dds);

#15. Put the results into variable.
res <- results(dr);

#16. Look for the column descriptions in the results
mcols(res, use.names=TRUE)


#17. DESeq will compute the probability that a gene is differentially 
# expressed (DE) for ALL genes in the table. It outputs both a 
# nominal and a multiple hypothesis corrected p-value (padj). 
# Because we are testing DE for over 20,000 genes, we do need to correct 
# for multiple hypothesis testing. To find genes that are significantly 
# DE select the ones has lower padj values. # higher fold changes and 
# visualize them on our scatter plot with different  # color. 
# padj values are corrected p-values which are multiplied by the number 
# of comparisons. Here we are going to use 0.01 for padj value 
# and > 1 log2foldchange.  (1/2 < foldChange < 2)

f1 <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]
res_selected <- f1[(f1$padj<0.01 & abs(f1$log2FoldChange)>1),]

#18. To Add a legend for all data we are going to add a text "Add" to 
# whole values.
Legend <- "All"
gdat1 <- cbind(gdat, Legend)
gdat_selected <- gdat[rownames(res_selected),]

#19. Add a legend for only significant ones that we selected.
Legend <- "Significant"
gdat_selected1 <- cbind(gdat_selected, Legend)

#20. We now need to merge selected and all data to draw them 
# together on a plot with different colors.
gdat2 <- rbind(gdat1, gdat_selected1)

#21. Make ggplot

ggplot() +
  geom_point(data=gdat2, aes_string(x="Treat", y="Control",
             colour="Legend"), alpha=6/10, size=3) +
  scale_colour_manual(values=c("All"="darkgrey","Significant"="red"))+
  scale_x_log10() +scale_y_log10()

#22. There is another way to visualize the same data using MAPlots.
# MA plots are basically a plot with log ratios(M) in the y axis vs.
# mean(A) scale in x axis.

gdatMA <- gdat2

#We can just add pseudocount value to prevent calculating infinite values 
# if the a gene quantification is 0. This transformation will not change 
# anything in data distribution.

gdatMA$Treat <- gdatMA$Treat+1
gdatMA$Control <- gdatMA$Control+1
g <- gdatMA
colnames(g) <- c("A", "M", "Legend")

# Here A <- ( log2(x) + log2(y) )/2  (Average in log scale)
g$A <- log2(gdatMA$Treat*gdatMA$Control/2)

# Here M <- log2(x) - log2(y) (Minus in log scale)
g$M <- log2(gdatMA$Treat/gdatMA$Control)

ggplot() +
  geom_point(data=g, aes_string(x="A", y="M", colour="Legend"),
             alpha=6/10, size=3) + ylab("Log2 Fold Change (M)") +
  xlab("Log2 mean normalized counts (A)") +
  scale_colour_manual(values=c("All"="black","Significant"="red"))+
  geom_abline(slope=0, linetype=2)

#23. If you want to save any ggplot as pdf use the command below
ggsave("~/workshop_data/transcriptomics/MAplot.pdf")

#24. For MA Plot there is another builtin function that you can use.
plotMA(dr,ylim=c(-2,2),main="DESeq2");


#25. The third way of visualizing the data is making a Volcano Plot.
# Here on the x axis you have log2foldChange values and y axis you 
# have your -log10 padj values. To see how significant genes are 
# distributed. Highlight genes that have an absolute fold change > 2 
# and a padj < 0.01

res$threshold = as.factor(abs(res$log2FoldChange) > 1 & res$padj < 0.01)
res$log10padj = -log10(res$padj)
dat<-data.frame(cbind(res$log2FoldChange, res$log10padj, res$threshold))

#Define your column names
colnames(dat)<-c("log2FoldChange", "log10padj", "threshold")

##Construct the plot object
ggplot(data=dat, aes_string(x="log2FoldChange", y="log10padj",
              colour="threshold")) + geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-2.5, 2.5)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")

#26. Adding more information to your dataset and save as an excel file
# You can select the significant genes from the dataset and add their 
# log2foldChange and padj values to the dataset

calc_cols<-res[, c("log2FoldChange", "padj")]
f1<-cbind(data[rownames(calc_cols), ], calc_cols)
f1<-f1[!is.na(res$padj), ]
res_selected<-f1[f1$padj<0.01 & abs(f1$log2FoldChange)>1, ]
write.csv(res_selected, "~/workshop_data/transcriptomics/selected_genes.csv")

# You can write all of your data to a csv file.
write.csv(data, "~/workshop_data/transcriptomics/all_data.csv")

#27. The forth way of visualizing the data that is widely used in this 
# type of analysis is clustering and Heatmaps.
# Here we usually add a pseudocount value 0.1 in this case.

ld <- log2(filtd[rownames(res_selected),]+0.1)
#Scaling the value using their mean centers can be good it the data is 
# uniformely distributed in all the samples.

cldt <- scale(t(ld), center=TRUE, scale=TRUE);
cld <- t(cldt)

# We can define different distance methods to calculate the distance
# between samples. Here we focus on euclidean distance and correlation
#a. Euclidean distance

distance<-dist(cldt, method = "euclidean")

# To plot only the cluster you can use the command below
plot(hclust(distance, method = "complete"),
     main="Euclidean", xlab="")
#The heatmap

heatmap.2(cld, Rowv=TRUE,dendrogram="column",
          Colv=TRUE, col=redblue(256),labRow=NA,
          density.info="none",trace="none", cexCol=0.8);

#b. Correlation between libraries
#Here we calculate the correlation between samples.
dissimilarity <- 1 - cor(cld)
# We define it as distance. This will create a square matrix that
# will include the dissimilarites between samples.
distance <- as.dist(dissimilarity)

# To plot only the cluster you can use the command below

plot(hclust(distance, method = "complete"),
     main="1-cor", xlab="")

#The heatmap
heatmap.2(cld, Rowv=TRUE,dendrogram="column",
          Colv=TRUE, col=redblue(256),labRow=NA,
          density.info="none",trace="none", cexCol=0.8)
