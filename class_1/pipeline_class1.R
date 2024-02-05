####
#loading the data and libraries
####
setwd("C:/Users/danie/OneDrive/Documents/bioinformatic_pipelines/class_1")
getwd()
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)


#install.packages("R.utils")
#install.packages("DESeq2")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")
#BiocManager::install("NOISeq")

library(R.utils)
library(BiocManager)
library(DESeq2)
library(GEOquery)
library(NOISeq)
BiocManager::install("apeglm")
library(apeglm)

######
## NOISEQ
######
###### https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf


## WE NEED THE METADATA (with GEOquery library)
#Se guardan solo las columnas de interes y los factores en el cual solo importa el estado de la enfermedad



## https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html

gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]

Factors_GSE198256 <- Meta_GSE198256[,c("disease state:ch1")]


## WE NEED BIOLOGICAL INFORMATION: GC, GENE LENGTH, CHROMOSOME,...

# Write the names
write.table(rownames(GSE198256_count),"gene_names.entrez.txt",
            col.names = FALSE,row.names = FALSE,quote=F)

# Additional Biological information.
# https://www.ensembl.org/biomart/martview/7f2a95d66853c3b8aea7639401e47aba

# Import the information
annotgene <- read.csv("mart_export.txt",sep="\t",header = T)
# How many genes do I get annotated?
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)

# Filter the information so only contains annotated gene in ch 1-22 and x, y
annotgene <- annotgene[annotgene$Chromosome %in% c(as.character(1:22) ,"X","Y"),]
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)

## Annotation... solving some issues...
#so when setting Entrezgene as rownames, the system does not allow them to have duplicate names and you get a warning message
#this can be used to detect duplicate entries.
rownames(annotgene) <- annotgene$Entrezgene
annotgene[annotgene$Entrezgene=="132989",]

##Filter the duplicate genes and double check there are no more duplicates
annotgene_filt <- annotgene[!duplicated(annotgene$Entrezgene),]
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)
sum(annotgene_filt$Entrezgene %in% rownames(GSE198256_count))
annotgene_filt[annotgene_filt$Entrezgene=="132989",]

## Overlap between annotation and gnes
#assign gene id as row names for the filtered df
rownames(annotgene_filt) <- as.character(annotgene_filt$Entrezgene)
sum(as.character(rownames(annotgene_filt)) %in% rownames(GSE198256_count))

##  Work with the annotated genes!
#for this code 3 df are created: 1) containing those genes in the count table that matches with the filtered annotated data
#and 2) the excluded genes from the count table and 3) contains ordered data in filtered annotation that match with the 
#contained in my filtered count table.
GSE198256_count_filt <- GSE198256_count[rownames(GSE198256_count) %in% rownames(annotgene_filt),]
GSE198256_count_exc <-GSE198256_count[!(rownames(GSE198256_count) %in% rownames(annotgene_filt)),]
annotgene_ord <- annotgene_filt[rownames(GSE198256_count_filt ),]

sum(rownames(annotgene_ord)==rownames(GSE198256_count_filt))


######
library(NOISeq)
#BiocManager::install("NOISeq",force = TRUE)

##this codes generate a df that contains all the the health states and rename the column of the generated df as "group"

Factors_GSE198256 <- data.frame(Meta_GSE198256 [ colnames(GSE198256_count_filt),c("disease state:ch1")])
colnames(Factors_GSE198256)[1]<- "Group"

data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=abs(annotgene_ord$end-annotgene_ord$start),
                        gc=annotgene_ord$GC,
                        biotype= annotgene_ord$type ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)


#not working because i have not assigned biotype, etc.
##myexplodata <- dat(data_NOISEQ, type = "biotype")
##explo.plot(myexplodata, plottype = "persample")
##mynicedata <- dat2save(myexplodata)
##mybiodetection <- dat(data_NOISEQ, k = 0, type = "biodetection", factor = NULL)


####
##assigning values to length, gc content and chr # to use NOISEQ
####
lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
names(lengthuse) <- rownames(annotgene_ord)
gc <- annotgene_ord$GC
names(gc) <- rownames(annotgene_ord)
biotype <-annotgene_ord$type
names(biotype) <- rownames(annotgene_ord)

chromosome <- annotgene_ord[,c("Chromosome","start","end")]
browseVignettes("NOISeq")

data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)

#to determine the presence of contaminant RNA
myexplodata <- dat(data_NOISEQ, type = "biodetection")
explo.plot(myexplodata, plottype = "persample")

#same as above but comparing 2 samples
par(mfrow = c(1, 2))
explo.plot(myexplodata, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")

#to check the number of counts per biotype and its distribution
mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")

## 1) saturation curve showing the behavior of sequence depth related to the discovery of new features
## 2) Comparison of newly discovered features between 4 differente samples
mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:2, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:4)

#shows the distribution of counts of the biotype protein coding
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")

#Visualize the features with low count per sample to decide a potential removal threshold
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

#explore if lenght of the features introduce bias to the analysis
mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

#explore if GC content of the features introduce bias to the analysis
myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")

#comparison of read distribution per feature of 12 samples against a reference 
mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd,samples = 1:12)

#create a PCA and general QC
myPCA = dat(data_NOISEQ, type = "PCA")
explo.plot(myPCA, factor = "Group")

QCreport(data_NOISEQ, samples = NULL, factor = "Group", norm = FALSE)

########
##### SAVE JUST IN CASE
########

save(data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="GSE198256_step1.Rda")

############
## STEP 3: NORMALIZATION & DIFF EXPRESSION
############
#normalization to correct for lenght of features

myRPKM = rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA = uqua(assayData(data_NOISEQ)$exprs, long = lengthuse, lc = 0.5, k = 0)
myTMM = tmm(assayData(data_NOISEQ)$exprs, long = 1000, lc = 0)

############
## STEP 3.1: DESEQ2
############

browseVignettes("DESeq2")
#BiocManager::install("DESeq2")
# https://mac.r-project.org/tools/
# sudo xcode-select --install
library(DESeq2)
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
setwd("~/Dropbox/WORKING_KAUST/LECTURING/B3XX_SettinngPipelines/DATA_Sets/GSE198256")

load("GSE198256_step1.Rda") # data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="GSE198256_step1.Rda")

############
# STEP 3.1.1: SET THE CLASS
############

GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pData(data_NOISEQ),
                                           design = ~ Group)

# Warning
#renaming the different groups (because of the previous warning) and defining data to work with
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="Covid19: Acute infection"] <- "Covid19AI"
pDataUSE[pDataUSE=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
pDataUSE[pDataUSE=="Covid19: Recovery 6Mo"] <- "Covid196Mo"
str(pDataUSE)
#this last line of code converts the first column (containing the groups) into a factors or categories
pDataUSE[,1] <- as.factor(pDataUSE[,1])



##DESeq2 is provided with the data. The design is specified, meaning the categories to have into account.
#a "-1 + group" means that the intercept is excluded. Including the intercept may be appropriate when you 
#want to compare each group to a baseline condition, while excluding the intercept may be suitable for 
#comparing groups directly without a specific reference level.

GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ -1 + Group)
resultsNames(GSE198256_DESeq2)
GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ Group)

str(GSE198256_DESeq2)


############
# STEP 3.1.2: WITH WHICH GENES TO WORK?
############

## Do we use all the genes?
## How do we select which ones?

#The code below filters the count matrix and keep only genes that have counts greater than or 
#equal to 10 in at least six samples
smallestGroupSize <- 6
keep <- rowSums(counts(GSE198256_DESeq2) >= 10) >= smallestGroupSize
#Once the filter criteria is ready, we can create a modified phenodata containing genes of interest.
GSE198256_DESeq2_F <- GSE198256_DESeq2[keep,]

#how many genes do i have after the filter
length(rownames(counts(GSE198256_DESeq2_F)))

############
# STEP 3.1.3: DIFFERENTIAL EXPRESSION?
############

GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
GSE198256_res <- results(GSE198256_DESeq2_F)
GSE198256_res
resultsNames(GSE198256_DESeq2_F)
str(GSE198256_DESeq2_F)


############
# STEP 3.1.4: WE NEED TO UNDERSTAND MORE...
############

## Questions in my mind:
# How do I define the question?
# How the differential expression is done?
# How to interpret the results?
# Technical replicates?

## STEP 3.1.4: plot MA

plotMA(GSE198256_res, ylim=c(-2,2))
#Interpretation?


lfcShrink(GSE198256_DESeq2_F,coef=c("Group_Healthy_vs_Covid193Mo"))
res_lfcShrink <- lfcShrink(GSE198256_DESeq2_F,coef=c("Group_Covid196Mo_vs_Covid193Mo"))

plotMA(res_lfcShrink, ylim=c(-2,2))
# Why to shrink: it looks at the largest fold changes that are not due 
# to low counts and uses these to inform a prior distribution. 
# So the large fold changes from genes with lots of statistical information are 
# not shrunk, while the imprecise fold changes are shrunk. This allows you to 
# compare all estimated LFC across experiments, for example, which is not really
# feasible without the use of a prior. 
# Michael Love https://support.bioconductor.org/p/77461/



## STEP 3.1.4: Define questions

GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
#res <- results(GSE198256_DESeq2_F, contrast=c('factorName','numeratorLevel','denominatorLevel'))
res <- results(GSE198256_DESeq2_F, contrast=c("Group","Covid19AI","Healthy"))
res
resultsNames(GSE198256_DESeq2_F)


## STEP 3.1.4: How differential expression is conducted...

# DESeq2 offers two kinds of hypothesis tests: 
#   the Wald test, 
#        where we use the estimated standard error of a log2 fold 
#        change to test if it is equal to zero, 
#   the likelihood ratio test (LRT). 
#        The LRT examines two models for the counts, a full model 
#        with a certain number of terms and a reduced model, in 
#        which some of the terms of the full model are removed. 
#        The test...



#####
#Code to compare the with the article
####
acutevshealthy <- results(GSE198256_DESeq2_F, contrast=c("Group","Covid19AI","Healthy"))
mo3vshealthy <- results(GSE198256_DESeq2_F, contrast=c("Group","Covid193Mo","Healthy"))
mo6vshealthy <- results(GSE198256_DESeq2_F, contrast=c("Group","Covid196Mo","Healthy"))

#plot acute vs control
plotMA(acutevshealthy, ylim=c(-2,2))

#plot early recovery vs control
plotMA(mo3vshealthy, ylim=c(-2,2))

#plot late recovery vs control
plotMA(mo6vshealthy, ylim=c(-2,2))
