# Script for identifying TFs putative locations

# Adapted from: 
# https://ivanek.github.io/analysisOfGenomicsDataWithR/13_Motifs_html.html

# Set up 
rm(list = ls())
library(universalmotif)
library(JASPAR2020)
library(TFBSTools)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(glmnet)

# loading binding sites for GATA2
hGATA2_motif_seqs <- readRDS('hGATA2_binding_sites.rds')

# creating a motif
hGATA2_motif <- create_motif(hGATA2_motif_seqs, type="PCM", name = "GATA2")
hGATA2_motif

# converting to Position Probability Matrix
convert_type(hGATA2_motif, type="PPM", pseudocount = 1)

# converting to Position Weight Matrix
convert_type(hGATA2_motif, type="PWM", pseudocount = 1) # log2(PPM/bg)

# visualization
view_motifs(hGATA2_motif, use.type = "ICM")

# converting to Position Weight Matrix
hGATA2_motif <- convert_type(hGATA2_motif, type="PWM", pseudocount = 1)
hGATA2_motif['motif']

# how is PWM used?
test_seq <- "CAAGAT"
score_match(hGATA2_motif, test_seq)
# 0.8875253 + (-4.321928) + 1.944858 + (-4.321928) + 1.944858 + (-0.6214884)

# extracting TFs from JASPAR
opts <- list(
  tax_group    = "vertebrates",
  collection   = "CORE",
  matrixtype   = "PWM",
  all_versions = FALSE)

# run query
jaspar_pwms <- getMatrixSet(JASPAR2020, opts) 
length(jaspar_pwms)

# keeping only motifs that are longer than 10 bps
tmp <- sapply(jaspar_pwms, length)
hist(tmp, breaks=64)
jaspar_pwms <- jaspar_pwms[tmp > 10]
length(jaspar_pwms)

# mapping vector of motifs ID -> name
motif2tf <- sapply(jaspar_pwms, function(pwm) {pwm@name})
head(motif2tf)

# loading peaks
peaks <- import('refererence.bed')

# considering only chr2 peaks
peaks <- peaks[seqnames(peaks) == 'chr2']
seqlevels(peaks) <- 'chr2'

# computing motif scores
system.time(
  matches_rse <- matchMotifs(
    jaspar_pwms, peaks, p.cutoff = 5e-05,
    genome = BSgenome.Hsapiens.UCSC.hg38, 
    out="scores"
  ))

# inspecting
matches_rse
motifCounts(matches_rse)[1:3, 1:3]
rowRanges(matches_rse)[1:3]
motif2tf['MA0029.1']

# question: how can we use motif location with the results
# of the differential analysis?

# load differential analysis results
da_res <- read.csv('differential_analysis_res.csv')
rownames(da_res) <- da_res$id
da_res$id <- NULL
head(da_res)

# extracting the motif counts
motif_counts <- as.matrix(motifCounts(matches_rse))
rownames(motif_counts) <- paste0(seqnames(rowRanges(matches_rse)), '_', 
                                 start(rowRanges(matches_rse)), '_', 
                                 end(rowRanges(matches_rse)))
motif_counts <- scale(motif_counts, center = TRUE, scale = TRUE)
motif_counts[1:3, 1:3]

# ensuring da results are ordered in the same way as motif_counts
da_res <- da_res[rownames(motif_counts), ]
da_res[1:3, 1:3]

# let's now compute the impact of the motif presence on the log2FC!

# first we find the best hyper parameter...
cv.res <- cv.glmnet(x = motif_counts, 
                    y = da_res$log2FoldChange, 
                    family = gaussian, alpha = 0.99)
cv.res$lambda.min

# and here is the model
m <- glmnet(x = motif_counts, 
            y = da_res$log2FoldChange, 
            family = gaussian, alpha = 0.99, 
            lambda = cv.res$lambda.min)

# extracting the coefficients
coef(m)
m_coef <- as.numeric(coef(m)) # as vector
names(m_coef) <- rownames(coef(m)) # motif names
m_coef <- m_coef[abs(m_coef) > 0] # keeping only coeffs different from zero
names(m_coef) <- motif2tf[names(m_coef)]
m_coef

# positive coefficients: chromating open in Mac at 3h
head(sort(m_coef, decreasing = TRUE), 10)

# negative coefficients: chromating open in HL60
# the na is the intercept
tail(sort(m_coef, decreasing = TRUE), 10)

