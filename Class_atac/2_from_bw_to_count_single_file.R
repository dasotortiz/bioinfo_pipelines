# Script for coverting bigWig + bed files to counts

# Partially adapted from https://lcolladotor.github.io/protocols/bigwig_DEanalysis/

# set up
rm(list = ls())
library(rtracklayer)
library(GenomicRanges)

# control panel
raw_data_folder <- 'bed_bigwig_files'
read_length <- 36
chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))
reference_file <- 'refererence.bed'

# bw files
bw_files <- file.path(raw_data_folder, dir(raw_data_folder, pattern = '*.bw'))
bw_files <- bw_files[1:6]

# loading reference bed
peaks <- import(reference_file)

# current file
i <- 1
print(paste0('sample ', i, ' out of ', length(bw_files)))
bw_file <- bw_files[i]

# loading and downsizing the bigwigfile
bw_file_list <- BigWigFileList(bw_file)
coverage <- import(bw_file_list[[1]], as = 'RleList')
coverage <- coverage[names(coverage) %in% chromosomes]
coverage[1]

# split the peaks across chromosomes
peaks_list <- split(peaks, seqnames(peaks))
peaks_list

# coverage per peak
coverage <- coverage[names(peaks_list)]
peaks_coverage <- Views(coverage, ranges(peaks_list))
peaks_coverage

# count values
counts <- sapply(peaks_coverage, sum)

# ensuring to have the right peak information
chrs <- rep(names(peaks_coverage), sapply(peaks_coverage, length))
starts <- sapply(peaks_coverage, start)
ends <- sapply(peaks_coverage, end)

# converting to vector
counts <- unlist(counts)
names(counts) <- paste0(chrs, '_', unlist(starts), '_', unlist(ends))
head(counts)

# rounding up
counts <- round(counts / read_length)
head(counts)

# fractions of reads in peaks
sum(counts)
all_counts <- sapply(coverage, sum)
all_counts <- sum(all_counts)
all_counts <- round(all_counts / read_length)
all_counts
sum(counts) / all_counts
