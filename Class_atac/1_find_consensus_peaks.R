# Script for finding consensus peaks

# set up
rm(list = ls())
library(rtracklayer)
library(GenomicRanges)
setwd("C:/Users/danie/OneDrive/Documents/bioinformatic_pipelines/class_3/data_code_students/bed_bigwig_files")
# control panel
raw_data_folder <- "C:/Users/danie/OneDrive/Documents/bioinformatic_pipelines/class_3/data_code_students/bed_bigwig_files"
min_overlap <- 0.3
chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))

# bed files 
bed_files <- file.path(raw_data_folder, 
                       dir(raw_data_folder, pattern = '*.bed'))
bed_files <- bed_files[1:6]
bed_files
import

# load bed files
bed_1 <- import(bed_files[1])

# inspecting
bed_1

# useful method for GenomicRanges
length(bed_1)
as.character(seqnames(bed_1))[1:5]
seqlevels(bed_1)
start(bed_1)[1:5]
end(bed_1)[1:5]
width(bed_1)[1:5]

# appending information to bed_1
bed_1$num_bps
bed_1$num_bps <- width(bed_1)

# subsetting bed_1 (Filters all the peaks that are not annotated to according to chromosomes list provided above)
bed_1 <- bed_1[seqnames(bed_1) %in% chromosomes]
seqlevels(bed_1) <- chromosomes
length(bed_1)

# loading and subsetting the second and third bed file
bed_2 <- import(bed_files[2])
bed_2 <- bed_2[seqnames(bed_2) %in% chromosomes]
seqlevels(bed_2) <- chromosomes
length(bed_2)
#bed_3 <- import('bed_bigwig_files/GSM2083756_HL60-Rep3.peaks.bed')
#bed_3 <- bed_3[seqnames(bed_3) %in% chromosomes]
#seqlevels(bed_3) <- chromosomes
#length(bed_3)

# finding overlap between the two first files (It happens that peaks does not overlap between samples in some cases.
# even when they are replicates)
hits <- findOverlaps(bed_1, bed_2)

# inspecting hits...
hits

# what are the overlaps with at least min_overalp for both peaks?

# quantify the overlap (This creates basically an index, where the peaks from the first sample are matched with 
# a second sample)
overlaps <- pintersect(bed_1[queryHits(hits)], bed_2[subjectHits(hits)])
overlaps

# overlap fraction with respect to the original peaks
percent_overlap_on_1 <- width(overlaps) / width(bed_1[queryHits(hits)])
percent_overlap_on_2 <- width(overlaps) / width(bed_2[subjectHits(hits)])
hits <- hits[percent_overlap_on_1 > min_overlap & 
               percent_overlap_on_2 > min_overlap]



# Lets check the number of overlaping peaks
length(hits)

# subsetting the bed files
bed_1 <- bed_1[queryHits(hits)]
head(bed_1)
length(bed_1)
bed_2 <- bed_2[subjectHits(hits)]
head(bed_2)
length(bed_2)

# "reducing" the peaks
start_1 <- start(bed_1)
end_1 <- end(bed_1)
start_2 <- start(bed_2)
end_2 <- end(bed_2)
start_3 <- start(bed_3)
end_3 <- end(bed_3)
reduced_start <- pmin(start_1, start_2, start_3)
reduced_end <- pmax(end_1, end_2, start_3)
reference_bed <- bed_1
start(reference_bed) <- reduced_start
end(reference_bed) <- reduced_end
reference_bed
length(reference_bed)
#### apply for all bed files ####

# your code here!

for (i in 2:length(bed_files)) {
  # Load BED file
  bed_data <- import(bed_files[i])
  
  # Processing data
  bed_data <- bed_data[seqnames(bed_data) %in% chromosomes]
  seqlevels(bed_data) <- chromosomes
  
  # finding peak overlaps
  hits <- findOverlaps(bed_1, bed_data)
  overlaps <- pintersect(bed_1[queryHits(hits)], bed_data[subjectHits(hits)])
  
  # overlap fraction with respect to the original peaks
  percent_overlap_on_1 <- width(overlaps) / width(bed_1[queryHits(hits)])
  percent_overlap_on_2 <- width(overlaps) / width(bed_data[subjectHits(hits)])
  hits <- hits[percent_overlap_on_1 > min_overlap & 
                 percent_overlap_on_2 > min_overlap]
  length(hits)
  
  # subsetting bed files
  bed_data <- bed_data[subjectHits(hits)]
  head(bed_data)
  length(bed_data)
  
  
  # Save processed BED data in variables
  assign(paste0("bed_", i), bed_data)
}




reference_bed <- bed_1














# reference_bed <- ...


# how many peaks?
length(reference_bed)

# width?
summary(width(reference_bed))
summary(width(bed_1))

#### Discussion: when is a peak present? ####

# overlapping peaks in HL60

# your code here
HL60_reference_bed <- reference_bed
  
# how many peaks?
length(HL60_reference_bed)

# overlapping peaks in 3h-Mac

# your code here
start_1 <- start(bed_4)
end_1 <- end(bed_4)
start_2 <- start(bed_5)
end_2 <- end(bed_5)
start_3 <- start(bed_6)
end_3 <- end(bed_6)
reduced_start <- pmin(start_1, start_2, start_3)
reduced_end <- pmax(end_1, end_2, start_3)
Mac_reference_bed <- bed_1
start(Mac_reference_bed) <- reduced_start
end(Mac_reference_bed) <- reduced_end
Mac_reference_bed
length(Mac_reference_bed)


# how many peaks?
length(Mac_reference_bed)

# peak union!
union_reference_bed <- c(HL60_reference_bed, Mac_reference_bed)
length(union_reference_bed)

# reducing: concatenating intervals that are overlapping
union_reference_bed <- reduce(union_reference_bed)
length(union_reference_bed)

# loading black listed regions
# https://www.encodeproject.org/annotations/ENCSR636HFF/
black_listed_bed <- import('ENCFF356LFX.bed')

# any hit?
hits <- findOverlaps(union_reference_bed, black_listed_bed)
hits

# what about the length of the overlap?
overlaps <- pintersect(union_reference_bed[queryHits(hits)], 
                       black_listed_bed[subjectHits(hits)])
summary(width(overlaps)/width(union_reference_bed[queryHits(hits)]))

# eliminating the blacklisted regions
union_reference_bed <- union_reference_bed[-queryHits(hits)]

# writing the reference
export.bed(union_reference_bed, con = 'refererence.bed')
