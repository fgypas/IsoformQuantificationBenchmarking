#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 25, 2014
### Modified: May 28, 2015
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
#==================#
#    HEADER END    #
#==================#


#<--- COMMAND-LINE ARGUMENTS --->

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
outFile <- args[2]


#<--- PARAMETERS --->#

param_trx_biotype_filter <- c("antisense", "processed_pseudogene", "lincRNA", "nonsense_mediated_decay", "retained_intron", "processed_transcript", "protein_coding")
param_pseudocount <- 1/32
param_pas_expr_bins <- 4
param_trx_biotype_bins <- c("protein_coding", "retained_intron", "processed_transcript", "nonsense_mediated_decay", "antisense", "lincRNA")
param_gene_expr_bins <- 4
param_gene_biotype_bins <- c("protein_coding", "lincRNA", "antisense")
param_gene_expr_filter <- 2:3
param_gene_trx_no_bins_max <- c(1, 3, 6)
param_gene_pas_no_bins_max <- c(1, 2)


#<--- FUNCTIONS --->#

## Generate a fixed number of bins; sorted values are distributed over (approximately) equal-sized bins, from lowest to highest values
equal_bins <- function(data_vector, bin_number) {
        data_vector <- sort(data_vector)
        group_size <- ceiling(length(data_vector) / bin_number)
        bin_start <- seq(1, length(data_vector), group_size)
        bins <- lapply(bin_start, function(start) {
                end <- start + group_size - 1
                if (end > length(data_vector)) end <- length(data_vector)
                data_vector[start:end]
        })
        return(bins)
}

## Generate bins based on a vector of maximum allowed values per bin; values exceeding the last max value are added to a final bin
discrete_bins <- function(data_vector, max_vector) {
	data_vector <- sort(data_vector)
        bins <- lapply(max_vector, function(max) {
                return <- data_vector[data_vector <= max]
                data_vector <<- data_vector[data_vector > max]
                return(return)
        })
        bins$last <- data_vector
        return(bins)
}

## Generates pretty names for gene/transcript biotypes
prettyBiotypes <- function(character, prefix) {
	character <- gsub("protein_coding", "protein-coding", character)
	gsub("_", " ", character)
	paste(prefix, "biotype:", character)
}

## Generates pretty names for expression bins
prettyExpressionBins <- function(list) {
	min <- round(sapply(lapply(list, log2), min), digits=1)
	max <- round(sapply(lapply(list, log2), max), digits=1)
	paste("Expression:", min, "to", max, "log2 TPM")
}

## Generates pretty names for discrete bins
prettyDiscreteBins <- function(list, prefix) {
	prefix <- paste0(prefix, ":")
	min <- round(sapply(list, min))
	max <- round(sapply(list, max))
	ifelse(min == max, paste(prefix, min), paste(prefix, min, "to", max))
}

## Adds lengths to filter names
addLengthToNames <- function(list, names) {
	lengths <- sapply(list, length)
	lengths <- format(lengths, big.mark=",", trim=TRUE, scientific=FALSE)
	paste0(names, " (", lengths, ")")
}

## Make names file system friendly
friendlyFilenames <- function(list) {
	names(list) <- gsub("(", "--op--", names(list), fixed=TRUE)
	names(list) <- gsub(")", "--cl--", names(list), fixed=TRUE)
	names(list) <- gsub(",", "--com--", names(list), fixed=TRUE)
        names(list) <- gsub(".", "--dot--", names(list), fixed=TRUE)
	names(list) <- gsub(":", "--col--", names(list), fixed=TRUE)
	gsub(" ", "_", names(list), fixed=TRUE)
}


#<--- LOAD PACKAGES --->#

# Load 'rtracklayer' package
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }


#<--- LOAD ANNOTATIONS, LOOKUP TABLES & GC CONTENT INFO --->#

# Import GTF file
gtf <- import(file.path(root, "resources", "gencode.vM2.annotation.ENS_compatible.gtf"), asRangedData=FALSE)

## Load lookup tables
lookup_trx_2_gene <- read.table(file.path(root, "resources", "gencode.vM2.annotation.ENS_compatible.filtered.trx_gene_lookup_table"), col.names=c("trx_id", "gene_id"), colClasses=c("character", "character"))
lookup_trx_2_pas  <- read.table(file.path(root, "resources", "gencode.vM2.annotation.ENS_compatible.filtered.trx_pas_lookup_table"), col.names=c("trx_id", "pas_id"), colClasses=c("character", "character"))
lookup_pas_2_gene <- read.table(file.path(root, "resources", "gencode.vM2.annotation.ENS_compatible.filtered.pas_gene_lookup_table"), col.names=c("pas_id", "gene_id"), colClasses=c("character", "character"))

# Subset & filter transcripts
trx <- gtf[gtf$type == "transcript" & gtf$transcript_type %in% param_trx_biotype_filter]

# Subset & filter genes
genes <- gtf[gtf$type == "gene" & gtf$gene_id %in% lookup_trx_2_gene$gene_id]


#<--- LOAD / EXTRACT GROUND TRUTH PAS ABUNDANCES --->#

## Load true expression data into named, name-sorted vector
pas_df_1 <- read.table(file.path(root, "estimates", "mmu_1.processing_sites.A-seq-2.reference"), col.names=c("pas_id", "value"), colClasses=c("character", "numeric"), row.names=1)
pas_df_2 <- read.table(file.path(root, "estimates", "mmu_2.processing_sites.A-seq-2.reference"), col.names=c("pas_id", "value"), colClasses=c("character", "numeric"), row.names=1)
pas_df <- merge(pas_df_1, pas_df_2, by=0)
pas_vec <- setNames((pas_df$value.x + pas_df$value.y) / 2, pas_df$Row.names)

# Get all unique poly(A) sites that overlap with transcript ends
pas_with_trx_all <- sort(unique(lookup_trx_2_pas$pas_id))

# Subset poly(A) sites that have an annotated transcript
pas_vec <- pas_vec[names(pas_vec) %in% pas_with_trx_all]

## Add non-expressed (according to A-seq) poly(A) sites that have an annotated transcript (set param_pseudocount)
add <- pas_with_trx_all[! pas_with_trx_all %in% names(pas_vec)]
pas_vec <- c(pas_vec, setNames(rep(param_pseudocount, length(add)), add))
pas_vec <- pas_vec[order(names(pas_vec))]


#<--- EXTRACT GROUND TRUTH GENE ABUNDANCES --->#

## Sum poly(A) site abundances for each gene
gene_vec <- merge(lookup_pas_2_gene, pas_vec, by.x=1, by.y=0)
gene_vec <- aggregate(y~gene_id, gene_vec, sum)
gene_vec <- setNames(gene_vec$y, gene_vec$gene_id)


#<--- INITIALIZE CONTAINER FOR SUBSETS --->#

# Generate container for expression subsets
subsets <- list()


#<--- Poly(A) sites --->#

## Basic
pas_all <- pas_vec[names(pas_vec) %in% pas_with_trx_all]
pas_expr <- pas_vec[log2(pas_vec) > log2(param_pseudocount)]
pas_not_expr <- pas_all[setdiff(names(pas_all), names(pas_expr))]

## Transcript number
pas_all_trx_no <- table(lookup_trx_2_pas[lookup_trx_2_pas$pas_id %in% names(pas_all), 2])
pas_expr_trx_no <- table(lookup_trx_2_pas[lookup_trx_2_pas$pas_id %in% names(pas_expr), 2])
pas_all_one_trx <- names(pas_all_trx_no[pas_all_trx_no == 1])
pas_expr_one_trx <- names(pas_expr_trx_no[pas_expr_trx_no == 1])

## Expression bins
pas_expr_bins <- equal_bins(pas_expr, param_pas_expr_bins)
names(pas_expr_bins) <- prettyExpressionBins(pas_expr_bins)

## Add subsets to container
subsets$pas$all <- names(pas_all)
subsets$pas$`All poly(A) sites overlapping annotated transcript ends` <- names(pas_all)
subsets$pas$`Not expressed` <- names(pas_not_expr)
subsets$pas$`Expressed poly(A) sites overlapping annotated transcript ends` <- names(pas_expr)
subsets$pas$`All poly(A) sites overlapping exactly one annotated transcript end` <- pas_all_one_trx
subsets$pas$`Expressed poly(A) sites overlapping exactly one annotated transcript end` <- pas_expr_one_trx
subsets$pas <- c(subsets$pas, lapply(pas_expr_bins, names))

## Make pretty names
names(subsets$pas)[2:length(subsets$pas)] <- addLengthToNames(subsets$pas[2:length(subsets$pas)], names(subsets$pas)[2:length(subsets$pas)])

## Make pretty names
names(subsets$pas)[2:length(subsets$pas)] <- addLengthToNames(subsets$pas[2:length(subsets$pas)], names(subsets$pas)[2:length(subsets$pas)])


#<--- Transcripts --->#

## Basic
trx_ids_all <- sort(unique(trx$transcript_id))
trx_ids_anno <- sort(unique(lookup_trx_2_pas$trx_id))
trx_ids_expr <- sort(unique(lookup_trx_2_pas[lookup_trx_2_pas$pas_id %in% names(pas_expr), ]$trx_id))
trx_ids_not_expr <- setdiff(trx_ids_all, trx_ids_expr)

## Transcript biotypes
trx_ids_biotypes <- sapply(param_trx_biotype_bins, function(biotype) {
	trx[trx$transcript_type == biotype & trx$transcript_id %in% trx_ids_expr]$transcript_id
})
names(trx_ids_biotypes) <- prettyBiotypes(param_trx_biotype_bins, "Transcript")

## Add subsets to container
subsets$trx$all <- trx_ids_all
subsets$trx$`All annotated transcripts` <- trx_ids_all
subsets$trx$`Not expressed` <- trx_ids_not_expr
subsets$trx$`Transcripts ending in annotated poly(A) sites` <- trx_ids_anno
subsets$trx$`Transcripts ending in expressed poly(A) sites` <- trx_ids_expr
subsets$trx <- c(subsets$trx, trx_ids_biotypes)

## Make pretty names
names(subsets$trx)[2:length(subsets$trx)] <- addLengthToNames(subsets$trx[2:length(subsets$trx)], names(subsets$trx)[2:length(subsets$trx)])


#<--- Genes --->#

## Basic
gene_ids_all <- sort(unique(lookup_trx_2_gene$gene_id))
gene_ids_all_with_pas <- sort(unique(lookup_pas_2_gene[lookup_pas_2_gene$pas_id %in% names(pas_all), ]$gene_id))
gene_ids_expr <- sort(unique(lookup_pas_2_gene[lookup_pas_2_gene$pas_id %in% names(pas_expr), ]$gene_id))
gene_ids_not_expr <- setdiff(gene_ids_all, gene_ids_expr)
gene_ids_pas_not_expr <- setdiff(gene_ids_all_with_pas, gene_ids_expr)

# Expression bins
gene_expr <- gene_vec[names(gene_vec) %in% gene_ids_expr]
gene_expr_bins <- equal_bins(gene_expr, param_gene_expr_bins)
names(gene_expr_bins) <- prettyExpressionBins(gene_expr_bins)

## Gene biotypes
gene_biotype_bins <- sapply(param_gene_biotype_bins, function(biotype) {
	genes[genes$gene_type == biotype & genes$gene_id %in% gene_ids_expr]$gene_id
})
names(gene_biotype_bins) <- prettyBiotypes(param_gene_biotype_bins, "Gene")

# Apply expression level filter for further subsetting
gene_expr_filter <- setNames(unlist(gene_expr_bins[param_gene_expr_filter], use.names=FALSE), unlist(lapply(gene_expr_bins[param_gene_expr_filter], names), use.names=FALSE))

## Transcript number
gene_trx_no <- table(lookup_trx_2_gene$gene_id)
gene_trx_no <- gene_trx_no[names(gene_trx_no) %in% names(gene_expr_filter)]
gene_trx_no_bins <- discrete_bins(gene_trx_no, param_gene_trx_no_bins_max)
names(gene_trx_no_bins) <- prettyDiscreteBins(gene_trx_no_bins, "Transcripts per gene")

## PAS number
gene_pas_no <- table(lookup_pas_2_gene$gene_id)
gene_pas_no <- gene_pas_no[names(gene_pas_no) %in% names(gene_expr_filter)]
gene_pas_no_bins <- discrete_bins(gene_pas_no, param_gene_pas_no_bins_max)
names(gene_pas_no_bins) <- prettyDiscreteBins(gene_pas_no_bins, "Poly(A) sites per gene")

## Add subsets to container
subsets$gene$all <- gene_ids_all
subsets$gene$`All annotated genes` <- gene_ids_all
subsets$gene$`Not expressed` <- gene_ids_not_expr
subsets$gene$`Genes containing annotated poly(A) sites` <- gene_ids_all_with_pas
subsets$gene$`Genes containing annotated but not expressed poly(A) sites` <- gene_ids_pas_not_expr
subsets$gene$`Genes containing expressed poly(A) sites` <- gene_ids_expr
subsets$gene <- c(subsets$gene, lapply(gene_expr_bins, names))
subsets$gene <- c(subsets$gene, gene_biotype_bins)
subsets$gene <- c(subsets$gene, lapply(gene_trx_no_bins, names))
subsets$gene <- c(subsets$gene, lapply(gene_pas_no_bins, names))

## Make pretty names
names(subsets$gene)[2:length(subsets$gene)] <- addLengthToNames(subsets$gene[2:length(subsets$gene)], names(subsets$gene)[2:length(subsets$gene)])


#<--- MAKE FILE SYSTEM FRIENDLY NAMES --->#
for (level in names(subsets)) {
	names(subsets[[level]]) <- friendlyFilenames(subsets[[level]])
}


#<--- SAVE SUBSETS --->#
save(subsets, file=outFile)
