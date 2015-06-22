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
expression_bins_no_trx <- 4
expression_bin_names_trx <- c(
	"Expression--col--_-2--dot--3_to_0--dot--0_Log2_TPM_--op--4--com--751--cl--",
	"Expression--col--_0--dot--0_to_3--dot--0_Log2_TPM_--op--4--com--751--cl--",
	"Expression--col--_3--dot--0_to_5--dot--5_Log2_TPM_--op--4--com--751--cl--",
	"Expression--col--_5--dot--5_to_13--dot--2_Log2_TPM_--op--4--com--751--cl--"
)
expression_filter_trx <- expression_bin_names_trx[2:3]
transcript_length_bins_no_trx <- 4
transcript_length_bin_names_trx <- c(
	"Transcript_length--col--_40_to_570_nts_--op--2--com--376--cl--",
	"Transcript_length--col--_570_to_840_nts_--op--2--com--376--cl--",
	"Transcript_length--col--_840_to_1990_nts_--op--2--com--376--cl--",
	"Transcript_length--col--_1990_to_29630_nts_--op--2--com--374--cl--"
)
exon_no_bins_max_trx <- c(1,2,3,4,6,10)
exon_no_bin_names_trx <- c(
	"Exons_per_transcript--col--_1_--op--788--cl--",
	"Exons_per_transcript--col--_2_--op--1--com--585--cl--",
	"Exons_per_transcript--col--_3_--op--1--com--359--cl--",
	"Exons_per_transcript--col--_4_--op--1--com--203--cl--",
	"Exons_per_transcript--col--_5_to_6_--op--1--com--747--cl--",
	"Exons_per_transcript--col--_7_to_10_--op--1--com--385--cl--",
	"Exons_per_transcript--col--_11_to_107_--op--1--com--435--cl--"
)
gc_content_bins_no_trx <- 4
gc_content_bin_names_trx <- c(
	"GC_percentage--col--_20_to_44_--op--2--com--376--cl--",
	"GC_percentage--col--_44_to_50_--op--2--com--376--cl--",
	"GC_percentage--col--_50_to_57_--op--2--com--376--cl--",
	"GC_percentage--col--_57_to_78_--op--2--com--374--cl--"
)

expression_bins_no_pas <- 4
expression_bin_names_pas <- c(
	"Expression--col--_-2--dot--3_to_0--dot--3_Log2_TPM_--op--1--com--100--cl--",
	"Expression--col--_0--dot--3_to_3--dot--4_Log2_TPM_--op--1--com--100--cl--",
	"Expression--col--_3--dot--4_to_5--dot--7_Log2_TPM_--op--1--com--100--cl--",
	"Expression--col--_5--dot--7_to_13--dot--2_Log2_TPM_--op--1--com--100--cl--"
)
expression_filter_pas <- expression_bin_names_pas[2:3]
trx_no_bins_max_pas <- c(1,3)
trx_no_bin_names_pas <- c(
	"Transcripts_per_poly--op--A--cl--site--col--_1_--op--759--cl--",
	"Transcripts_per_poly--op--A--cl--site--col--_2_to_3_--op--822--cl--",
	"Transcripts_per_poly--op--A--cl--site--col--_4_to_17_--op--619--cl--"
)

expression_bins_no_gen <- 4
expression_bin_names_gen <- c(
	"Expression--col--_-2--dot--3_to_1--dot--1_Log2_TPM_--op--3--com--232--cl--",
	"Expression--col--_1--dot--1_to_4--dot--1_Log2_TPM_--op--3--com--232--cl--",
	"Expression--col--_4--dot--1_to_6--dot--2_Log2_TPM_--op--3--com--232--cl--",
	"Expression--col--_6--dot--2_to_13--dot--2_Log2_TPM_--op--3--com--229--cl--"
)
expression_filter_gen <- expression_bin_names_gen[2:3]
trx_no_bins_max_gen <- c(1,4,8,12)
trx_no_bin_names_gen <- c(
	"Transcripts_per_gene--col--_1_--op--1--com--322--cl--",
	"Transcripts_per_gene--col--_2_to_4_--op--1--com--217--cl--",
	"Transcripts_per_gene--col--_5_to_8_--op--1--com--428--cl--",
	"Transcripts_per_gene--col--_9_to_12_--op--1--com--134--cl--",
	"Transcripts_per_gene--col--_13_to_59_--op--1--com--363--cl--"	
)
pas_no_bins_max_gen <- c(1,2)
pas_no_bin_names_gen <- c(
	"Poly--op--A--cl--_sites--col--_1_--op--2--com--611--cl--",
	"Poly--op--A--cl--_sites--col--_2_--op--1--com--113--cl--",
	"Poly--op--A--cl--_sites--col--_3_to_9_--op--530--cl--"
)


#<--- FUNCTIONS --->#

## Generate a fixed number of bins; sorted values are distributed over (approximately) equal-sized bins, from lowest to highest values
equal_bins <- function(data_vector, bin_number, names) {
        data_vector <- sort(data_vector)
        group_size <- ceiling(length(data_vector) / bin_number)
        bin_start <- seq(1, length(data_vector), group_size)
        bins <- lapply(bin_start, function(start) {
                end <- start + group_size - 1
                if (end > length(data_vector)) end <- length(data_vector)
                data_vector[start:end]
        })
        names(bins) <- names
        return(bins)
}

## Generate bins based on a vector of maximum allowed values per bin; values exceeding the last max value are added to a final bin; supply names for length(max_vector) + 1
discrete_bins <- function(data_vector, max_vector, names) {
        bins <- lapply(max_vector, function(max) {
                return <- data_vector[data_vector <= max]
                data_vector <<- data_vector[data_vector > max]
                return(return)
        })
        bins$last <- data_vector
        names(bins) <- names
        return(bins)
}


#<--- LOAD PACKAGES --->#

# Load 'rtracklayer' package
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }


#<--- LOAD ANNOTATIONS, LOOKUP TABLES & GC CONTENT INFO --->#

# Import GTF file
gtf <- import(file.path(root, "resources", "gencode.v19.annotation.ENS_compatible.gtf"), asRangedData=FALSE)

## Load lookup tables
trx_2_gene <- read.table(file.path(root, "resources", "gencode.v19.annotation.ENS_compatible.trx_gene_lookup_table"), col.names=c("trx_id", "gene_id"), colClasses=c("character", "character"), row.names=1)
trx_2_pas  <- read.table(file.path(root, "resources", "gencode.v19.annotation.ENS_compatible.trx_pas_lookup_table"), col.names=c("trx_id", "pas_id"), colClasses=c("character", "character"), row.names=1)
pas_2_gene <- read.table(file.path(root, "resources", "gencode.v19.annotation.ENS_compatible.pas_gene_lookup_table"), col.names=c("pas_id", "gene_id"), colClasses=c("character", "character"))

## Load GC content info
gc_content <- read.table(file.path(root, "resources", "Homo_sapiens.GRCh37.74.cdna_ncrna.GC_compatible.gc_content.tab"), header=TRUE)[,1:2]
gc_content <- setNames(gc_content$Percent_GC, gc_content$ID)


#<--- LOAD GROUND TRUTH TRANSCRIPT ABUNDANCES --->#

## Load true expression data into named, name-sorted vector
trx_df <- read.table(file.path(root, "estimates", "sim_1.transcripts.Ground_truth.reference"), col.names=c("trx_id", "value"), colClasses=c("character", "numeric"), row.names=1)
trx_vec <- setNames(trx_df$value, rownames(trx_df))
trx_vec <- trx_vec[order(names(trx_vec))]


#<--- INITIALIZE CONTAINER FOR SUBSETS --->#

# Generate container for expression subsets
subsets <- list()


#<--- TRANSCRIPT SUBSETS --->#

# Initialize subsets transcript container
subsets$trx <- list()

        #<-- Expression levels -->#

        ## Generate bins
        not_expressed <- trx_vec[trx_vec == 0]
        expressed <- trx_vec[trx_vec != 0]
	expression_bins <- equal_bins(expressed, expression_bins_no_trx, expression_bin_names_trx)

# Apply expression level filter for further subsetting
expression_level_filter <- setNames(unlist(expression_bins[expression_filter_trx], use.names=FALSE), unlist(lapply(expression_bins[expression_filter_trx], names), use.names=FALSE))

        #<--- Transcript lengths --->#

        # Subset only exons
        exons <- gtf[gtf$type == "exon"]

        # Filter only transcripts in the specified expression range
        exons_filtered <- exons[exons$transcript_id %in% names(expression_level_filter)]

        ## Get transcript lengths
        exon_lengths <- data.frame(transcript_id=exons_filtered$transcript_id, width=width(exons_filtered))
        trx_lengths <- aggregate(exon_lengths$width, by=list(exon_lengths$transcript_id), sum)
        trx_lengths <- setNames(trx_lengths[,2], trx_lengths[,1])
        transcript_length_bins <- equal_bins(trx_lengths, transcript_length_bins_no_trx, transcript_length_bin_names_trx)
        
        #<--- Number of exons --->#

        ## Get number of exons per transcript
        exon_no <- table(exons_filtered$transcript_id)
	exon_no_bins <- discrete_bins(exon_no, exon_no_bins_max_trx, exon_no_bin_names_trx)

        #<--- GC content --->#
        gc_content_filtered <- gc_content[names(gc_content) %in% names(expression_level_filter)]
	gc_content_bins <- equal_bins(gc_content_filtered, gc_content_bins_no_trx, gc_content_bin_names_trx)

## Add bins to subsets
subsets$trx$all <- names(trx_vec)
subsets$trx$`All_transcripts_--op--187--com--176--cl--` <- names(trx_vec)
subsets$trx$`Not_expressed_--op--168--com--172--cl--` <- names(not_expressed)
subsets$trx$`Expressed_--op--19--com--004--cl--` <- names(expressed)
subsets$trx <- c(subsets$trx, lapply(expression_bins, names))
subsets$trx <- c(subsets$trx, lapply(transcript_length_bins, names))
subsets$trx <- c(subsets$trx, lapply(exon_no_bins, names))
subsets$trx <- c(subsets$trx, lapply(gc_content_bins, names))


#<--- PAS SUBSETS --->#

# Initialize subsets transcript container
subsets$pas <- list()

# Get poly(A) sites from transcripts
pas_vec <- merge(trx_2_pas, trx_df, by=0)
pas_vec <- aggregate(value~pas_id, pas_vec, sum)
pas_vec <- setNames(pas_vec$value, pas_vec$pas_id)

        #<-- Expression levels -->#

        ## Generate bins
        not_expressed <- pas_vec[pas_vec == 0]
        expressed <- pas_vec[pas_vec != 0]
        expression_bins <- equal_bins(expressed, expression_bins_no_pas, expression_bin_names_pas)

# Apply expression level filter for further subsetting
expression_level_filter <- setNames(unlist(expression_bins[expression_filter_pas], use.names=FALSE), unlist(lapply(expression_bins[expression_filter_pas], names), use.names=FALSE))	

        #<--- Number of transcripts --->#

        # Get transcripts per PAS
        trx <- table(trx_2_pas$pas_id)

        # Filter only PAS in the specified expression range
        trx_filtered <- trx[names(trx) %in% names(expression_level_filter)]

        ## Get number of transcripts per gene
        trx_no <- sort(trx_filtered)
	trx_no_bins <- discrete_bins(trx_no, trx_no_bins_max_pas, trx_no_bin_names_pas)
		
## Add bins to subsets
subsets$pas$all <- names(pas_vec)
subsets$pas$`All_poly--op--A--cl--sites_--op--25--com--139--cl--` <- names(pas_vec)
subsets$pas$`Not_expressed_--op--20--com--739--cl--` <- names(not_expressed)
subsets$pas$`Expressed_--op--4--com--400--cl--` <- names(expressed)
subsets$pas <- c(subsets$pas, lapply(expression_bins, names))
subsets$pas <- c(subsets$pas, lapply(trx_no_bins, names))


#<--- GENE SUBSETS --->#

# Initialize subsets transcript container
subsets$gene <- list()

# Get genes from transcripts
gen_vec <- merge(trx_2_gene, trx_df, by=0)
gen_vec <- aggregate(value~gene_id, gen_vec, sum)
gen_vec <- setNames(gen_vec$value, gen_vec$gene_id)

        #<-- Expression levels -->#

        ## Generate bins
        not_expressed <- gen_vec[gen_vec == 0]
        expressed <- gen_vec[gen_vec != 0]
        expression_bins <- equal_bins(expressed, expression_bins_no_gen, expression_bin_names_gen)

# Apply expression level filter for further subsetting
expression_level_filter <- setNames(unlist(expression_bins[expression_filter_gen], use.names=FALSE), unlist(lapply(expression_bins[expression_filter_gen], names), use.names=FALSE))

        #<--- Number of transcripts --->#

        # Get transcripts per gene
        trx <- table(trx_2_gene$gene_id)

        # Filter only genes in the specified expression range
        trx_filtered <- trx[names(trx) %in% names(expression_level_filter)]

        ## Get number of transcripts per gene
        trx_no <- sort(trx_filtered)
	trx_no_bins <- discrete_bins(trx_no, trx_no_bins_max_gen, trx_no_bin_names_gen)
	
        #<--- Number of poly(A) sites --->#

        # Get PAS per gene
        pas <- table(pas_2_gene$gene_id)

        # Filter only genes in the specified expression range
        pas_filtered <- pas[names(pas) %in% names(expression_level_filter)]

        ## Get number of PAS per gene
        pas_no <- sort(pas_filtered)
	pas_no_bins <- discrete_bins(pas_no, pas_no_bins_max_gen, pas_no_bin_names_gen)

## Add bins to subsets
subsets$gene$all <- names(gen_vec)
subsets$gene$`All_genes_--op--48--com--755--cl--` <- names(gen_vec)
subsets$gene$`Not_expressed_--op--35--com--830--cl--` <- names(not_expressed)
subsets$gene$`Expressed_--op--12--com--925--cl--` <- names(expressed)
subsets$gene <- c(subsets$gene, lapply(expression_bins, names))
subsets$gene <- c(subsets$gene, lapply(trx_no_bins, names))
subsets$gene <- c(subsets$gene, lapply(pas_no_bins, names))


#<--- SAVE SUBSETS --->#
save(subsets, file=outFile)
