log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("dplyr")
library("purrr")
library("tibble")


#Read in the first kraken file to get the tax_ids of bacteria...kraken must have been run with  "--report-zero-counts"
kraken_data <- read.table(snakemake@input[["kraken"]][1], header=FALSE, sep = "\t",quote = "", comment.char="")

#Read in the bracken data
bracken_data <- lapply(snakemake@input[["bracken"]], read.table, header=TRUE, sep = "\t", quote = "",comment.char="")


#This function will get the TaxIDs of all nodes under bacteria
get_bac_taxIds <- function(x){
    # X is a kraken report
    x<-x[,4:6]
    colnames(x) <- c("Rank_code", "Tax_ID", "ID")
    x$ID<- trimws(x$ID, which = "left") 
    domain_indices <- subset(x, Rank_code == "D") #get indexes of domain ranks 
    domain_indices$index <-c(1:nrow(domain_indices)) #make a vector for indexing
    bac_index <- domain_indices$index[domain_indices$ID =="Bacteria"] #Find which row is bacteria 
    next_domain <- bac_index+1
    bac_index <- as.numeric(row.names(domain_indices)[bac_index])
    if(next_domain == nrow(domain_indices)){
      stop_index <- nrow(x)
    }
    else{ 
      stop_index <-as.numeric(row.names(domain_indices)[next_domain])-1
    }
    x <- x[bac_index:stop_index,]
    return(x$Tax_ID)
}

#Get a vector of bacterial tax_ids 
bac_tax_ids <- get_bac_taxIds(kraken_data)

#Filter only the bacterial taxa in bracken_data
bracken_data_bacteria <- lapply(bracken_data, subset, taxonomy_id %in% bac_tax_ids)

#Recalculate the fraction of total reads. Extract name and fraction
bracken_bac_abundance <- lapply(bracken_data_bacteria, function(x){
    total_reads <- sum(x$new_est_reads)
    x$fraction_total_reads <- x$new_est_reads/total_reads
    out <- x[,c("name", "fraction_total_reads")]
    return(out)
})

bracken_bac_counts <- lapply(bracken_data_bacteria, function(x){
    total_reads <- sum(x$new_est_reads)
    x$fraction_total_reads <- x$new_est_reads/total_reads
    out <- x[,c("name", "new_est_reads")]
    return(out)
})


#Make Abundance Table 
# abund_table <- 
#     bracken_bac_abundance %>%
#         reduce(left_join, by = "name") %>%
#         column_to_rownames(var="name")

abund_table <- 
    bracken_bac_abundance %>%
        reduce(left_join, by = "name")

counts_table <- bracken_bac_counts %>%
	reduce(left_join, by = "name")

colnames(abund_table)[2:ncol(abund_table)] <- gsub(".bracken.out", "", basename(snakemake@input[["bracken"]]))
colnames(counts_table)[2:ncol(counts_table)] <- gsub(".bracken.out", "", basename(snakemake@input[["bracken"]]))

write.table(abund_table, file=snakemake@output[[1]], sep='\t', row.names = F)
write.table(counts_table, file=snakemake@output[[2]], sep='\t', row.names = F)
