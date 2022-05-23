## This script runs cameraPR, an enrichment function from the limma package, on one user input for one database.
## Potential error sources were an empty term list or empty pathway indices because cameraPR does not take care of this itself. 
## Therefore I needed to check those myself
library(tictoc)
## write all stdout and stderr to a log file
log_file <- snakemake@log[["log_file"]]
log_fh <- file(log_file, open = "wt")
sink(log_fh ,type = "output")
sink(log_fh, type = "message")

## concrete example, for debugging
## N1nQ0x555LUz.243159.KEGG
# infile <- "data/interim/filtered_deduplicated_user_inputs/N1nQ0x555LUz.243159.input.tsv"
## dedup_id_file <- "data/interim/deduplicated_dataIds/243159.tsv"
# database_file <- "data/interim/annotation_terms/243159.term_list.KEGG.rds"
# output_file   <- "data/results/cameraPR/enrichment/N1nQ0x555LUz.243159.KEGG.tsv"

#
# infile <- "data/interim/filtered_deduplicated_user_inputs/4sQkh9tb9a1E.9606.input.tsv"
## dedup_id_file <- "data/interim/deduplicated_dataIds/9606.tsv"
# database_file <- "data/interim/annotation_terms/9606.term_list.GO_BP.rds"

## input & output files
infile <- snakemake@input[["infile"]]
database_file <- snakemake@input[["database_file"]]
output_file   <- snakemake@output[["output_file"]]
# this file is not a necessary input, it is just here to produce a correctly connected dag
aggregated_infile <- snakemake@input[["aggregated_infile"]] 

## set minimum and maximum overlap parameters
# MIN_OVERLAP <- 3
# MAX_OVERLAP <- 200

MIN_OVERLAP <- as.numeric(snakemake@params[['min_overlap']])
MAX_OVERLAP <- as.numeric(snakemake@params[['max_overlap']])

tic("whole script")

empty_result <-  data.frame(list("term_shorthand" = vector(mode = "character"),
                                 "overlap"        = vector(mode = "numeric"),
                                 "direction"      = vector(mode = "character"),
                                 "p_value"        = vector(mode = "numeric"), 
                                 "q_value"        = vector(mode = "numeric"),
                                 "overlap_genes"  = vector(mode = "character")
                                 ))


filter_terms <- function(term_members, input_proteins, mi, ma){
    term_members_filtered <- term_members
    intersection <- intersect(term_members, input_proteins)
    overlap <- length(intersection)
    if(overlap < mi || overlap > ma){
      term_members_filtered <- NULL
    }
    return(term_members_filtered)
}


## read in RDS object which is a named list of shorthands, names are pathways
print("Reading gene sets.")
term_list <- readRDS(database_file)

if ((term_list == '') || (length(term_list) == 0)){
    print("No gene sets available, returning empty results.")
    ## create an empty result dataframe in case there are no terms for this database
    enrich_result <- empty_result
}else{
    print("Reading user input.")
    ## read in the user input file
    gene_value <- read.table(infile, col.names = c('shorthand', 'value'), stringsAsFactors = FALSE)
    gene_value_vector <- setNames(gene_value$value, gene_value$shorthand)
    

    ## Filter the term list to remove terms that have a too large or too small 
    ## overlap with the user input data, analogous to what STRING's global enrichment does.
    ## Only gene sets that have any overlap within the size constraints (MIN_OVERLAP 
    ## and MAX_OVERLAP) get included in enrichment analysis. This filtering does not change 
    ## the set size yet though because this is done by limma's enrichment function.

    if ((MIN_OVERLAP == 0) && (MAX_OVERLAP == Inf)){
       
        print("Not filtering gene sets because there are no boundaries on the overlap with user input.")
        # save time by skipping filtering
        term_list_filtered <- term_list

    }else{

        print("Filtering gene sets by input.")
        print(paste0("(Keep gene sets that overlap with input by minimum of ", 
                      as.character(MIN_OVERLAP), " and ", 
                      as.character(MAX_OVERLAP), ".)"))
        tic("filter term list")
        term_list_filtered <- lapply(term_list, filter_terms,
                                     input_proteins = names(gene_value_vector), 
                                     mi = MIN_OVERLAP, 
                                     ma = MAX_OVERLAP)
        # remove the NULL elements (these list elements should not be counted in multiple testing)
        term_list_filtered <- term_list_filtered[!sapply(term_list_filtered, is.null)]
        toc()
    }
  
    print("Running ids2indices to get genes actually present in the input.")
    tic("run ids2indices")
    ## get a named list of indices that correspond to the shorthands that are actually present in the given input
    pathway_indices <- limma::ids2indices(gene.sets = term_list_filtered, 
                                          identifiers = names(gene_value_vector),
                                          remove.empty = TRUE)
    toc()

    if (length(pathway_indices) == 0){
        print("No overlap between pathways and user input. Returning empty result.")
        # no overlap between the pathways and the user input
        enrich_result <- empty_result
    }else{
        print("Running cameraPR.")
        tic("run cameraPR")
        ## do camera preranked enrichment
        enrich_result_without_genes <- limma::cameraPR(statistic = gene_value_vector, 
                                         index = pathway_indices, 
                                         use.ranks = TRUE, 
                                         inter.gene.cor=0.01, 
                                         sort = TRUE)
        toc()

        ## FDR column will be missing if there is only one term.
        if(! "FDR" %in% colnames(enrich_result_without_genes)){
            FDR = enrich_result_without_genes$PValue[[1]]
            enrich_result_without_genes = cbind(enrich_result_without_genes, FDR)
        }
        
        tic("add overlap genes")
        print("Getting the enriched genes")
        # get the genes that are the overlap between the pathway and the input
        filter_by_existing_and_paste <- function(x, genes_in_input){
            filtered_genes <- genes_in_input[x]
            pasted_genes <- paste(filtered_genes, collapse = ' ')
            return(pasted_genes)
        }
        genes_in_input <- names(gene_value_vector)
        overlap_genes <- sapply(pathway_indices,
                                     filter_by_existing_and_paste, 
                                     genes_in_input = genes_in_input)

        print("Adding enriched genes to the enriched pathway result table.")
        # sort pathways by the order they are in in the enrich_result_without_genes table
        overlap_genes <- overlap_genes[rownames(enrich_result_without_genes)]
        
        # add the overlapping genes to the enrichment result dataframe
        overlapping_genes_per_term <- data.frame(overlap_genes)
        enrich_result <- merge(enrich_result_without_genes,
                                            overlapping_genes_per_term, 
                                            by = "row.names",
                                            sort = FALSE)
        toc()
        
        # rename the columns to something sensible
        colnames(enrich_result) <- c("term_shorthand",
                                     "overlap",
                                     "direction",
                                     "p_value",
                                     "q_value",
                                     "overlap_genes")

    }

}
    

## write table with enriched pathways and their pvals to file
print("Writing enrichment to file.")
tic("write table to file")
write.table(enrich_result, file = output_file, 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".",
            row.names = FALSE, col.names = TRUE)
toc()

toc()

## close connection to log file
sink()






# # I could use this in case I want to merge the enrichment stats with the overlap genes
# paste_genes <- function(x){return(paste(x, collapse = ' '))}
# overlap_genes <- sapply(term_list_filtered, paste_genes)
# overlapping_genes_per_term <- data.frame(overlap_genes)
# enrich_result_with_overlap <- merge(enrich_result, overlapping_genes_per_term, by = "row.names")


# # We can get the actual filtered shorthands from the pathway indices like so
# shorthands <- names(gene_value_vector)[pathway_indices[[1]]]


# write(shorthands, file = "test",
#       ncolumns = length(shorthands),
#       append = TRUE, 
#       sep = "\t")


# genes_in_list <- names(gene_value_vector)
# lapply(pathway_indices[1:3], function(x) genes_in_list[x])
       
       
# filter_by_existing_and_paste <- function(x){
#     filtered_genes <- genes_in_list[x]
#     pasted_genes <- paste(filtered_genes, collapse = ' ')
#     return(pasted_genes)
# }
       
# flattened_filtered <- sapply(pathway_indices[1:3],filter_by_existing_and_paste)
# pathways_and_contents <- paste(names(flattened_filtered), flattened_filtered, sep = '\t')

# writeLines(pathways_and_contents, con = "test")
