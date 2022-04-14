# This script will read the terms_members file into an R list and save it into one RDS file per category in a tmp folder.
# It is important that the terms_members file is sorted by entity type for this script to work.
# If the number of terms in the category exceeds chunk_size, an intermediate RDS is written every chunk_size lines.
# In the end, the tmp files are read back into R, combined into one list per category, and written out.
# I am doing the temporary saving in chunks because R gets so slow while building lists iteratively.
# With a few tweeks like appending the term names and contents separately and only combining them in the end, the
# chunking and tmp file writing could possibly be avoided, but it works, so that's good enough for me at the moment. 


###### Convert term members file to list of annotations.
read_term_list <- function(members_file, out_dir, taxid, nrows = -1L, chunk_size = 1e5L) {

        ## Build output file stub
        #rds_file <- sub(pattern = '.rds', replacement = '', x = rds_file, fixed = TRUE)

        tmp_dir <- file.path(out_dir, paste0('tmp_', taxid))
        if ( ! dir.exists(tmp_dir)){
            dir.create(tmp_dir)
        }
        rds_file <- file.path(tmp_dir, paste0(taxid, ".term_list"))


        ## Define an empty term list
        empty_list <- vector(mode = 'list')

        ## Define which etype corresponds to which term category 
        translation <- list(
                '-78' = 'STRINGclusters',
                '-57' = 'Reactome',
                '-56' = 'PubMed',
                '-55' = 'Pfam',
                '-54' = 'InterPro',
                '-53' = 'SMART',
                '-52' = 'KEGG',
                '-51' = 'UniProt',
                '-23' = 'GO_MF',
                '-22' = 'GO_CC',
                '-21' = 'GO_BP'
                )

        print(paste("Reading all lines of", members_file))
        line_list  <- readLines(con = members_file, n = nrows)
        term_list <- empty_list

        number_lines <- length(line_list)
        print(paste("Going to parse", number_lines, "lines."))
        
        categories <- c()
        i <- 0
        j <- 0
        is_first_line <- TRUE

        for (line in line_list) {
            fields <- strsplit(line, split = '\t')[[1]]

                term_id   <- fields[[1]]
                etype     <- fields[[2]]
                term_size <- fields[[3]]
                members   <- fields[4:length(fields)]

                # Skip any etypes (term categories) which we don't recognize
                if (! etype %in% names(translation)){
                    next
                }
                category <- translation[[etype]]

                ## We cannot check the previous line's category if we are on the first line,
                ## therefore we are setting it manually here.
                if (is_first_line) {
                    old_category <- category
                        is_first_line <- FALSE
                }

            ## Save the old category's term list to file
            if (category != old_category){
                partfile <- paste(rds_file, old_category, sprintf("%04d", j), 'rds', sep = '.')
                saveRDS(term_list, partfile)
                print(paste('Wrote', partfile))
                term_list <- empty_list
                j <- 0

                categories <- c(categories, old_category)    
                old_category <- category
            }

            ## Start the list for the new category
            term_list[[term_id]] <- members

                i <- i + 1
                if (i %% 1000 == 0) {
                    ## Report progress
                    percent <- round( i/number_lines, 2) * 100
                    print(paste(i, "of", number_lines, "terms parsed, or", percent, '%.'))
                }

            ## Write out a partial file if there are many terms
            ## or if it's the last line.
            if ((i %% chunk_size == 0) | ((i == number_lines) & (length(term_list) != 0)))   {
                partfile <- paste0(rds_file, '.', old_category, '.', sprintf("%04d", j), '.rds')
                saveRDS(term_list, partfile)
                print(paste('Wrote', partfile))
                term_list <- empty_list
                j <- j + 1
            }
        }
    ## Save the last category
    categories <- c(categories, category)
    print('All parsed.')

    print(paste("Reading the tmp files back in, combining them, and write them to", out_dir)) 
    for (cat in categories){
        print(cat)
        file_names <- list.files(tmp_dir, pattern=paste(taxid,"term_list", cat, "*", "rds", sep = "."), full.names=TRUE)
        rds_list <- lapply(file_names, readRDS)
        cat_term_list <- unlist(rds_list, recursive = FALSE, use.names = TRUE)
        cat_rds_file <- file.path(out_dir, paste(taxid, "term_list", cat, "rds", sep = "."))
        saveRDS(cat_term_list, cat_rds_file)
    }

    print("Removing tmp dir")
    unlink(tmp_dir, recursive = TRUE)
    print("Done.")
}




# ## input params
# members_file <- "data/raw/global_enrichment_annotations/9606.terms_members.tsv"
# nrows <- -1L # i.e. read all
# taxid <- 9606
# out_dir <- "data/interim/annotation_terms"
# chunk_size <- 1e4L


## write all stdout and stderr to a log file
log_file <- snakemake@log[["log_file"]]
log_fh <- file(log_file, open = "wt")
sink(log_fh ,type = "output")
sink(log_fh, type = "message")

## read other input parameters
members_file <- snakemake@input[["members_file"]]
flag_file <- snakemake@input[["flag_file"]]
nrows <- -1L # i.e. read all
taxid <- snakemake@wildcards[["taxId"]]
out_dir <- snakemake@params[["output_dir"]] #"data/interim/annotation_terms"
out_files <- snakemake@output
print(out_files)
chunk_size <- 1e4L

print(members_file)
print(out_dir)
print(taxid)

read_term_list(members_file = members_file, out_dir = out_dir, taxid = taxid, nrows = nrows, chunk_size = chunk_size)

## Since my function only produces output files for the database term categories which exist, I need to make empty .rds
## files for the rest of them so that I can assume a fixed number of term list rds files for each taxId.
for (file in out_files){
    if (! file.exists(file)){
        saveRDS('', file)
    }
}
#and to close connections
sink()
#sink()
