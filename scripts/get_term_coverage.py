# Count the number of terms each gene/protein is annotated with, per annotation source.
# Only terms smaller than 250 genes are counted in order to count specific terms only.
# An alternative is counting the terms and dividing each term count by the number of genes that
# are annotated by this term, but since this is harder to understand and interpret, we are 
# using a simple cutoff instead.
#
# In a second step, we are also counting the number of protein-coding genes covered by terms
# smaller than term size 250 (covered zero times, once, twice, or three+ times).


import os
import sys
import numpy as np
import pandas as pd
import utils

_, proteins_to_shorthands_file, protein_info_file, terms_members_file, output_dir, taxId, term_size_threshold = sys.argv

# taxId=9606
# lo, hi = 0, 250
# proteins_to_shorthands_file = "data/raw/proteins_to_shorthands.v11.tsv"
# protein_info_file = f"data/raw/{taxId}.protein.info.v11.0.txt.gz"
# terms_members_file = f"data/raw/global_enrichment_annotations/{taxId}.terms_members.tsv"
# output_dir = f"data/results/database_stats"

lo = 0
hi = int(term_size_threshold)
term_size_limits = (lo, hi)


# output files
term_coverage_file = os.path.join(output_dir, f"term_coverage_per_gene_{lo}-{hi}.tsv") 
genome_coverage_file = os.path.join(output_dir, f"genome_coverage_{lo}-{hi}.tsv")
percent_genome_coverage_file = os.path.join(output_dir, f"percent_genome_coverage_{lo}-{hi}.tsv")

etype_to_database = utils.etype_to_database
dbColors = utils.dbColors


def add_database_column(df, col, etype_to_database):
    
    col = str(col)
    etypes = df.loc[:,col].values
    
    # make sure the etype is given as a negative integer
    negative_etypes = [-abs(e) for e in etypes]
    
    database = [etype_to_database[e] for e in negative_etypes]
    df = df.assign(database = database)

    return df


def count_term_coverage_per_protein_by_database(protein_info_shorthands, terms_members_file,
                                                etype_to_database, term_size_limits):
    '''
    Count the occurrence of each protein in the term members file.

    protein_info_shorthands: A dataframe containing the complete set of proteins for 
                            the organism, with at least a 'shorthand' column and
                            optional other columns.
    terms_members_file: tsv file with one line per annotation term for the organism 
                        in question. Fields: [term_id, etype, term_size, 
                        term_member_shorthands]
    etype_to_database: dictionary with the translation of each entity type code aka etype to
                        the name of the annotation system aka database
    term_size_limits: Tuple. The minimum and maximum term size to be included in the coverage 
                      count.
    
    Returns dataframe with shorthands as index and 'normalized_term_count', 'term_count', 'database',
    and 'etype' as columns, joined with the protein info table.
    '''
    

    
    # get the term counts per protein in the annotation file
    term_coverage_dict = dict()
    
    with open(terms_members_file) as f:
        lines = f.readlines()
    
    for string in lines:
        fields = string.split()
        term_id   = int(fields[0])
        etype     = int(fields[1])
        term_size = int(fields[2])
        
        if (term_size > term_size_limits[1]) or (term_size < term_size_limits[0]):
            continue
            
        proteins  = fields[3:]
        
        if etype not in term_coverage_dict:
            term_coverage_dict[etype] = dict()
            
        for p in proteins:
            p = int(p)
            if p not in term_coverage_dict[etype]:
                term_coverage_dict[etype][p] = {'normalized_term_count': 0,
                                                'term_count': 0}
            term_coverage_dict[etype][p]['normalized_term_count'] +=  1/term_size
            term_coverage_dict[etype][p]['term_count'] += 1
                
    # turn each etype's subdict into a dataframe
    etype_frames = []
    for e in term_coverage_dict:
        df = pd.DataFrame.from_dict(term_coverage_dict[e], orient='index')
        df = df.merge(protein_info_shorthands, left_index=True, right_on='shorthand', how = 'right')
        df = df.assign(etype = e)
        etype_frames.append(df)
        
    # combine dataframes of all etypes
    term_coverage_df = pd.concat(etype_frames)
    term_coverage_df.loc[:,"etype"] = pd.Categorical(term_coverage_df.etype)
    
    
    # replace nan counts with zero
    term_coverage_df.loc[term_coverage_df.term_count.isna(), ['normalized_term_count', 'term_count']] = 0
    term_coverage_df = term_coverage_df.astype(dtype = {"term_count":"int64"}, copy = False)
    
    # add translation of etype into database
    term_coverage_df = add_database_column(df=term_coverage_df, 
                                           col="etype", 
                                           etype_to_database=etype_to_database)
    # make database column categorical
    term_coverage_df.loc[:,'database'] = pd.Categorical(term_coverage_df.database,
                                                       categories = dbColors.index,
                                                       ordered = True)
    return term_coverage_df



def get_stratified_genome_coverage(term_coverage_df):
    """
    Turn the dataframe with term count (i.e. how many terms of this database can this
    gene be found in) information into a wide-form dataframe with the number of genes
    that are covered by 0, 1, 2, or 3+ terms in a given database.
    
    term_coverage_df: dataframe with at least columns: [shorthand, term_count, database]
    """
    
    def classify_gene_coverage(n):
        group = str(int(n)) if n < 3 else '3+'
        return group

    gene_coverage_group = term_coverage_df.term_count.apply(classify_gene_coverage)
    term_coverage_df = term_coverage_df.assign(gene_coverage_group = gene_coverage_group)


    genome_coverage_df = term_coverage_df.groupby(['database', 'gene_coverage_group']).apply(len)
    genome_coverage_df = genome_coverage_df.unstack()

    # Sort columns
    sorted_columns = genome_coverage_df.columns.to_list()
    reverse_columns = genome_coverage_df.columns.sort_values(ascending = False).to_list()

    genome_coverage_df = genome_coverage_df.loc[:, reverse_columns]
    
    return genome_coverage_df








# Load protein info
protein_info = pd.read_table(protein_info_file, 
                             header = 0,
                             names = ['stringId', 
                                     'preferred_name', 
                                     'protein_size', 
                                     'annotation'])
proteins_to_shorthands = pd.read_table(proteins_to_shorthands_file, 
                                       header = None,
                                       names = ["stringId", "shorthand"])
protein_info_shorthands = protein_info.merge(proteins_to_shorthands)


# calculate how many terms each gene is annotated in in each database
term_coverage_per_gene = count_term_coverage_per_protein_by_database(protein_info_shorthands,
                                                                      terms_members_file,
                                                                      etype_to_database,
                                                                      term_size_limits )


# Check that for each gene/protein, there is one row per database
assert len(term_coverage_per_gene)/ len(dbColors) == len(protein_info_shorthands)


# Save to file
term_coverage_per_gene.to_csv(term_coverage_file, sep = "\t", index = False)



### Percentage of protein-coding genome covered by annotation
### Create coverage categories: covered 0 times, once, twice, or three+ times.

genome_coverage_df = get_stratified_genome_coverage(term_coverage_df = term_coverage_per_gene)
genome_coverage_df.to_csv(genome_coverage_file, sep = "\t", index = True)

percent_genome_coverage_df = 100*genome_coverage_df/len(protein_info_shorthands)
percent_genome_coverage_df.to_csv(percent_genome_coverage_file, sep = "\t", index = True)


