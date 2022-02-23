import numpy as np
import pandas as pd
import sys
import scipy.stats
import sklearn.cluster
import os
import multiprocessing as mp


def _spearman(data_matrix, i, j, min_overlap, absolute):
    '''
    Internal function, called by calculate_spearman. Get the two rows of the data matrix
    that correspond to the two user inputs, get their overlap, and calculate Spearman's rho
    if the overlap is above our threshold.
    Return the tuple (i,j,R2). i and j are the numeric indices of the two user inputs.
    '''

    
    # Get all finite protein values that are present in both inputs
    x_y = data_matrix.iloc[[i,j],:]
    x_y_finite = x_y.dropna(axis = 1, how = 'any')
    x, y = x_y_finite.values
    
    if absolute:
        x = abs(x)
        y = abs(y)
    
    #import pdb; pdb.set_trace()
    
    overlap = len(x)
    
    if (overlap < min_overlap): 
        # Need to filter out small array lengths first because e.g. on an array of length 0 we cannot
        # calculate standard deviation.
        R2 = 0.0
    else: 
        # Check if any of the x and y are constant because this can happen for datasets with e.g. many
        # zeros, where the non-zeros are excluded because they are not present in the other input.
        # Constant input to the spearman calculation would result in NaN. Although technically incorrect, I
        # will replace this NaN by a zero because there is no evidence of the two user inputs being related.  
        
        std_x = np.std(x) # this can produce RuntimeWarnings when the mean is 0. Cannot be caught while running in parallel.
        std_y = np.std(y)
        
        # Set a more or less arbitrary tolerance threshold for the standard deviation under which we
        # consider the input to be constant.
        tolerance = 1e-12
        
        if (std_x < tolerance) or (std_y < tolerance):
            R2 = 0.0
        else:
            R, pval = scipy.stats.spearmanr(x,y) # any RuntimeWarnings will have been prevented
            R2 = R**2 # we want positive and negative correlations
    
    return (i, j, R2)


def _sym_diff(data_matrix, i, j):
    '''
    Internal function to calculate the symmetric difference between two user inputs.
    '''
    # Turn the two pandas rows into a boolean numpy array of the same dimension, 
    # stating whether they have a protein entry or not
    isnan = np.isnan(data_matrix.iloc[[i,j],:].values)
    
    # Each protein gets a True value only if it is in one row but not the other
    xor = np.logical_xor(isnan[0], isnan[1])
    
    # Sum the number of proteins that occur in one row but not the other
    symmetric_difference = np.sum(xor)
    
    return i, j, symmetric_difference





def calculate_metric(data_matrix,
                     metric_type = "spearman",
                     absolute = False,
                     min_overlap = 100, 
                     metric_matrix_dir = 'data/interim/metric_matrices/9606',
                     nr_cores = 10,
                     parallelization_threshold = 50):
    '''
    This function calculates correlation/ distance metrics between all user inputs based on the input values.
    Available are the (absolute) Spearman correlation R^2 and the Symmetric Difference.
    It does so in an incremental way, meaning that if a metric has already been calculated for some
    of the user inputs, these values can be provided to the function and don't need to be re-calculated.
    
    data_matrix:         A dataframe with dataIds as index, and protein shorthands as columns, filled with the
                         user-provided values (each row is one user input) or NaN for proteins which do not
                         appear in the particular input.
    
    metric_type:         One of "spearman" or "sym_diff". Whether to calculate Spearman's R^2 or Symmetric Difference.
    
    min_overlap:         Minimum number of proteins present in both dataIds in order to calculate a 
                         correlation coefficient. dataId pairs that have an overlap < min_overlap will get
                         an Rˆ2 of 0. Only for spearman.
    absolute:            Whether to calculate the correlation on the absolute input values or the original input values.
                         Only for spearman.
    metric_matrix_dir:   The directory where pre-calculated metrics dataId x dataId dataframes will be stored as tsv. 
    
    nr_cores:            The number of cores to be used in parallel computation.
    
    parallelization_threshold: The minimum number of user inputs for running the computation in parallel. 
                               Of course, parallel computation may not be required if there is only one new
                               dataId being added and all other values are already calculated. But at
                               parallelization_threshold = 50, this means that at least 49 correlations/distance
                               values will have to be calculated, which can for sure be parallelized.
    
    
    Returns: A dataId x dataId dataframe with the correlation / distance metric values for all dataId pairs. 
    '''
    
    if metric_type not in ["spearman", "sym_diff"]:
        sys.exit("metric_type has to be either 'spearman' or 'sym_diff'.")
        
    dataIds = data_matrix.index.to_list()
    dim = len(dataIds)
    metric_matrix = np.zeros((dim, dim)) + np.nan
    
    
    if (metric_type == "spearman") and not absolute:
        existing_metric_file = f"spearman.overlap_{min_overlap}.tsv"
    elif (metric_type == "spearman") and absolute:
        existing_metric_file = f"spearman_abs.overlap_{min_overlap}.tsv" 
    else:
        existing_metric_file = "sym_diff.tsv"
        
    existing_metric_file = os.path.join(metric_matrix_dir, existing_metric_file)
    
    
    if os.path.exists(existing_metric_file):
        existing_metric_df = pd.read_table(existing_metric_file, index_col=0, header=0)
        existing_dataIds = existing_metric_df.index.to_list()
    else:
        existing_metric_df = None
        existing_dataIds = []
  

    dataId_indices_to_calculate = [i for i, d in enumerate(dataIds) if d not in existing_dataIds] 
    
    
    # Enumerate all the combinations that have to be calculated and remove duplicates
    # Since we're not dealing with symmetric coordinates, we can't use the old "if i > j" trick for
    # calculating only one triangle.
    coords = [(i,j) for i in range(dim) for j in dataId_indices_to_calculate]
    coords = [tuple(sorted(pair)) for pair in coords]
    coords = list(set(coords))
     
    ## Check if it's worth it to parallelize. I could also make this dependent on how many coordinates
    ## have to be calculated.
    if len(data_matrix)  > parallelization_threshold :
        ## parallel computation
        
        print("Opening pool")
        pool = mp.Pool(nr_cores)
        
        try:
            
            if metric_type == "spearman":
                results = pool.starmap_async(_spearman, 
                                             [(data_matrix, i, j, min_overlap, absolute) for i,j in coords]).get()
            else:
                results = pool.starmap_async(_sym_diff, 
                                             [(data_matrix, i, j) for i,j in coords]).get()
                
            
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            pool.join()
            sys.exit()
            
        else:
            pool.close()
            print("Pool closed")
            pool.join()

        for i, j, metric in results:

            metric_matrix[i,j] = metric
            metric_matrix[j,i] = metric
    
    else:
        ## sequential_computation
        for i, j in coords:
            
            if metric_type == "spearman":
                i,j,metric = _spearman(data_matrix, i, j, min_overlap, absolute)
            elif metric_type == "sym_diff":
                i,j,metric = _sym_diff(data_matrix, i, j)
            else:
                sys.exit("metric_type has to be either 'spearman' or 'sym_diff'.")

            metric_matrix[i,j] = metric
            metric_matrix[j,i] = metric
 
    ## Convert matrix to dataframe
    metric_df = pd.DataFrame(metric_matrix, index = dataIds, columns = dataIds)
    

    ## Put the previously calculated data in the right place
    if existing_metric_df is not None:
        # Reducing the pre-calculated dataIds down to the ones that are actually in the 
        # current data to deduplicate.
        desired_existing_dataIds = list(set(existing_dataIds).intersection(set(dataIds)))
        
        idx = (desired_existing_dataIds, desired_existing_dataIds)

        assert(all(metric_df.loc[idx].isna()))
        
        metric_df.loc[idx] = existing_metric_df.loc[idx].values
        
#         assert(all(metric_df.loc[existing_dataIds, existing_dataIds].isna()))
        
#         metric_df.loc[existing_dataIds, existing_dataIds] = \
#                             existing_metric_df.loc[existing_dataIds, existing_dataIds].values
    
    if metric_type == "sym_diff":
        metric_df = metric_df.astype(int)
    
    # Write metrics to file so that the next run can find these data.
    os.makedirs(metric_matrix_dir, exist_ok=True)
    metric_df.to_csv(existing_metric_file, sep = '\t')
    
    return metric_df





def remove_similar_inputs_of_one_species(data_matrix, 
                                         min_overlap = 100, 
                                         r2_threshold = 0.8, 
                                         symm_diff_threshold = 2, 
                                         dedup_id_file = None, 
                                         metric_matrix_dir = 'data/interim/metric_matrices/9606',
                                         nr_cores = 10,
                                         parallelization_threshold = 50):
    """
    dedup_id_file: a file containing the taxId where all dataIds will be written to. 

    """
    nr_all_inputs = len(data_matrix)
    nr_of_inputs = {'full': nr_all_inputs}
    print(f"Number of user inputs: {nr_all_inputs}:")
    
    r2_dist_threshold = 1 - r2_threshold
    
    
    
    ### Calculcate distance metrics
    
    print("Calculating spearman")
    spearman_df = calculate_metric(data_matrix,
                                     metric_type = "spearman",
                                     min_overlap = min_overlap, 
                                     absolute = False,
                                     metric_matrix_dir = metric_matrix_dir,
                                     nr_cores = nr_cores,
                                     parallelization_threshold = parallelization_threshold)
    
    print("Calculating absolute spearman")
    abs_spearman_df = calculate_metric(data_matrix,
                                         metric_type = "spearman",
                                         min_overlap = min_overlap, 
                                         absolute = True,
                                         metric_matrix_dir = metric_matrix_dir,
                                         nr_cores = nr_cores,
                                         parallelization_threshold = parallelization_threshold)
    
    print("Calculating symmetric difference")
    sym_diff_df = calculate_metric(data_matrix,
                                   metric_type = "sym_diff",
                                   min_overlap = min_overlap, 
                                   metric_matrix_dir = metric_matrix_dir,
                                   nr_cores = nr_cores,
                                   parallelization_threshold = parallelization_threshold)
        
        
    
    metric_dict = dict(spearman_r2_dist     = dict(matrix    = 1 - spearman_df,
                                                   threshold = r2_dist_threshold),
                       
                       abs_spearman_r2_dist = dict(matrix    = 1 - abs_spearman_df,
                                                   threshold = r2_dist_threshold),
                       
                       symm_diff            = dict(matrix    = sym_diff_df,
                                                   threshold = symm_diff_threshold))
    
    
    
    if nr_all_inputs <= 1:
        non_duplicate_dataIds = sorted(data_matrix.index.to_list())
        # Create a fake cluster label dataframe with zero as cluster label for all similarity metrix 
        zero_array = np.zeros((len(data_matrix),len(metric_dict)))
        clustered_data = pd.DataFrame(zero_array, columns = metric_dict.keys(), index = data_matrix.index)
        
    else:

        ### Agglomerative clustering
        label_dict = dict()
        for metric in metric_dict:
            
            # Initialize a single linkage clustering
            sl_clustering = sklearn.cluster.AgglomerativeClustering(affinity='precomputed', 
                                                                    compute_full_tree='auto', 
                                                                    linkage='single', 
                                                                    distance_threshold=metric_dict[metric]['threshold'], 
                                                                    n_clusters = None)
            

            # Cluster
            sl_clustering.fit(X= metric_dict[metric]['matrix'])
            
            #import pdb; pdb.set_trace()
            
            # Get the cluster labels
            label_dict[metric] = sl_clustering.labels_


        ### Calculate number of proteins in input
        input_sizes = (~data_matrix.isna()).sum(axis = 1)
        assert( np.all(input_sizes.index.values == data_matrix.index.values))

        ### Combine input sizes and cluster labels into one dataframe
        clustered_data = pd.DataFrame(label_dict, index = data_matrix.index)
        clustered_data = clustered_data.assign(input_sizes = input_sizes)

        ### Sort by size and drop duplicate clusters
        ### We want the largest input to be kept. The other columns are just sorted by in order to get a deterministic outcome.
        ### The cluster labels per se do not provide any information so the sort order does not matter, only whether they are the same or not.
        deduplicated_data = clustered_data.sort_values(['input_sizes'] + list(metric_dict.keys()), ascending = False)

        ### Save the number of inputs after each duplication step
        for metric in metric_dict:
            deduplicated_data.drop_duplicates(subset = metric, inplace=True)
            print(f"After {metric}, {len(deduplicated_data) = }")
            nr_of_inputs[metric] = len(deduplicated_data)

        ### Return the non-duplicated IDs
        non_duplicate_dataIds = deduplicated_data.index.to_list()
        non_duplicate_dataIds = sorted(non_duplicate_dataIds)


    ## Write deduplicated IDs to file. 
    if dedup_id_file:
        with open(dedup_id_file, 'w') as out:
            id_string = '\n'.join(non_duplicate_dataIds)+'\n'
            out.write(id_string)
            
            
    clustered_inputs_file = os.path.join(metric_matrix_dir, "single_linkage_clusters.tsv")
    clustered_data.to_csv(clustered_inputs_file, sep = '\t')

        
    return clustered_data, nr_of_inputs, non_duplicate_dataIds


    
    




_, species_matrix_file, metric_matrix_dir, nr_cores, min_overlap, dedup_id_file  = sys.argv 
nr_cores = int(nr_cores)
min_overlap = int(min_overlap)

species_matrix = pd.read_table(species_matrix_file, index_col=0, header=0)

## For each species, calculate the metrics between inputs, cluster similar inputs, and keep only one cluster representative.  
clustered_data, nr_of_inputs, non_duplicate_dataIds = \
                        remove_similar_inputs_of_one_species(data_matrix = species_matrix, 
                                                                 min_overlap = min_overlap, 
                                                                 r2_threshold = 0.8, 
                                                                 symm_diff_threshold = 2, 
                                                                 dedup_id_file = dedup_id_file, 
                                                                 metric_matrix_dir = metric_matrix_dir,
                                                                 nr_cores = nr_cores,
                                                                 parallelization_threshold = 40)



#######################
##### TESTS ###########
#######################


## Read input matrix (columns are proteins, rows are inputs)
# species_matrix_file = "data/interim/species_matrices/9606.tsv"
# species_matrix = pd.read_table(species_matrix_file, index_col=0, header=0)

# ## Verify that there are no infinite values
# all(species_matrix.apply(lambda row: all(np.isfinite(row.dropna())), axis = 1))

# ## Verify that all inputs fit into float64 dtype
# all(species_matrix.apply(lambda row: row.dtype, axis = 1) == "float64")

# # Check min and max values
# all(species_matrix.apply(lambda row: row.max(), axis = 1) < 1e200)
# all(species_matrix.apply(lambda row: row.min(), axis = 1) > -1e200)

#ids = "DFwDoizXLnN3  FXLCF05mc7Gl  L00P8PBTTqNk  fZwwiHVs7Glc".split()


# clustered_data, nr_of_inputs, non_duplicate_dataIds = \
#                         remove_similar_inputs_of_one_species(data_matrix = species_matrix.iloc[0:700,:], 
#                                                                  min_overlap = 100, 
#                                                                  r2_threshold = 0.8, 
#                                                                  symm_diff_threshold = 2, 
#                                                                  dedup_id_file = None, 
#                                                                  metric_matrix_dir = 'test',#'data/interim/metric_matrices/9606',
#                                                                  nr_cores = 10,
#                                                                  parallelization_threshold = 10)



# spearman_df = calculate_metric(data_matrix = species_matrix.iloc[0:500,:],
#                                  metric_type = "spearman",
#                                  min_overlap = 100, 
#                                  absolute = False,
#                                  metric_matrix_dir = 'test',
#                                  nr_cores = 10,
#                                  parallelization_threshold = 50)

# sym_diff_df = calculate_metric(data_matrix = species_matrix.iloc[0:200,:],
#                                  metric_type = "sym_diff",
#                                  min_overlap = 100, 
#                                  absolute = False,
#                                  metric_matrix_dir = 'test',#'data/interim/metric_matrices/9606',
#                                  nr_cores = 2,
#                                  parallelization_threshold = 2)
