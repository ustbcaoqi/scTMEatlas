from cytotrace2_py.cytotrace2_py import *

input_path = './cytotrace2_Neu_input.txt'
annotation_path = './cytotrace2_Neu_clusters_annotation.txt'
species = 'human'

results =  cytotrace2(input_path,
                     annotation_path=annotation_path,
                     species=species,
                     batch_size = 10000,
                     smooth_batch_size = 1000,
                     disable_plotting = False,
                     disable_parallelization = False,
                     max_cores = None,
                     seed = 42
                     )