from util.metrics_count import count

# %% read data
data_dir = '../../../datasets/simulate/sampling/'
methods = ["seurat", "harmony", 'scanorama', 'desc', 'scVI', 'cell_blast', 'scidrl', 'scDisinfact', "scDisco"]
method_name = 'simulate_'

# %% count metrics
count(data_dir, method_name, methods)

