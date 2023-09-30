from util.metrics_count_single import count

# %% data dir
data_dir = '../../../datasets/ductal/'
methods = ["harmony", 'scanorama', 'desc', 'scVI', 'cell_blast', 'scidrl', 'scDisinfact', "scDisco"]
method_name = 'ductal_'

# %% count metrics
count(data_dir, method_name, methods, single=True, sample=False)

