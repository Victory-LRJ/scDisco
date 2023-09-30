from util.metrics_count_single import count

# ### #### read data
data_dir = '../../datasets/human_lung/'

methods = ["scDisco", "scDisco_ae"]
method_name = 'lung_'
count(data_dir, method_name, methods, single=False)

