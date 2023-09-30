from util.metrics_count_single import count

# ### #### read data
data_dir = '../../datasets/human_ductal/'

methods = ["scDisco", "scDisco_ae"]
method_name = 'ductal_'
count(data_dir, method_name, methods, single=False)

