from util.metrics_count_single import count
#%% 1. read data
data_dir = '../../datasets/'
data_dir_names = [data_dir + 'pancreas2/', data_dir + 'lung/', data_dir + 'pancreas/',
                  data_dir + 'mouse_mucosa/', data_dir + 'epithelium/', data_dir + 'ductal/']
methods = ["scDisco", "scDisco_ae"]
method_names = ['pancreas2_', 'lung_', 'pancreas_', 'mouse_mucosa_', 'epithelium_', 'ductal_']

#%% 2. count metrics for scDisco-AE
for i in range(6):
    count(data_dir_names[i], method_names[i], methods, single=False)



