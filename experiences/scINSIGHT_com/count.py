from metrics_com import count

# %% read data
data_dir = '../../datasets/'
data_dir_ = '../../datasets/scINSIGHT_to_scDisco/'
datasets = ["simulate", "pancreas2", 'lung', 'pancreas', 'mouse_mucosa']
method_names = ['scINSIGHT', 'scDisco']

# %% count metrics
count(data_dir, data_dir_, datasets, method_names)

