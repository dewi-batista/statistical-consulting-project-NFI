import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

# load data
csv_string = 'Mixtures'
data = pd.read_csv(f'{csv_string}.csv')

# replace fluids column with its one-hot encoding
fluids_col = data.columns[0]
one_hotted_fluids = pd.get_dummies(data[fluids_col], dtype=int)
data = data.drop(columns=[fluids_col])
data = pd.concat([one_hotted_fluids, data], axis=1)

# replace marker vals above threshold of 150 with 1 and 0 otherwise
num_fluids = one_hotted_fluids.shape[1]
for col_index in range(num_fluids, data.shape[1] - 1):
    data.iloc[:, col_index] = (data.iloc[:, col_index] >= 150).astype(int)

# save data ('%.0f' ensures that everything is int)
data.to_csv(f'preprop_{csv_string}.csv', index=False, float_format='%.0f')
