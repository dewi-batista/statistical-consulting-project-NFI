import pandas as pd

# load data
csv_string = 'mixtures'
data = pd.read_csv(f'data/{csv_string}.csv')

# replace fluids column with its one-hot encoding
fluids_col = data.columns[0]

# how we preprocess differs slightly between individuals.csv and mixtures.csv
if csv_string == 'individuals':
    one_hotted_fluids = pd.get_dummies(data[fluids_col], dtype=int)
else:
    one_hotted_fluids = data[fluids_col].str.get_dummies(sep='+')

# replace the single fluids column with the one-hot encodings
data = data.drop(columns=[fluids_col])
data = pd.concat([one_hotted_fluids, data], axis=1)

# replace NaNs with 0s
data.fillna(0, inplace=True)

# FOR MIXTURES ONLY - sanitise! conditioned on each fluid combination
# get rid of these random values, e.g. conditioned on Semen.fertile and Vaginal.mucosa
# CD93 is 0 for 61/62 samples and is ~155 for 1/62, stupid anomaly!
# to combat this, get rid of those values, thresholded at 10%
if csv_string == 'mixtures':

    # all six fluid combinations present in mixtures.csv
    fluid_combinations = [
        ('Semen.fertile', 'Vaginal.mucosa'),
        ('Saliva', 'Vaginal.mucosa'),
        ('Blood', 'Nasal.mucosa'),
        ('Nasal.mucosa', 'Saliva'),
        ('Blood', 'Vaginal.mucosa'),
        ('Blood', 'Menstrual.secretion')
    ]
    
    # set starting index of markers in csv
    marker_start_index = list(data.columns).index("HBB")
    markers = data.columns[marker_start_index:]
    for (fluid_1, fluid_2) in fluid_combinations:
        restricted_data = data[(data[fluid_1] == 1) & (data[fluid_2] == 1)].copy()
        for col in markers:
            non_zero_ratio = (restricted_data[col] != 0).sum() / len(restricted_data)
            # print(non_zero_ratio)
            if non_zero_ratio <= 0.1:
                data.loc[(data[fluid_1] == 1) & (data[fluid_2] == 1), col] = 0

# drop all columns right of PRM1
markers_to_remove_start_index = list(data.columns).index("RPS4Y1")
markers_to_remove = data.columns[markers_to_remove_start_index:]
data = data.drop(columns=markers_to_remove)

# save data ('%.0f' ensures that everything is int)
data.to_csv(f'data/preproc_{csv_string}.csv', index=False)#, float_format='%.0f')
