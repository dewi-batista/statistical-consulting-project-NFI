import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy  as np
import pandas as pd

def compute_dict_theta(data):
    
    # for assigning the fluid/marker columns
    marker_start_index = list(data.columns).index("HBB")

    fluid_cols = data.columns[: marker_start_index]
    marker_cols = data.columns[marker_start_index : -1] # disclude repl_val

    # extract how often each marker is present for each fluid, store in dict
    dict_theta = {}
    for fluid_name in fluid_cols:

        # strip dataframe down to just this fluid's rows
        fluid_rows = data[data[fluid_name] == 1]
        for marker_name in marker_cols:
            dict_theta[(fluid_name, marker_name)] = fluid_rows[marker_name].mean()
    
    return dict_theta

# essentially computes the probability of the fluids not appearing (prod of probs)
# then returns the probability of at least observing the marker 
def compute_bernoulli_param(dict_theta, fluids, marker):
    
    prod = 1
    for fluid in fluids:
        prod *= (1 - dict_theta[(fluid, marker)])
    return 1 - prod

if __name__ == "__main__":

    # load data
    csv_string = 'Individuals'
    data = pd.read_csv(f'preprop_{csv_string}.csv')

    # compute dictionary of theta parameters
    dict_theta = compute_dict_theta(data)

    # example usage
    fluids = ['Semen.fertile', 'Vaginal.mucosa']
    marker = 'ALAS2'
    marker_prob = compute_bernoulli_param(dict_theta, fluids, marker)
    print("Probability of at least one ", np.round(marker_prob, 3))

    # print(dict_theta[('Blood', 'ALAS2')], dict_theta[('Saliva', 'ALAS2')])
    # print(1 - (1 - 0.9603) * (1 - 0.0093))