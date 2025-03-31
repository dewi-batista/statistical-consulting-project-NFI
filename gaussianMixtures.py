# Good for all
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# For stats things
from scipy.stats import multivariate_normal, ks_2samp
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import train_test_split

# Python quality of life
from pdb import set_trace
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

def gauss_mixture_fit(data_set, idx_marker, idx_sample, 
                       n_components=1, covariance_type="diag",
                       split_proportion = 0.999,
                       plot_=True, do_model_selection=False, return_BIC = False, 
                       verbose = 2):
    """
    This function fits a mixture of Gaussians. If enabled, it will plot the histogram overlapped
    with the fitted distribution, as well as perform model selection for the optimal number of components.

    The documentation for GaussianMixture: 
    https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture

    --- Inputs:

    data_set = string. If it equals "individuals", then it uses the data set from individual fluids.
    Otherwise it uses the mixtures

    idx_marker = integer. The index of the relevant marker in the list of markers from the data
    The list: [HBB,ALAS2,CD93,HTN3,STATH,BPIFA1,MUC4,MYOZ1,CYP2B7P1,MMP10,MMP7,MMP11,SEMG1,KLK3,PRM1]

    idx_sample = integer. The index of the relevant fluid type in the list of fluids. The list:
    - individuals: ['Blood' 'Saliva' 'Vaginal.mucosa' 'Menstrual.secretion' 'Semen.fertile'
    'Semen.sterile' 'Nasal.mucosa' 'Blank_PCR' 'Skin' 'Skin.penile']
    - mixtures: ['Semen.fertile+Vaginal.mucosa' 'Saliva+Semen.fertile'
    'Saliva+Vaginal.mucosa' 'Blood+Nasal.mucosa' 'Nasal.mucosa+Saliva'
    'Vaginal.mucosa+Blood' 'Menstrual.secretion+Blood']

    split_proportion = float in the range (0, 1). It describes the proportion of the samples to be alocated
    to the training of the Gaussian mixture

    n_components = integer. The number of components to be used for the Gaussian mixture, set by the user.
    Relevant only if 'do_model_selection' is set to False

    covariance_type = string. Choose between ["diag", "full", "spherical", "tied"]. It gives the 
    type of matrix to be used in the mixture definitions. 

    plot_ = boolean. If true, it also plots a histogram of the relevant data together with the best mixture fit based on the
    training data 

    do_model_selection = boolean. If true, it performs model selection to find the best number of mixture components.
    **Note**: it uses BIC, which is the stricter criterion. It also has a maximum mixture number of 10. If those need to 
    be updated, this is an avenue of improvement

    return_BIC = boolean. If true, it returns the BIC for the best model computed on the training data

    verbose = 0, 1 or 2. If 2, the function prints everything. If 1, it prints only the convergence status. If 0, it prints nothing

    --- Output:

    model = model class of Gaussian Mixtures. It returns the model with the best number of components (user specified or via BIC)
    (Optional) best_BIC = float. The BIC for the best model, on the training data
    """
    
    # - Load data
    if data_set == "individuals":
        data = pd.read_csv('data/individuals.csv').fillna(0)
    else:
        data = pd.read_csv('data/mixtures.csv').fillna(0)

    markers = data.columns[1:-5]
    samples = data.iloc[:, 0].unique()

    data_sample = data[data.iloc[:, 0] == samples[idx_sample]]
    X = data_sample.loc[:, [markers[idx_marker]]]
    X_train, _ = train_test_split(X, train_size=split_proportion, random_state=42, shuffle=True)

    if do_model_selection == False:
        model = GaussianMixture(n_components=n_components, covariance_type=covariance_type, 
                                n_init=3, random_state=0).fit(X_train) 
        if verbose == 2:
            print(f"Number of components: {n_components}")
    else:
        max_components = 10
        BIC_list = np.zeros(max_components)

        for i in range(1, max_components + 1):
            model_temp = GaussianMixture(n_components=i, covariance_type=covariance_type, 
                                n_init=3, random_state=0).fit(X_train) 
            BIC_list[i-1] = model_temp.bic(X_train)

        # Select the model with the smallest BIC

        best_no_comp = np.argmin(BIC_list) + 1
        model = GaussianMixture(n_components=best_no_comp, covariance_type=covariance_type, 
                                n_init=3, random_state=0).fit(X_train)
        if verbose == 2:
            print(f"Number of components: {best_no_comp}") 

    if verbose == 1 or verbose == 2:
        print(f"Convergence status: {model.converged_}")

        if verbose == 2:
            print(f"Weights: {model.weights_}")
            print(f"Means: {model.means_}")
            print(f"Covariances: {model.covariances_}")

    # - Plot the fitted curve:
    if plot_ == True:
        # Generate the curve

        x_val = np.arange(data_sample[markers[idx_marker]].min(), data_sample[markers[idx_marker]].max(), 1) 
        y_val = np.zeros((len(model.weights_), x_val.shape[0]))
        y_all = np.zeros(x_val.shape)

        # Compute each individual contribution:
        for i in range(len(model.weights_)):
            weight_ = model.weights_[i]
            mean_ = model.means_[i]

            if covariance_type == "full":
                covar_ = model.covariances_[i] 
            elif covariance_type == "spherical":
                covar_ = model.covariances_[i]
            elif covariance_type == "tied":
                covar_ = model.covariances_
            else:
                covar_ = np.diag(model.covariances_[i])
                
            y_val[i, :] = weight_ * multivariate_normal.pdf(x_val, mean=mean_, cov=covar_)
            y_all = y_all + y_val[i, :]


        plt.figure(figsize=(6,4))

        hist_vals, _, _ = plt.hist(data_sample[markers[idx_marker]], bins=30, label="True Data", density=True)
        plt.plot(x_val, y_all, color="blue", label="Sum Contrib")

        for i in range(len(model.weights_)):
            plt.plot(x_val, y_val[i, :], label=f"Component {i+1}")

        plt.title(f'{samples[idx_sample]} - {markers[idx_marker]}\nTrue data and fit')
        plt.xlabel('Value')
        plt.ylabel('Freq.')
        plt.ylim(0, 1.3 * hist_vals.max())
        plt.legend()
        plt.show()
    
    if return_BIC == False:
        return model
    else:
        return model, BIC_list[best_no_comp - 1]

def data_generator(n, model, covariance_type="diag", threshold=150, seed=42):
    """
    This function generates data artificially, following the Gaussian mixture model from 'model'.
    The number of samples generated is 'n'
    The function also implements a thresolding mechanic: if a generated sample would be below 150, it is 
    returned as zero
    """

    np.random.seed(seed)

    # Initialize
    weights = model.weights_
    new_samples = np.zeros(n)
    temp1 = 0

    # Sample mixture components to take from 
    temp2 = np.random.choice(range(len(weights)), size=n, p=weights)
    mixture_number = np.bincount(temp2, minlength=len(weights)) 

    # Sample values from mixtures
    for i in range(len(weights)):
        # Extract curve parameters for relevant mixture component
        mean_ = model.means_[i]

        if covariance_type == "diag":
            covar_ = np.diag(model.covariances_[i])
        elif covariance_type == "full":
            covar_ = model.covariances_[i] 
        elif covariance_type == "spherical":
            covar_ = model.covariances_[i]
        else:
            covar_ = model.covariances_
        
        # Sample
        new_samples[temp1:(temp1+mixture_number[i])] = multivariate_normal.rvs(mean=mean_, cov=covar_, size=mixture_number[i])
        temp1 = temp1 + mixture_number[i]

    # Do thresholding
    new_samples[new_samples <= threshold] = 0

    return new_samples

def fit_evaluator(data_set, idx_marker, idx_sample, n_reps, split_proportion = 0.66,
                n_components=1, covariance_type="diag", do_model_selection=False,
                threshold=150, seed=42, plot_fit = False, plot_data = False): 
    """
    This function tests the adequacy of the generated data against the true data in few different ways.
    
    Firstly, it outputs the BIC of the best fit on the training data, which is useful to compare with other 
    distribution families (e.g. Gamma mixtures). 

    Secondly, it performs the two-sample KS test between the true test data and generated data sets of the same
    size as the test data (the KS test works better between samples of comparable size). This is done repeatedly 
    for multiple seeds, to obtain a distribution of p-values against different generation mechanisms.

    Lastly, it plots the histogram of the entire data overlapped with the best fit from the training data. 
    It helps with the visual comparison.

    --- Inputs:

    data_set = string. If it equals "individuals", then it uses the data set from individual fluids.
    Otherwise it uses the mixtures

    idx_marker = integer. The index of the relevant marker in the list of markers from the data
    The list: [HBB,ALAS2,CD93,HTN3,STATH,BPIFA1,MUC4,MYOZ1,CYP2B7P1,MMP10,MMP7,MMP11,SEMG1,KLK3,PRM1]

    idx_sample = integer. The index of the relevant fluid type in the list of fluids. The list:
    - individuals: ['Blood' 'Saliva' 'Vaginal.mucosa' 'Menstrual.secretion' 'Semen.fertile'
    'Semen.sterile' 'Nasal.mucosa' 'Blank_PCR' 'Skin' 'Skin.penile']
    - mixtures: ['Semen.fertile+Vaginal.mucosa' 'Saliva+Semen.fertile'
    'Saliva+Vaginal.mucosa' 'Blood+Nasal.mucosa' 'Nasal.mucosa+Saliva'
    'Vaginal.mucosa+Blood' 'Menstrual.secretion+Blood']

    n_reps = integer. The number iterations artificial samples are generated and compared to the true data

    split_proportion = float in the range (0, 1). It describes the proportion of the samples to be alocated
    to the training of the Gaussian mixture

    n_components = integer. The number of components to be used for the Gaussian mixture, set by the user.
    Relevant only if 'do_model_selection' is set to False.

    covariance_type = string. Choose between ["diag", "full", "spherical", "tied"]. It gives the 
    type of matrix to be used in the mixture definitions. 

    do_model_selection = boolean. If true, it performs model selection to find the best number of mixture components.
    **Note**: it uses BIC, which is the stricter criterion. It also has a maximum mixture number of 10. If those need to 
    be updated, this is an avenue of improvement

    plot_fit = boolean. If true, it also plots a histogram of the relevant data together with the best mixture fit
    based on only the full data

    plot_data = boolean. If true, it plots a histogram like the one in plot_fit, but with the generated data added as
    well. Furthermore, only training data is used for those mixtures. The data is generated based on the seed 'seed'

    threshold = float. The threshold set for the data generation. If a value is below the threshold, it is 
    assigned zero

    --- Outputs:
    """
    
    # --- Train the model on the training data
    model, best_BIC = gauss_mixture_fit(data_set, idx_marker, idx_sample, split_proportion=split_proportion, 
                              n_components=n_components, covariance_type=covariance_type, plot_=False, 
                              do_model_selection=do_model_selection, return_BIC=True, verbose=1)
    
    print(f"The BIC of the best model on training data: {best_BIC}")

    # --- Do the KS test for multiple seeds

    # Split the data just like gauss_mixture_fit()

    if data_set == "individuals":
        data = pd.read_csv('data/individuals.csv').fillna(0)
    else:
        data = pd.read_csv('data/mixtures.csv').fillna(0)

    markers = data.columns[1:-5]
    samples = data.iloc[:, 0].unique()
    data_sample = data[data.iloc[:, 0] == samples[idx_sample]]
    X = data_sample.loc[:, [markers[idx_marker]]]
    _, X_test = train_test_split(X, train_size=split_proportion, 
                                       random_state=42, shuffle=True)

    # Generate repeated artificial data sets. Record the p-values

    p_val_list = np.zeros(n_reps)
    seeds = np.random.randint(1, 3 * n_reps, n_reps)

    for i in range(n_reps):
        temp_artif_data = data_generator(X_test.size, model, threshold=threshold, seed=seeds[i])
        OUT = ks_2samp(X_test.to_numpy().flatten(), temp_artif_data)
        p_val_list[i] = OUT.pvalue

    print(f"The minimum p-value: {np.round(p_val_list.min(), 3)}")
    print(f"The median p-value: {np.round(np.median(p_val_list), 3)}")

    # --- Plot the data against the fit
    
    if plot_fit == True:
        _ = gauss_mixture_fit(data_set, idx_marker, idx_sample, split_proportion=0.999, 
                              n_components=n_components, covariance_type=covariance_type, plot_=True, 
                              do_model_selection=do_model_selection, return_BIC=False, verbose=1)

    if plot_data == True:
        artificial_data = data_generator(X_test.size, model, threshold=threshold, seed=seed)
        
        # Create the mixture curves
        x_val = np.arange(X.to_numpy().min(), X.to_numpy().max(), 1) 
        y_val = np.zeros((len(model.weights_), x_val.shape[0]))
        y_all = np.zeros(x_val.shape)

        # Compute each individual contribution:
        for i in range(len(model.weights_)):
            weight_ = model.weights_[i]
            mean_ = model.means_[i]
            covar_ = np.diag(model.covariances_[i])

            y_val[i, :] = weight_ * multivariate_normal.pdf(x_val, mean=mean_, cov=covar_)
            y_all = y_all + y_val[i, :]

        # The actual figure
        plt.figure(figsize=(6,4))

        hist_vals1, _, _ = plt.hist(artificial_data, bins=100, label="Generated Data", density=True, color="red", alpha=0.5)
        hist_vals2, _, _, = plt.hist(X, bins=50, label="True Data", density=True, color="blue", alpha=0.5)
        hist_max = max(hist_vals1.max(), hist_vals2.max())

        plt.plot(x_val, y_all, color="black", label="Sum Contrib")

        for i in range(len(model.weights_)):
            plt.plot(x_val, y_val[i, :], label=f"Component {i+1}")

        plt.title(f'{samples[idx_sample]} - {markers[idx_marker]}\nTrue data and generated data')
        plt.xlabel('Value')
        plt.ylabel('Freq.')
        plt.ylim(0, 1.3 * hist_max)
        plt.legend()
        plt.show()

if __name__ == "__main__":
    # --- Check each fluid-marker pair

    fit_evaluator("mixtures", 0, 6, 100, do_model_selection=True, seed=42, plot_data=True, plot_fit=True)
    
