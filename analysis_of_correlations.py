import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# load data
data = pd.read_csv('data/preproc_mixtures.csv')

# all six fluid combinations present in mixtures.csv
fluid_combinations = [
    ('Semen.fertile', 'Vaginal.mucosa'),
    # ('Saliva', 'Vaginal.mucosa'),
    # ('Blood', 'Nasal.mucosa'),
    # ('Nasal.mucosa', 'Saliva'),
    # ('Blood', 'Vaginal.mucosa'),
    # ('Blood', 'Menstrual.secretion')
]

# set starting index of markers in csv
marker_start_index = list(data.columns).index("HBB")
markers = data.columns[marker_start_index :]

# compute greatest correlates between markers for all fluid combinations
def greatest_correlates(corr_threshold, max_num_correlates):
    
    for (fluid_1, fluid_2) in fluid_combinations:
        
        # restrict dataset to just the rows where these two fluids appear
        restricted_data = data[(data[fluid_1] == 1) & (data[fluid_2] == 1)]

        # compute correlation matrix
        corr = restricted_data[markers].loc[:, (restricted_data[markers] != 0).any(axis=0)].corr()

        # only consider markers that actually exist in the correlation matrix
        non_zero_markers = corr.index.tolist()
        print(fluid_1, fluid_2)
        for marker in non_zero_markers:
            top_n_correlates = corr.loc[marker].abs().nlargest(max_num_correlates + 1)
            top_n_correlates = top_n_correlates[top_n_correlates > corr_threshold].index.tolist()
            if marker in top_n_correlates:
                top_n_correlates.remove(marker)
            print(marker, top_n_correlates)
        print()

        # plot that correlation matrix boiiii
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap="coolwarm", vmin=-1, vmax=1)
        # plt.title(f"Correlation matrix of markers for {fluid_1} and {fluid_2}")
        plt.tight_layout()
        plt.savefig('figures/corr_matrix_semen.fertile_and_vaginal.mucosa.pdf')

# greedy clustering algorithm
def get_clusters(corr_threshold):
    
    clusters_dict = {}
    for (fluid_1, fluid_2) in fluid_combinations:
        restricted_data = data[(data[fluid_1] == 1) & (data[fluid_2] == 1)]
        corr = restricted_data[markers].loc[:, (restricted_data[markers] != 0).any(axis=0)].corr()

        non_zero_markers = corr.index.tolist()
        pairs = []
        for i in range(len(non_zero_markers)):
            for j in range(i + 1, len(non_zero_markers)):
                val = corr.iloc[i, j]
                pairs.append((non_zero_markers[i], non_zero_markers[j], val))

        pairs.sort(key=lambda x: abs(x[2]), reverse=True)
        used = set()
        clusters = []
        for m1, m2, c in pairs:
            if m1 not in used and m2 not in used and abs(c) > corr_threshold:
                clusters.append((m1, m2, c))
                used.add(m1)
                used.add(m2)
        # print(fluid_1, fluid_2, "clusters:", clusters)
        clusters_dict[f"{fluid_1}+{fluid_2}"] = clusters

    return clusters_dict

if __name__ == "__main__":
    corr_threshold = 0.2
    print(get_clusters(corr_threshold))
    greatest_correlates(corr_threshold, max_num_correlates=3)
