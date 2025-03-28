import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# load data
data = pd.read_csv('data/preproc_mixtures.csv')

# all six fluid combinations present in mixtures.csv
fluid_combinations = [
    # ('Semen.fertile', 'Vaginal.mucosa'),
    # ('Saliva', 'Vaginal.mucosa'),
    ('Blood', 'Nasal.mucosa'),
    ('Nasal.mucosa', 'Saliva'),
    # ('Blood', 'Vaginal.mucosa'),
    ('Blood', 'Menstrual.secretion')
]

# set starting index of markers in csv
marker_start_index = list(data.columns).index("HBB")
markers = data.columns[marker_start_index:]
for col in markers:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # Adjusted to 1 by 3
    # fig.suptitle(f'Distributions of {col} conditioned on fluid combinations', fontsize=16, fontweight='bold')
    
    for i, (fluid_1, fluid_2) in enumerate(fluid_combinations):
        ax = axes[i]  # Since it's a 1D array now
        restricted_data = data[(data[fluid_1] == 1) & (data[fluid_2] == 1)]
        sns.histplot(restricted_data[col].dropna(), bins=30, kde=False, ax=ax)
        ax.set_title(f'{fluid_1} & {fluid_2}', fontweight='bold')
        ax.set_xlabel('Value')
        ax.set_ylabel('Freq.')
    plt.tight_layout()
    # plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"figures/2D_hists/distributions_of_{col}.pdf")

