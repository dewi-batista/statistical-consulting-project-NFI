import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# load data
data = pd.read_csv('data/preproc_mixtures.csv')

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
for col in markers:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    for i, (fluid_1, fluid_2) in enumerate(fluid_combinations):
        row, col_idx = divmod(i, 3)
        ax = axes[row][col_idx]
        restricted_data = data[(data[fluid_1] == 1) & (data[fluid_2] == 1)]
        sns.histplot(restricted_data[col].dropna(), bins=30, kde=False, ax=ax)
        ax.set_title(f'{fluid_1} & {fluid_2}', fontweight='bold')
        ax.set_xlabel('Value')
        ax.set_ylabel('Freq.')

    plt.tight_layout()
    plt.savefig(f"figures/2D_histograms/distributions_of_{col}.pdf")
