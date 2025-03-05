import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# load data
data = pd.read_csv('Individuals.csv')

# histograms and density plots to check distribution shape
for col in data.columns[1 : -1]:
    plt.figure(figsize=(6,4))
    sns.histplot(data[col].dropna(), bins=30, kde=True)
    plt.title(f'{col}')
    plt.xlabel('Value')
    plt.ylabel('Freq.')
    plt.show()
    