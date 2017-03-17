import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sns.set_context("paper", font_scale=2)

#FILE = 'corr_res.csv'
FILE = 'HERA_350_core_pm.csv'

df = pd.DataFrame.from_csv(FILE)
df['peak'] /= np.amax(df['peak'])
df['weighted_peak'] = df['mult']*df['peak']
df['weighted_peak'] /= np.amax(df['weighted_peak'])


g = sns.pairplot(df)
# g = sns.PairGrid(df)
# g = g.map_diag(sns.kdeplot, lw=3)
# g = g.map_offdiag(sns.kdeplot, lw=1)
plt.show()