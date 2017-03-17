import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sns.set_context("paper", font_scale=2)

#FILE = 'corr_res.csv'
FILES = ['HERA_243_p.csv', 'HERA_128_p.csv']
LABELS = ['HERA243', 'HERA128']

dflist = []
for i, file in enumerate(FILES):
	df = pd.read_csv(file)
	df['label'] = LABELS[i]
	df['peak'] /= np.amax(df['peak'])
	df['Theta'] = np.sqrt(df['mult'])*df['peak']
	df['Theta'] /= np.amax(df['Theta'])
	dflist.append(df)

df = pd.concat(dflist)

g = sns.pairplot(df,hue='label',vars=['dT','peak','Theta','bl1','bl2'],
	plot_kws={'alpha':0.2, "s":30})
# g = sns.PairGrid(df)
# g = g.map_diag(sns.kdeplot, lw=3)
# g = g.map_offdiag(sns.kdeplot, lw=1)
plt.show()


