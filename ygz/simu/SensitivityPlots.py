import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from numpy.random import random

sns.set_context("paper")
#sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)
plt.rc('axes', linewidth=2.5)

#FILE = 'corr_res.csv'
FILES = ['HERA_350_pm.csv', 'HERA_243_pm.csv', 'HERA_128_pm.csv', 'HERA_37_pm.csv','PAPER_128_pm.csv']
FILES = ['HERA_350_all.csv', 'HERA_243_all.csv', 'HERA_128_all.csv', 'HERA_37_all.csv','PAPER_128_all.csv']
FILES = ['HERA_350_pm.csv', 'HERA_128_pm.csv', 'PAPER_128_pm.csv']
LABELS = ['HERA350', 'HERA243', 'HERA128', 'HERA37', 'PAPER128']
LABELS = ['HERA350', 'HERA128', 'PAPER128']

#FILES = FILES[::-1]; LABELS = LABELS [::-1]
def gen_color(l=1):
	colors = []
	for i in range(l): colors.append((random(),random(),random()))
	return np.array(colors)
COLORS = gen_color(len(FILES))

def pairplot(Theta_min=0):

	sns.set(style='ticks', font_scale=1.5,font='DejaVu Serif')
	dflist = []
	for i, file in enumerate(FILES):
		df = pd.read_csv(file)
		df['Array'] = LABELS[i]
		df['$\Theta$'] = df['peak']
		df['$\Theta$'] /= np.amax(df['peak'])
		df['$L$'] = df['bl1']
		df['rho0'] = 0.001*40/df['bl1']
		df['$\widetilde{\Theta}$'] = np.sqrt(df['mult'])*df['peak']/np.sqrt(1+df['rho0']*2*np.sqrt(df['mult']))
		df['$\widetilde{\Theta}$'] /= np.amax(df['$\widetilde{\Theta}$'])
		df = df.loc[df['$\widetilde{\Theta}$']>Theta_min]
		df['$dT$'] = df['dT']
		dflist.append(df)

	df = pd.concat(dflist)
	plt.locator_params(axis='x', nticks=10)
	g = sns.pairplot(df,hue='Array',vars=['$dT$','$\Theta$','$\widetilde{\Theta}$','$L$'],
		plot_kws={'alpha':0.4, "s":30}, diag_kws={'histtype':"step", "linewidth":3})

	# g = sns.PairGrid(df,hue='Array',vars=['$dT$','$\Theta$','$\widetilde{\Theta}$','$L$'])
	# g = g.map_diag(plt.hist, histtype="step", linewidth=3)
	# #g = g.map_upper(sns.kdeplot)
	# g = g.map_upper(plt.scatter, alpha=0.5, s=10)
	# g = g.map_lower(plt.scatter, alpha=0.5, s=30)
	# g = g.add_legend()

	# ax0 = g.axes[0,0]
	# start, end = ax0.get_xlim()
	# start+=0.01; end+=0.01
	# ax0.set_xticks(np.linspace(start, end, 3))
	# ax1 = g.axes[0,1]
	# start, end = ax1.get_xlim()
	# ax1.set_xticks(np.linspace(start, end, 3))
	# ax2 = g.axes[0,2]
	# start, end = ax2.get_xlim()
	# ax2.set_xticks(np.linspace(start, end, 3))
	# ax3 = g.axes[0,3]
	# start, end = ax3.get_xlim()
	# ax3.set_xticks(np.linspace(start, end, 3))

	#import IPython; IPython.embed()
	plt.gcf().subplots_adjust(right=0.9)
	#g.legend.remove()
	# g = sns.PairGrid(df)
	# g = g.map_diag(sns.kdeplot, lw=3)
	# g = g.map_offdiag(sns.kdeplot, lw=1)

def get_imp(df, Theta_min=0.0):
	dft = df.loc[df['Theta']>Theta_min]
	dfeq = dft.loc[dft['sep']==dft['sep2']]
	dfnq = dft.loc[dft['sep']!=dft['sep2']]
	
	totalsumsq = np.sum(dft['Theta']**2)
	eqsumsq = np.sum(dfeq['Theta']**2)
	totalsum = np.sum(dft['Theta'])
	eqsum = np.sum(dfeq['Theta'])
	totalsens = totalsumsq#/totalsum*np.sqrt(len(dft.index))
	eqsens = eqsumsq#/eqsum*np.sqrt(len(dfeq.index))
	improve = (totalsens-eqsens)/eqsens
	return totalsens, eqsens

def sensplot():
	print "========= Statistics of sensitibity contribution =========="
	sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)

	plt.figure()
	plt.rc('axes', linewidth=4)
	for i, file in enumerate(FILES):
		df = pd.read_csv(file)
		df['rho0'] = 0.001*40/df['bl1']
		df['Theta'] = np.sqrt(df['mult'])*df['peak']/np.sqrt(1+df['rho0']*2*np.sqrt(df['mult']))
		df['Theta'] /= np.amax(df['Theta'])
		#df['Theta'] /= 1000

		TL = np.arange(0,1,0.01)
		 
		imp = np.array([get_imp(df, l) for l in TL])
		totalsensL, eqsensL = imp.T
		plt.plot(TL, totalsensL, label=LABELS[i], color=COLORS[i], linewidth=3)
		plt.plot(TL, eqsensL,'--', color=COLORS[i], linewidth=3)
	plt.legend()
	plt.xlabel(r'$\widetilde{\Theta}_{min}$')
	plt.ylabel(r'$\rho$')
	plt.gcf().subplots_adjust(bottom=0.2)
	plt.rc('axes', linewidth=4)


if __name__=="__main__":
	#sensplot()
	pairplot(0.05)

	plt.show()
	#import IPython; IPython.embed()

