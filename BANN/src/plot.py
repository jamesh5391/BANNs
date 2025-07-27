import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot(x, y, title, out_name, clip=None):
	plt.figure(figsize=(6,6))
	if clip:
		plt.title(title + ' (Clipped)')
		plt.scatter(x.clip(upper=clip), y.clip(upper=clip), alpha=0.6, edgecolors='k')
	else:
		plt.title(title) 
		plt.scatter(x, y, alpha=0.6, edgecolors='k')
	plt.xlabel('-log10(p) PRSet')
	plt.ylabel('-log10(p) BANN') 
	plt.grid(True)
	plt.axline((0, 0), slope=1, color='gray', linestyle='--')
	plt.tight_layout()

	plt.savefig(out_name, dpi=300)

#RESULTS FOR MAIN EFFECTS 

bann_main_df = pd.read_csv('main_effect_relu_res_df.csv')
prset_main_df = pd.read_csv('main_effects_res_df.csv')
prset_main_df = prset_main_df.iloc[2:]

print(bann_main_df['pathway'].head)
print(prset_main_df['pathway'].head)

main_effects_bann = bann_main_df[bann_main_df['term'] == 'ppgs'] 
main_effects_prset = prset_main_df[prset_main_df['term'] == 'ppgs']

print(len(main_effects_bann))
print(len(main_effects_prset))

bann_logp = -np.log(main_effects_bann['p.value'].clip(lower=1e-200))
prset_logp = -np.log(main_effects_prset['p.value'].clip(lower=1e-200))

main_title = 'Comparison of Main Genetic Effect -log10(p) Values'
plot(prset_logp, bann_logp, main_title, "main_effects_scatter.png")
plot(prset_logp, bann_logp, main_title, "main_effects_scatter_clip30.png", clip=30)



corr = np.corrcoef(main_effects_bann['estimate'], main_effects_prset['estimate'])
np.savetxt('beta_corr.txt', corr)



#RESULTS FOR INTERACTION EFFECTS
bann_df = pd.read_csv('bann_ppgs_relu_res_df.csv') 
prset_df = pd.read_csv('prset_ppgs_res_df.csv')
prset_df = prset_df.iloc[3:]

bann_df = bann_df[bann_df['term'] == 'ppgs:bmi']
prset_df = prset_df[prset_df['term'] == 'ppgs:bmi']

print(bann_df['pathway'].head)
print(prset_df['pathway'].head)

bann_logp = -np.log10(bann_df['p.value'].clip(lower=1e-200))
prset_logp = -np.log10(prset_df['p.value'].clip(lower=1e-200))

interaction_title = 'Comparison of Interaction Effect -log10(p) Values'
plot(prset_logp, bann_logp, interaction_title, 'interaction_effects_scatter.png')
plot(prset_logp, bann_logp, interaction_title, 'interaction_effects_scatter_clip5.png', clip=5)


# RESULTS FOR SIGNIFICANCE RATIOS

def plotBox(x, title, out_name): 
	plt.figure(figsize=(6,6))
	plt.boxplot(x)
	plt.ylabel('Log (-log(p)) Ratios (BANN/PRset)')
	plt.title(title)
	plt.tight_layout()

	plt.savefig(out_name, dpi=300)

#ratios = bann_logp / prset_logp
#print(ratios.head(50))
#print(ratios.describe())
#log_ratios = np.log10(ratios.replace(0, np.nan).dropna())



z_ratios = bann_df['statistic'] / prset_df['statistic'].replace(0, pd.NA)
print(z_ratios.describe())
plotBox(z_ratios.dropna(), "Interaction Z-Ratios", "interaction_boxplot.png")
