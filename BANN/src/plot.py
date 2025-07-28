import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 


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
print("#######HI########")
print("CWD:", os.getcwd())
bann_main_df = pd.read_csv('BANNS/BANN/src/main_effect_relu_res_df.csv')
prset_main_df = pd.read_csv('BANNS/BANN/src/main_effects_res_df.csv')
prset_main_df = prset_main_df.iloc[2:]


main_effects_bann = bann_main_df[bann_main_df['term'] == 'ppgs'] 
main_effects_prset = prset_main_df[prset_main_df['term'] == 'ppgs']

#Filter by PRSet p-value < 0.05
main_effects_prset_significant = main_effects_prset[main_effects_prset['p.value'] < 0.05]

# keeps only the pathways that are significant in PRSet AND also exist in the BANN results.
merged_df = pd.merge(main_effects_bann, main_effects_prset_significant, on='pathway', suffixes=('_bann', '_prset'))

bann_logp = -np.log10(merged_df['p.value_bann'].clip(lower=1e-200))
prset_logp = -np.log10(merged_df['p.value_prset'].clip(lower=1e-200))

main_title = 'Comparison of Main Genetic Effect -log10(p) Values (PRSet p < 0.05)'
plot(prset_logp, bann_logp, main_title, "main_effects_scatter.png")
plot(prset_logp, bann_logp, main_title, "main_effects_scatter_clip30.png", clip=30)

# Calculate, sort, and print top main effect ratios
main_ratios = bann_logp / prset_logp.replace(0, np.nan)
main_ratios.index = merged_df['pathway'] # Assign pathway names for readability
top_main_ratios = main_ratios.sort_values(ascending=False)

print("--- Top 10 Pathways with Highest BANN/PRSet Main Effect Ratios ---")
print(top_main_ratios.head(10))
print("-" * 60)

bann_better_main = np.sum(merged_df['p.value_bann'] < merged_df['p.value_prset'])
total_main = len(merged_df)
percent_bann_better_main = (bann_better_main / total_main) * 100
print(f"--- Main Effect Significance ---")
print(f"BANN had greater statistical significance in {percent_bann_better_main:.2f}% of pathways.")
print(f"Compared a total of {total_main} pathways.")
print("-" * 30)

corr = np.corrcoef(main_effects_bann['estimate'], main_effects_prset['estimate'])
np.savetxt('beta_corr.txt', corr)



#RESULTS FOR INTERACTION EFFECTS
bann_df = pd.read_csv('BANNS/BANN/src/bann_ppgs_relu_res_df.csv') 
prset_df = pd.read_csv('BANNS/BANN/src/prset_ppgs_res_df.csv')
prset_df = prset_df.iloc[3:]

bann_df = bann_df[bann_df['term'] == 'ppgs:bmi']
prset_df = prset_df[prset_df['term'] == 'ppgs:bmi']


# CHANGE 1: Filter only the PRSet dataframe for p-values < 0.05
prset_df_significant = prset_df[prset_df['p.value'] < 0.05]

# CHANGE 2: Perform an inner merge to align the dataframes
# This keeps only the pathways that are significant in PRSet AND also exist in the BANN results.
merged_df = pd.merge(bann_df, prset_df_significant, on='pathway', suffixes=('_bann', '_prset'))

# CHANGE 3: Create logp values from the new merged dataframe
bann_logp = -np.log10(merged_df['p.value_bann'].clip(lower=1e-200))
prset_logp = -np.log10(merged_df['p.value_prset'].clip(lower=1e-200))


interaction_title = 'Comparison of Interaction Effect -log10(p) Values (PRSet p < 0.05)'
plot(prset_logp, bann_logp, interaction_title, 'interaction_effects_scatter.png')
plot(prset_logp, bann_logp, interaction_title, 'interaction_effects_scatter_clip5.png', clip=5)

# CHANGE 2: Calculate, sort, and print top interaction effect ratios
interaction_ratios = bann_logp / prset_logp.replace(0, np.nan)
interaction_ratios.index = merged_df['pathway'] # Assign pathway names
top_interaction_ratios = interaction_ratios.sort_values(ascending=False)

print("\n--- Top 10 Pathways with Highest BANN/PRSet Interaction Effect Ratios ---")
print(top_interaction_ratios.head(10))
print("-" * 60)

# CHANGE 2: Calculate and print the percentage for interaction effects
bann_better_interaction = np.sum(merged_df['p.value_bann'] < merged_df['p.value_prset'])
total_interaction = len(merged_df)
percent_bann_better_interaction = (bann_better_interaction / total_interaction) * 100
print(f"\n--- Interaction Effect Significance ---")
print(f"BANN had greater statistical significance in {percent_bann_better_interaction:.2f}% of pathways.")
print(f"Compared a total of {total_interaction} pathways.")
print("-" * 30)


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
