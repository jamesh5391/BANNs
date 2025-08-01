import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import seaborn as sns
from scipy.stats import spearmanr

def plot(x, y, title, out_name, clip=None):
	sns.set_context("talk", font_scale=1.2)
	sns.set_style("whitegrid")
	plt.figure(figsize=(8,8))
	if clip:
		plt.title(f"{title} (Clipped)", fontsize=30, fontweight='bold')
		plt.scatter(x.clip(upper=clip), y.clip(upper=clip), alpha=0.7, linewidths=0.5, s=80, edgecolors='black')
	else:
		plt.title(title, fontsize=20, fontweight='bold') 
		plt.scatter(x, y, alpha=0.7, edgecolors='black', linewidths=0.5, s=80)
	
	plt.xlabel(r'$-\log_{10}(p)$ PRSet', fontsize=20)
	plt.ylabel(r'$-\log_{10}(p)$ BANN', fontsize=20)
	plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=2)
	
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.grid(True, linestyle='--', alpha=0.3)

	sns.despine()

	plt.tight_layout()

	plt.savefig(out_name, dpi=300, bbox_inches='tight')

'''
MAIN EFFECTS
'''

bann_main_df = pd.read_csv('BANNS/BANN/src/main_effect_relu_res_df.csv')
prset_main_df = pd.read_csv('BANNS/BANN/src/main_effects_res_df.csv')
prset_main_df = prset_main_df.iloc[2:]


main_effects_bann = bann_main_df[bann_main_df['term'] == 'ppgs'] 
main_effects_prset = prset_main_df[prset_main_df['term'] == 'ppgs']

#Filter by PRSet p-value < 0.05
main_effects_bann_significant = main_effects_bann[main_effects_bann['p.value'] < 0.05]

# keeps only the pathways that are significant in PRSet AND also exist in the BANN results.
merged_df = pd.merge(main_effects_bann_significant, main_effects_prset, on='pathway', suffixes=('_bann', '_prset'))


main_z_ratios = (merged_df['statistic_bann'] / merged_df['statistic_prset'].replace(0, pd.NA)).dropna()
main_z_ratios.index = merged_df['pathway'] 


bann_logp = -np.log10(merged_df['p.value_bann'].clip(lower=1e-200))
prset_logp = -np.log10(merged_df['p.value_prset'].clip(lower=1e-200))

main_title = 'Main Effect -log(p) Values'
plot(prset_logp, bann_logp, main_title, "main_effects_scatter.png")
plot(prset_logp, bann_logp, main_title, "main_effects_scatter_clip30.png", clip=30)

# Calculate, sort, and print top main effect ratios
main_ratios = bann_logp / prset_logp.replace(0, np.nan)
main_ratios.index = merged_df['pathway'] # Assign pathway names for readability
top_main_ratios = main_ratios.sort_values(ascending=False)

print("--- Top 10 Pathways with Highest BANN/PRSet Main Effect Ratios ---")
print(top_main_ratios.head(20))
print("-" * 60)

bann_better_main = np.sum(merged_df['p.value_bann'] < merged_df['p.value_prset'])
total_main = len(merged_df)
percent_bann_better_main = (bann_better_main / total_main) * 100
print(f"--- Main Effect Significance ---")
print(f"BANN had greater statistical significance in {percent_bann_better_main:.2f}% of pathways.")
print(f"Compared a total of {total_main} pathways.")
print("-" * 30)

print()
merged_all_df = pd.merge(main_effects_bann, main_effects_prset, on='pathway', suffixes=('_bann', '_prset'))
interaction_corr_all = np.corrcoef(merged_all_df['estimate_bann'], merged_all_df['estimate_prset'])
print(f"The correlation for ALL main effect sizes is: {interaction_corr_all[0, 1]:.4f}")
print()

interaction_corr_sig = np.corrcoef(merged_df['estimate_bann'], merged_df['estimate_prset'])
print("Number of signiificant pathways: ", len(main_effects_bann_significant))
print(f"The correlation for SIGNIFICANT main effect sizes is: {interaction_corr_sig[0, 1]:.4f}")
print()



#RESULTS FOR INTERACTION EFFECTS
bann_df = pd.read_csv('BANNS/BANN/src/bann_ppgs_relu_res_df.csv') 
prset_df = pd.read_csv('BANNS/BANN/src/prset_ppgs_res_df.csv')
prset_df = prset_df.iloc[3:]

bann_df = bann_df[bann_df['term'] == 'ppgs:bmi']
prset_df = prset_df[prset_df['term'] == 'ppgs:bmi']


# Filter only the BANN dataframe for p-values < 0.05
bann_df_significant = bann_df[bann_df['p.value'] < 0.05]

# This keeps only the pathways that are significant in PRSet AND also exist in the BANN results.
merged_df = pd.merge(bann_df_significant, prset_df, on='pathway', suffixes=('_bann', '_prset'))


merged_all_df = pd.merge(bann_df, prset_df, on='pathway', suffixes=('_bann', '_prset'))
interaction_corr_all = np.corrcoef(merged_all_df['estimate_bann'], merged_all_df['estimate_prset'])
print(f"The correlation for ALL interaction effect sizes is: {interaction_corr_all[0, 1]:.4f}")

interaction_corr_sig = np.corrcoef(merged_df['estimate_bann'], merged_df['estimate_prset'])
print("Number of signiificant pathways: ", len(bann_df_significant))
print(f"The correlation for SIGNIFICANT interaction effect sizes is: {interaction_corr_sig[0, 1]:.4f}")

# CHANGE 2: Calculate the z-ratios for the interaction effects
interaction_z_ratios = (merged_df['statistic_bann'] / merged_df['statistic_prset'].replace(0, pd.NA)).dropna()
interaction_z_ratios.index = merged_df['pathway'] # Set pathway as index


# CHANGE 3: Combine the ratios and calculate the correlation
# Combine the two ratio Series into a single DataFrame, aligning by pathway
combined_ratios_df = pd.concat([main_z_ratios.rename('main_effect_z_ratio'),
                                interaction_z_ratios.rename('interaction_z_ratio')],
                               axis=1).dropna()

# Calculate the Pearson correlation between the two columns
correlation, p_val = spearmanr(combined_ratios_df['main_effect_z_ratio'], combined_ratios_df['interaction_z_ratio'])
#correlation = combined_ratios_df['main_effect_z_ratio'].spearmanr(combined_ratios_df['interaction_z_ratio'])

print("--- Correlation between Main and Interaction Effect Improvements ---")
#print(f"The Spearman correlation between the z-ratios is: {correlation}")
print(f"The Spearman correlation is: {correlation:.4f}")
print(f"The p-value is: {p_val:.4f}")


print(f"Analysis was based on {len(combined_ratios_df)} overlapping pathways.")
print("-" * 60)

# CHANGE 3: Create logp values from the new merged dataframe
bann_logp = -np.log10(merged_df['p.value_bann'].clip(lower=1e-200))
prset_logp = -np.log10(merged_df['p.value_prset'].clip(lower=1e-200))


interaction_title = 'Interaction Effect -log(p) Values'
plot(prset_logp, bann_logp, interaction_title, 'interaction_effects_scatter.png')
plot(prset_logp, bann_logp, interaction_title, 'interaction_effects_scatter_clip5.png', clip=5)

# CHANGE 2: Calculate, sort, and print top interaction effect ratios
interaction_ratios = bann_logp / prset_logp.replace(0, np.nan)
interaction_ratios.index = merged_df['pathway'] # Assign pathway names
top_interaction_ratios = interaction_ratios.sort_values(ascending=False)

print("\n--- Top 10 Pathways with Highest BANN/PRSet Interaction Effect Ratios ---")
print(top_interaction_ratios.head(20))
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


z_ratios = bann_df['statistic'] / prset_df['statistic'].replace(0, pd.NA)
print(z_ratios.describe())
plotBox(z_ratios.dropna(), "Interaction Z-Ratios", "interaction_boxplot.png")
