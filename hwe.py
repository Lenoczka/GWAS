import pandas as pd
import matplotlib.pyplot as plt


hwe_results=pd.read_csv("/home/lenoczka/BInf3_ITA.ID_1577890-1/src/hwe_results.hwe", sep="\s+")

plt.figure(figsize=(10, 6))
plt.hist(hwe_results['P'], bins=40, color='steelblue', edgecolor='black')
plt.xlabel('HWE p-value')
plt.ylabel('Number of SNPs')
plt.title('Distribution of Hardy-Weinberg Equilibrium p-values')
plt.savefig('task_2_hwe.jpg', dpi=300, bbox_inches='tight')
plt.close()