import pandas as pd
import matplotlib.pyplot as plt


results=pd.read_csv("/home/lenoczka/BInf3_ITA.ID_1577890-1/src/genome_results.genome", sep="\s+")

plt.figure(figsize=(10, 6))
plt.hist(results['PI_HAT'], bins=40, color='steelblue', edgecolor='black')
plt.xlabel('PI_HAT (coefficient of relatedness)')
plt.ylabel('Number of pairs')
plt.title('Distribution of Relatedness (PI_HAT) Among All Sample Pairs')
plt.savefig('task_3_pi_hat.jpg', dpi=300, bbox_inches='tight')
plt.close()