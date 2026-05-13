import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

maf_results=pd.read_csv("/home/lenoczka/BInf3_ITA.ID_1577890-1/src/maf_results.frq", sep="\s+")

plt.figure(figsize=(10, 6))
plt.hist(maf_results['MAF'], bins=40, color='steelblue', edgecolor='black')
plt.xlabel('Minor Allele Frequency (MAF)')
plt.ylabel('Number of SNPs')
plt.title('Distribution of Minor Allele Frequencies')
plt.savefig('maf.jpg', dpi=300, bbox_inches='tight')
plt.close()
