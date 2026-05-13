import pandas as pd
import matplotlib.pyplot as plt

snp= pd.read_csv("/home/lenoczka/BInf3_ITA.ID_1577890-1/src/miss_output.lmiss", sep="\s+") #\s+	один или несколько пробелов/табуляций подряд
patient= pd.read_csv("/home/lenoczka/BInf3_ITA.ID_1577890-1/src/miss_output.imiss", sep="\s+")

# Гистограмма по SNP
plt.figure(figsize=(10, 6))
plt.hist(snp['F_MISS'], bins=40, color='steelblue', edgecolor='black')
plt.xlabel('Proportion of missing genotypes (F_MISS)')
plt.ylabel('Number of SNPs')
plt.title('Distribution of missing genotypes by SNP')
plt.savefig('SNP_miss.jpg', dpi=300, bbox_inches='tight')
plt.close()

# Гистограмма по индивидам
plt.figure(figsize=(10, 6))
plt.hist(patient['F_MISS'], bins=40, color='steelblue', edgecolor='black')
plt.xlabel('Proportion of missing genotypes (F_MISS)')
plt.ylabel('Number of individuals')
plt.title('Distribution of missing genotypes by individual')
plt.savefig('individual_miss.jpg', dpi=300, bbox_inches='tight')
plt.close()
