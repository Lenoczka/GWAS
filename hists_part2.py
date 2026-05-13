import pandas as pd
import matplotlib.pyplot as plt

# Загружаем данные
het = pd.read_csv('het_results.het', sep='\s+')
# Рассчитываем гетерозиготность
# колонки в файле het_results.het, который создает PLINK:
# O(HOM)	Observed number of homozygous genotypes
# N(NM)	Number of non-missing genotypes
het['HET_rate'] = (het['N(NM)'] - het['O(HOM)']) / het['N(NM)']

# Статистика
mean = het['HET_rate'].mean()
std = het['HET_rate'].std()
lower = mean - 3*std
upper = mean + 3*std

# Находим выбросы
outliers = het[(het['HET_rate'] < lower) | (het['HET_rate'] > upper)]

# Строим гистограмму
plt.figure(figsize=(10,6))
plt.hist(het['HET_rate'], bins=50, color='skyblue', edgecolor='black')
plt.axvline(mean, color='red', linewidth=2, label=f'Mean = {mean:.4f}')
plt.axvline(lower, color='orange', linestyle='--', label=f'Mean-3SD = {lower:.4f}')
plt.axvline(upper, color='orange', linestyle='--', label=f'Mean+3SD = {upper:.4f}')



plt.xlabel('Heterozygosity Rate')
plt.ylabel('Number of Samples')
plt.title('Heterozygosity Distribution')
plt.legend()
plt.savefig('task_1_het.jpg', dpi=300)
plt.close()

if len(outliers) > 0:
    print("\nСписок выбросов:")
    print(outliers[['FID', 'IID', 'HET_rate']].to_string(index=False))
    # Сохраняем
    # Сохраняем с кавычками вокруг FID (как в FAM файле)
    with open('het_outliers.txt', 'w') as f:
        for _, row in outliers.iterrows():
            fid = int(row['FID'])
            iid = int(row['IID'])
            f.write(f'"{fid}"\t{iid}\n')
    print(f"Сохранено {len(outliers)} выбросов")
else:
    print("Выбросов нет")