import pandas as pd
import matplotlib.pyplot as plt

# Загружаем PCA
pca = pd.read_csv('pca_results.eigenvec', sep='\s+', header=None)
pca = pca[[0, 1, 2, 3]]
pca.columns = ['FID', 'IID', 'PC1', 'PC2']

# Загружаем собственные значения
eigenval = pd.read_csv('pca_results.eigenval', header=None)
eigenval.columns = ['Eigenvalue']

# Рассчитываем % дисперсии
total = eigenval['Eigenvalue'].sum()
var_pc1 = eigenval['Eigenvalue'].iloc[0] / total * 100
var_pc2 = eigenval['Eigenvalue'].iloc[1] / total * 100

# Загружаем фенотипы
fam = pd.read_csv('filtered_data_hwe_filtered.fam', sep='\s+', header=None, usecols=[0,1,5])
fam.columns = ['FID', 'IID', 'PHENOTYPE']

# Объединяем
pca = pca.merge(fam, on=['FID', 'IID'])

# Разделяем на группы по фенотипу
controls = pca[pca['PHENOTYPE'] == 1]    # Контроли (синий)
cases = pca[pca['PHENOTYPE'] == 2]       # Случаи (красный)
unknown = pca[pca['PHENOTYPE'] == -9]    # Неизвестно (серый)

# График с легендой
plt.figure(figsize=(10, 8))

# Рисуем каждую группу отдельно
plt.scatter(controls['PC1'], controls['PC2'],
            c='blue', s=30, alpha=0.7, edgecolors='black', linewidth=0.5,
            label='Controls')
plt.scatter(cases['PC1'], cases['PC2'],
            c='red', s=30, alpha=0.7, edgecolors='black', linewidth=0.5,
            label='Cases')
plt.scatter(unknown['PC1'], unknown['PC2'],
            c='gray', s=30, alpha=0.7, edgecolors='black', linewidth=0.5,
            label='Unknown')

# Подписи осей с % дисперсии
plt.xlabel(f'Principal Component 1 ({var_pc1:.2f}%)')
plt.ylabel(f'Principal Component 2 ({var_pc2:.2f}%)')

plt.title('PCA of Genetic Data')
plt.legend()  # Добавляем легенду
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('task_3_pca12.jpg', dpi=300)
plt.close()
