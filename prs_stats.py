import pandas as pd
from scipy.stats import mannwhitneyu

df = pd.read_csv('/home/lenoczka/BInf3_ITA.ID_1577890-1/src/prs_results.profile', sep=r'\s+')

cases = df[df['PHENO'] == 1]['SCORE'].dropna()
controls = df[df['PHENO'] == 2]['SCORE'].dropna()

stat, p_value = mannwhitneyu(cases, controls)

print(f"p-value = {p_value}")

with open('/home/lenoczka/BInf3_ITA.ID_1577890-1/src/task6_answers.txt', 'a') as f: # a от append -  дозапись текста в файл
    f.write(f"p = {p_value}\n")
