#!/bin/bash

INPUT_BFILE="/home/lenoczka/BInf3_ITA.ID_1577890-1/datasets/BInf3_ITA/penncath"
BASE_DIR="/home/lenoczka/BInf3_ITA.ID_1577890-1/src"

mkdir -p ${BASE_DIR}/QC_missing_geno_mind \
         ${BASE_DIR}/QC_heterozygosity \
         ${BASE_DIR}/QC_maf_hwe \
         ${BASE_DIR}/relatedness \
         ${BASE_DIR}/pca_stratification \
         ${BASE_DIR}/gwas_logistic \
         ${BASE_DIR}/prs_calculation

plink --bfile ${INPUT_BFILE} --missing --out ${BASE_DIR}/QC_missing_geno_mind/miss_output
python hists_part1.py

plink --bfile ${INPUT_BFILE} --geno 0.1 --mind 0.1 --make-bed --out ${BASE_DIR}/QC_missing_geno_mind/filtered_output

plink --bfile ${BASE_DIR}/QC_missing_geno_mind/filtered_output --het --out ${BASE_DIR}/QC_heterozygosity/het_results
# python script -> het_outliers.txt
python hists_part2.py

plink --bfile ${BASE_DIR}/QC_missing_geno_mind/filtered_output --remove ${BASE_DIR}/QC_heterozygosity/het_outliers.txt --make-bed --out ${BASE_DIR}/QC_heterozygosity/filtered_data_no_outliers

plink --bfile ${BASE_DIR}/QC_heterozygosity/filtered_data_no_outliers --freq --out ${BASE_DIR}/QC_maf_hwe/maf_results
python  maf.py

plink --bfile ${BASE_DIR}/QC_heterozygosity/filtered_data_no_outliers --maf 0.01 --make-bed --out ${BASE_DIR}/QC_maf_hwe/filtered_data_maf_filtered

plink --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_maf_filtered --hardy --out ${BASE_DIR}/QC_maf_hwe/hwe_results
python hwe.py

plink --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_maf_filtered --hwe 1e-6 --make-bed --out ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered

plink --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered --genome --out ${BASE_DIR}/relatedness/genome_results
python pi_hat.py

plink --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered --genome --min 0.2 --out ${BASE_DIR}/relatedness/genome_results_min_filtered

if [ -f ${BASE_DIR}/relatedness/genome_results_min_filtered.genome ]; then
    grep -v "FID" ${BASE_DIR}/relatedness/genome_results_min_filtered.genome | awk '{print $3, $4}' > ${BASE_DIR}/relatedness/samples_to_remove.txt
    plink --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered --remove ${BASE_DIR}/relatedness/samples_to_remove.txt --make-bed --out ${BASE_DIR}/relatedness/genome_data_filtered
    PCA_BFILE="${BASE_DIR}/relatedness/genome_data_filtered"
else
    PCA_BFILE="${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered"
fi

plink19 --bfile ${PCA_BFILE} --pca --out ${BASE_DIR}/pca_stratification/pca_results
python pca12.py

plink19 --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered --logistic --covar ${BASE_DIR}/pca_stratification/pca_results.eigenvec --covar-number 1-2 --out ${BASE_DIR}/gwas_logistic/gwas_results_corrected

plink19 --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered --logistic sex --covar ${BASE_DIR}/pca_stratification/pca_results.eigenvec --covar-number 1-2 --out ${BASE_DIR}/gwas_logistic/gwas_results_corrected_sex

awk 'NR==1 || ($9 < 5e-8 && NR>1)' ${BASE_DIR}/gwas_logistic/gwas_results_corrected_sex.assoc.logistic > ${BASE_DIR}/gwas_logistic/task4_significant.txt

awk 'NR>1 && $9 < 5e-8 {print $2}' ${BASE_DIR}/gwas_logistic/gwas_results_corrected_sex.assoc.logistic > ${BASE_DIR}/gwas_logistic/task4_snplist.txt

if [ -f "./draw_plots.py" ]; then
    ./draw_plots.py ${BASE_DIR}/gwas_logistic/gwas_results_corrected_sex.assoc.logistic
fi

# online FUMA aanotation of GWAS results FG_R12_CAD_remapped -> interpretation.txt

if [ -f "./PGS004595.txt" ]; then
    plink19 --bfile ${BASE_DIR}/QC_maf_hwe/filtered_data_hwe_filtered --score PGS004595.txt 1 4 5 --out ${BASE_DIR}/prs_calculation/prs_results
    head -5 ${BASE_DIR}/prs_calculation/prs_results.profile
fi

echo "Pipeline completed successfully"

python prs_stats.py
