#!/bin/bash
#PBS -l mem=10GB
#PBS -l walltime=01:00:00
#PBS -l ncpus=1

# This script uses LD score regression to calculate rg for the munged sumstats of female meta-analysis results with other traits
# Uses --print-delete-vals flag to get delete values of genetic covariance (systematically leaving out one chunk of data at a time)


### Environment ###

module load ldsc/20200724



### Preamble ###

directory=/path/04_LDSC/SNPrg/Females/

reference=/path/LDSC_reference/

Allcohorts_female_MDD=/path/04_LDSC/Munge/Females/Metaanalysis_MDD_female_AllCohorts_LDSC_munged.sumstats.gz

Allcohorts_male_MDD=/path/04_LDSC/Munge/Males/Metaanalysis_MDD_male_AllCohorts_LDSC_munged.sumstats.gz

AdamsMDD_sumstats=/path/Summary_statistics_published/Depression_Adams_2024/including_23andMe/Adams_Full_EUR_2025_FromIMB_munged.sumstats.gz

BloklandMDDfemale_sumstats=/path/Summary_statistics_published/Depression_sex_Blokland_2021/meta_STDERR_mdd_eur_auto_xchr_F1_sexstratified_06_gcOFF_pgc+ipsych_munged.sumstats.gz

BloklandMDDmale_sumstats=/path/Summary_statistics_published/Depression_sex_Blokland_2021/meta_STDERR_mdd_eur_auto_xchr_M1_sexstratified_06_gcOFF_pgc+ipsych_munged.sumstats.gz

SilveiraMDDfemale_sumstats=/path/Summary_statistics_published/Depression_sex_Silveira_2023/GWAS.BroadMDD_Females_munged.sumstats.gz

SilveiraMDDmale_sumstats=/path/Summary_statistics_published/Depression_sex_Silveira_2023/GWAS.BroadMDD_Males_munged.sumstats.gz

ADHD_sumstats=/path/Summary_statistics_published/ADHD_Demontis_2023/ADHD2022_iPSYCH_deCODE_PGC_munged.sumstats.gz

Anxiety_sumstats=/path/Summary_statistics_published/Anxiety_Friligkou_2024/ANX_EUR.txt.gz_munged.sumstats.gz

Bipolar_sumstats=/path/Summary_statistics_published/Bipolar_Mullins_2021/Bipolar_Neff_munged.sumstats.gz

PTSD_sumstats=/path/Summary_statistics_published/PTSD_Nievergelt_2024/PTSD_Neff_munged.sumstats.gz

Schizophrenia_sumstats=/path/Summary_statistics_published/Schizophrenia_Trubetskoy_2022/Schizophrenia_Neff_munged.sumstats.gz

BMI_sumstats=/path/Summary_statistics_published/BMI_Yengo_2018/Meta-analysis_Locke_et_al_plus_UKBiobank_2018_UPDATED_munged.sumstats.gz

BMI_female_sumstats=/path/Summary_statistics_published/BMI_sex_Pulit_2018/bmi.giant-ukbb.meta-analysis.females.23May2018_munged.sumstats.gz

BMI_male_sumstats=/path/Summary_statistics_published/BMI_sex_Pulit_2018/bmi.giant-ukbb.meta-analysis.males.23May2018_munged.sumstats.gz

Waist_to_hip_ratio_sumstats=/path/Summary_statistics_published/waist_to_hip_ratio_adjusted_for_BMI_Pulit_2019/fat-distn.giant.ukbb.meta-analysis.whradjbmi.combined_munged.sumstats.gz

Metabolic_syndrome_sumstats=/path/Summary_statistics_published/Metabolic_syndrome_Park_2024/GCST90444487_munged.sumstats.gz

EA_sumstats=/path/Summary_statistics_published/EA4_Okbay_2022/EA4_additive_excl_23andMe_munged.sumstats.gz

Drinks_per_week_sumstats=/path/Summary_statistics_published/Substance_use_Saunders_2022/GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR_munged.sumstats.gz

Ever_smoked_regularly_sumstats=/path/Summary_statistics_published/Substance_use_Saunders_2022/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_munged.sumstats.gz

output=rg_femaleMDD



### Submit script ###

cd ${directory}

# Run LDSC for the current combination of traits
ldsc.py \
  --rg ${Allcohorts_female_MDD},${Allcohorts_male_MDD},${AdamsMDD_sumstats},${BloklandMDDfemale_sumstats},${BloklandMDDmale_sumstats},${SilveiraMDDfemale_sumstats},${SilveiraMDDmale_sumstats},${ADHD_sumstats},${Anxiety_sumstats},${Bipolar_sumstats},${PTSD_sumstats},${Schizophrenia_sumstats},${BMI_sumstats},${Waist_to_hip_ratio_sumstats},${Metabolic_syndrome_sumstats},${EA_sumstats},${BMI_female_sumstats},${BMI_male_sumstats},${Drinks_per_week_sumstats},${Ever_smoked_regularly_sumstats} \
  --ref-ld-chr ${reference}eur_w_ld_chr/ \
  --w-ld-chr ${reference}eur_w_ld_chr/ \
  --print-delete-vals \
  --out ${output}




# Output table of results in a new file

# Extract the relevant part of the log
awk '/Summary of/{flag=1; next} /Analysis finished/{flag=0} flag' ${output}.log > ${output}_correlation_table.txt

# Define the identifiers and their variable names
declare -A identifiers=(
["$Allcohorts_female_MDD"]="Allcohorts_female_MDD"
["$Allcohorts_male_MDD"]="Allcohorts_male_MDD"
["$AdamsMDD_sumstats"]="AdamsMDD_sumstats"
["$BloklandMDDfemale_sumstats"]="BloklandMDDfemale_sumstats"
["$BloklandMDDmale_sumstats"]="BloklandMDDmale_sumstats"
["$SilveiraMDDfemale_sumstats"]="SilveiraMDDfemale_sumstats"
["$SilveiraMDDmale_sumstats"]="SilveiraMDDmale_sumstats"
["$ADHD_sumstats"]="ADHD_sumstats"
["$Anxiety_sumstats"]="Anxiety_sumstats"
["$Bipolar_sumstats"]="Bipolar_sumstats"
["$PTSD_sumstats"]="PTSD_sumstats"
["$Schizophrenia_sumstats"]="Schizophrenia_sumstats"
["$BMI_sumstats"]="BMI_sumstats"
["$Waist_to_hip_ratio_sumstats"]="Waist_to_hip_ratio_sumstats"
["$Metabolic_syndrome_sumstats"]="Metabolic_syndrome_sumstats"
["$EA_sumstats"]="EA_sumstats"
["$BMI_female_sumstats"]="BMI_female_sumstats"
["$BMI_male_sumstats"]="BMI_male_sumstats"
["$Drinks_per_week_sumstats"]="Drinks_per_week_sumstats"
["$Ever_smoked_regularly_sumstats"]="Ever_smoked_regularly_sumstats"
)

# Function to escape special characters for awk gsub (especially because Blokland sumstats file name includes '+')
escape_special_chars() {
    local value="$1"
    echo "$value" | sed -e 's/[\/&]/\\&/g' -e 's/+/\\+/g'
}

# Generate the replacement part of the awk script dynamically
replacement_script=""
for id in "${!identifiers[@]}"; do
    var_value="$id"
    escaped_value=$(escape_special_chars "$var_value")
    var_name="${identifiers[$id]}"
    replacement_script+="gsub(/${escaped_value}/, \"${var_name}\"); "
done

# Create a dynamic awk script
awk_script=$(cat <<EOF
BEGIN { OFS="\t" }
{
    $replacement_script
    print
}
EOF
)

# Run the dynamically created awk script
awk "$awk_script" "${output}_correlation_table.txt" > temp.txt && mv temp.txt "${output}_correlation_table.txt"
