#!/bin/bash
#PBS -l mem=10GB
#PBS -l walltime=01:00:00
#PBS -l ncpus=1
#PBS -J 1-11

# This script uses LD score regression to calculate rg for the munged sumstats of each cohort, of each sex
# To ensure all pairwise rg are calculated between sumstats this is a parallel script submitting the list of all sumstats, then from 2nd - last sumstats, then 3rd - last sumstats etc
# e.g. if there are 4 cohort there would be 3 jobs run in parallel:
        # cohort1,cohort2,cohort3,cohort4
        # cohort2,cohort3,cohort4
        # cohort3,cohort4

### Environment ###

module load ldsc/20190815



### Preamble ###

directory=/path/04_LDSC/SNPrg/Female_vs_male_all_cohorts/

reference=/path/LDSC_reference/

sumstats_files=(/path/01_Format_sumstats/Females/Munged_sumstats/*_LDSC_munged.sumstats.gz /path/01_Format_sumstats/Males/Munged_sumstats/*_LDSC_munged.sumstats.gz)

num_traits=${#sumstats_files[@]} # Define the number of files and the total number of traits

# Calculate the range of traits to process based on the job array index
start_index=$((PBS_ARRAY_INDEX - 1))  # Adjust index for 0-based array
end_index=$((PBS_ARRAY_INDEX + num_traits - 2))  # Adjust for 0-based index

# Create the list of traits for the current job
traits=("${sumstats_files[@]:$start_index:$num_traits}")

# Create the comma-separated list of traits for the current job Convert the list to a comma-separated string
traits_list=$(printf "%s\n" "${traits[@]}" | paste -sd, -)

output_file=rg_${PBS_ARRAY_INDEX}


### Submit script ###

cd ${directory}

echo "Processing file: ${sumstats_files[${PBS_ARRAY_INDEX}]} with index ${PBS_ARRAY_INDEX}"

echo "No. traits: ${num_traits}"

echo "Start index: ${start_index}"

echo "End index: ${end_index}"

echo "Traits: ${traits}"

echo "Traits list: ${traits_list}"



# Run LDSC for the current combination of traits
ldsc.py \
  --rg ${traits_list} \
  --ref-ld-chr ${reference}eur_w_ld_chr/ \
  --w-ld-chr ${reference}eur_w_ld_chr/ \
  --out ${output_file}




# Extract results from all files
awk '/Summary of Genetic Correlation Results/,0' rg_1.log | tail -n +2 | head -n -3 | \
    sed -e 's|/path/01_Format_sumstats/[^/]*/Munged_sumstats/||g' \
        -e 's|_sumstats_formatted_formetaanalysis_LDSC_munged.sumstats.gz||g' > ldsc_rg_MDD_all_female_male_cohorts_results.txt

for i in {2..11}
do
    awk '/Summary of Genetic Correlation Results/,0' rg_${i}.log | tail -n +3 | head -n -3 | \
        sed -e 's|/path/01_Format_sumstats/[^/]*/Munged_sumstats/||g' \
            -e 's|_sumstats_formatted_formetaanalysis_LDSC_munged.sumstats.gz||g' >> ldsc_rg_MDD_all_female_male_cohorts_results.txt
done
