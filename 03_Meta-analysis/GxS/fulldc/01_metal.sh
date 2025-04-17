#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l mem=64GB
#PBS -l ncpus=1
#PBS -J 0-6

# This script uses METAL to run a meta-analysis on multiple GWAS summary statistics
# It runs in parallel a meta-analysis for all cohorts, plus all possible LOO meta-analyses

### Environment ###

module load metal/20200505 


### Preanble ###

directory=/path/03_Metal/GxS/fulldc/

files=(/path/01_Format_sumstats/GxS/fulldc/*_MDD_GxS_fulldc_sumstats_formatted_formetaanalysis_transformed.txt)



### Run script ###

cd ${directory}

# Create the METAL configuration file dynamically for each job

# If PBS array is 0, create the config for the analysis with all cohorts
if [ $PBS_ARRAY_INDEX -eq 0 ]; then
    # Create the file for all cohorts
    config_file="01_metal_all.txt"
    outfile="Metaanalysis_MDD_GxS_fulldc_AllCohorts_"

    echo "Creating config for all cohorts..."
    {
	echo "SCHEME STDERR"
	echo "AVERAGEFREQ ON"
	echo "MINMAXFREQ ON"
	echo "GENOMICCONTROL OFF"
	echo "VERBOSE OFF"
	echo "COLUMNCOUNTING STRICT"
	echo "CUSTOMVARIABLE N_TOTAL"
	echo "CUSTOMVARIABLE NUM_STUDIES_METAANALYSED"
        echo "MARKER MARKER_build37"
        echo "ALLELE A1 A2"
        echo "EFFECT BETA_G_by_E_logistic"
        echo "PVALUE P_G_by_E"
        echo "WEIGHT N"
        echo "STDERR SE_of_BETA_G_by_E_logistic"
        echo "FREQLABEL FREQA1"
        echo "LABEL N_TOTAL as N"
        echo "LABEL NUM_STUDIES_METAANALYSED as STUDY"

        # Process all files
        for file in "${files[@]}"; do
            echo "PROCESS ${file}"
        done

        echo "OUTFILE ${outfile} .tbl"
        echo "ANALYZE HETEROGENEITY"
        echo "QUIT"
    } > $config_file

    # Run METAL with the created config file
    metal $config_file

else
    # Leave-one-out analysis
    leave_out=$(basename ${files[$PBS_ARRAY_INDEX - 1]} | sed 's|_MDD_GxS_fulldc_sumstats_formatted_formetaanalysis_transformed.txt||')
    config_file="01_metal_leaveout_${leave_out}.txt"
    outfile="Metaanalysis_MDD_GxS_fulldc_leaveout_${leave_out}_"

    echo "Creating config leaving out ${leave_out}..."
    {
	echo "SCHEME STDERR"
        echo "AVERAGEFREQ ON"
        echo "MINMAXFREQ ON"
        echo "GENOMICCONTROL OFF"
        echo "VERBOSE OFF"
        echo "COLUMNCOUNTING STRICT"
        echo "CUSTOMVARIABLE N_TOTAL"
        echo "CUSTOMVARIABLE NUM_STUDIES_METAANALYSED"
        echo "MARKER MARKER_build37"
        echo "ALLELE A1 A2"
        echo "EFFECT BETA_G_by_E_logistic"
        echo "PVALUE P_G_by_E"
        echo "WEIGHT N"
        echo "STDERR SE_of_BETA_G_by_E_logistic"
        echo "FREQLABEL FREQA1"
        echo "LABEL N_TOTAL as N"
        echo "LABEL NUM_STUDIES_METAANALYSED as STUDY"

        # Process all files except the one to leave out
	leave_out_index=$((PBS_ARRAY_INDEX - 1))

        for i in "${!files[@]}"; do
            if [ $i -ne $leave_out_index ]; then
                echo "PROCESS ${files[$i]}"
            fi
        done

        echo "OUTFILE ${outfile} .tbl"
        echo "ANALYZE HETEROGENEITY"
        echo "QUIT"
    } > $config_file

    # Run METAL with the created config file
    metal $config_file
fi

echo "Finished running job ${PBS_ARRAY_INDEX}"
