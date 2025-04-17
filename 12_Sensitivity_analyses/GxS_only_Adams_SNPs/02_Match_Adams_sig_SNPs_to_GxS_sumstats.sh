#!/bin/bash
#PBS -l walltime=01:30:00
#PBS -l mem=150GB
#PBS -l ncpus=1

# This script matches the Adams significant SNPs to my GxS meta-analysis sumststs file (both use build 37)
# Match on rsID, A1, A2 (allowing alleles to be switched) for more accurate SNP matching (just rsID won't account for multiallelic sites)


### Preamble ###

directory=/path/08_Sensitivity_analyses/Adams2024_topSNPs/

GxS=/path/03_Metal/GxS/fulldc/Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID.txt

Adams=Adams_Full_EUR_2025_FromIMB_sig_SNPs.txt

output=Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID_Adams_sig_SNPs.txt


### Run script ###

cd ${directory}


awk 'NR==FNR {
        marker[$1":"tolower($2)":"tolower($3)] = $1":"tolower($2)":"tolower($3);
        next
    }
    {
        if (FNR == 1) {
            print $0;
        } else if (($20":"$2":"$3) in marker) {
            # Check for exact match with composite key
            print $0;
        } else if (($20":"$3":"$2) in marker) {
            # Check for reverse match of $3 and $2
            print $0;
        }
    }' ${Adams} ${GxS} > ${output}
