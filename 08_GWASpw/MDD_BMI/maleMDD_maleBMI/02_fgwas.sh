#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=10GB
#PBS -l ncpus=1

# This script uses fgwas to calculate a posterior probability of association (ppa) for each SNP and each segment
# Note this is only for autosomes. Need bed file that includes X chromosome to include SNPs on X chromosome


### Environment ###

module load fgwas/0.3.6


### Preamble ###

directory=/path/06_GWAS_pw/MDD_BMI/maleMDD_maleBMI_update/

reference=/path/gwaspw_reference/ # location of reference directory where .bed file is located

input1=male_MDD_sumstats_formatted_fgwas.txt # 1st set of sum stats that are formatted for fgwas

input2=male_BMI_sumstats_formatted_fgwas.txt # 2nd set of sum stats that are formatted for fgwas

output1=maleMDD # name of fgwas output based on 1st set of sum stats

output2=maleBMI # name of fgwas output based on 2nd set of sum stats


### Submit script ###

cd ${directory}

# Run fgwas for first file

fgwas -i ${input1} \
 -bed ${reference}EUR/fourier_ls-all.bed \
 -print \
 -o ${output1}


# Run fgwas for second file

fgwas -i ${input2} \
 -bed ${reference}EUR/fourier_ls-all.bed \
 -print \
 -o ${output2}


# -print: to re-weight GWAS and output posterior probabilities of association
# -cc: when inputting summary statistics for a binary trait (case/control)
# documentation states 'if you include a column labeled SE in the header of the file, this column will be taken as a direct estimate of the standard errors of
# the log-odds ratio and will override the sample size and allele frequency columns'
# However, if only include SE get WARNING: detected SE in header, will override F and N and ERROR: cannot find F in header
# So included F and N to not get this error
# If include -cc flag get ERROR: cannot find NCASE in header (case/control should have NCASE and NCONTROL whereas quatitative just needs N)
# But if run without -cc flag and provide SE, F and N get WARNING: detected SE in header, will override F and N
# So it appears that -cc is just telling fgwas what headers to look for, and as we're using SE instead of F and N (NCASE and NCONTROL) this doesn't matter



# output: a .bfs and a .segbfs for each trait. 
# They will have a posterior probability of association (ppa) for each SNP (.bfs) or each segment (.segbfs). 
# segment PPA = the posterior probability that the segment contains a SNP associated with the trait
# SNP PPA = the posterior probability that the SNP is causal
