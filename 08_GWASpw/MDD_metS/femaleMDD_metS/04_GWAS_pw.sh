#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l mem=10GB
#PBS -l ncpus=1

# This script uses GWAS-pw to identify loci that influence female MDD and metabolic syndrome
# Note this is only for autosomes. Need bed file that includes X chromosome to include SNPs on X chromosome


### Environment ###

module load gwas-pw/0.21


### Preamble ###

directory=/path/06_GWAS_pw/MDD_metS/femaleMDD_metS/ # working directory

reference=/path/gwaspw_reference/ # location of reference directory where .bed file is located

input=femaleMDD_metS_sumstats_formatted_gwaspw.txt # summary stats of two traits combined and formatted for use in gwas-pw

correlation=0.004414508  # correlation of two traits using only SNPs not correlated with either trait (as output from previous step)

pheno1=femaleMDD # name of phenotype 1 = after "_" in input file

pheno2=metS # name of phenotype 2 = after "_" in input file

output=GWAS_pw_femaleMDD_metS # name of gwas-pw output files


### Submit script ###

cd ${directory}

# Run gwas-pw

gwas-pw -i ${input} \
-bed ${reference}EUR/fourier_ls-all.bed \
-phenos ${pheno1} ${pheno2} \
-cor ${correlation} \
-o ${output}



# -cor [float] if the two GWAS were performed using overlapping cohorts, use this flag to specify the expected correlation in summary statistics under the null (defaults to zero)
# use the genome-wide correlation calculated within -cor

###Output file format###



# -[output].segbfs.gz contains a line for each segment of the genome. The columns are:

# chunk: the internal numerical identifer for the segment
# NSNP: the number of SNPs in the segment
# chr: chromosome
# st: star position
# sp: end position
# max_abs_Z_[pheno1]: the maximum absolute value of the Z-score for phenotype 1 in the region
# max_abs_Z_[pheno2]: the maximum absolute value of the Z-score for phenotype 2 in the region
# logBF_1: ln(regional Bayes factor supporting model 1 [association only to phenotype 1] versus the null)
# logBF_2: ln(regional Bayes factor supporting model 2 [association only to phenotype 2] versus the null)
# logBF_3: ln(regional Bayes factor supporting model 3 [shared association to both phenotypes] versus the null)
# logBF_4: ln(regional Bayes factor supporting model 3 [two distinct associations, one to each phenotype] versus the null)
# pi_1: prior on model 1
# pi_2: prior on model 2
# pi_3: prior on model 3
# pi_4: prior on model 4
# PPA_1: posterior probability of model 1
# PPA_2: posterior probability of model 2
# PPA_3: posterior probability of model 3
# PPA_4: posterior probability of model 4

# -[output].bfs.gz contains a line for each SNP in the genome. The columns are:

# id: SNP identifier
# chr: chromosome 3: pos: position
# logBF_1: ln(Bayes factor measure the suppport for model 1 at the SNP)
# logBF_2: ln(Bayes factor measure the suppport for model 2 at the SNP)
# logBF_3: ln(Bayes factor measure the suppport for model 3 at the SNP)
# Z_[pheno1]: Z-score for association to phenotype 1
# V_[pheno1]: variance in the effect size estimate for phenotype 1
# Z_[pheno2]: Z-score for association to phenotype 2
# V_[pheno2]: variance in the effect size estimate for phenotype 2
# pi_1: prior on this SNP being the causal one under model 1
# pi_2: prior on this SNP being the causal one under model 2
# pi_3: prior on this SNP being the causal one under model 3
# PPA_1: posterior probability that this SNP is the causal one under model 1
# PPA_2: posterior probability that this SNP is the causal one under model 2
# PPA_3: posterior probability that this SNP is the causal one under model 3
# chunk: the internal numerical identifer for the segment this SNP falls in

# -[output].MLE contains the estimated regional prior probabilites of each model (same as in [output].segbfs.gz)

