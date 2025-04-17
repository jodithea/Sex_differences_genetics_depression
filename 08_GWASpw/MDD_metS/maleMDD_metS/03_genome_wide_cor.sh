#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l mem=30GB
#PBS -l ncpus=1

# This script is to determine the correlation in effect sizes between the two traits, only using SNPs that are not associated with either trait
# (i.e. SNPs from segments with a posterior probability of association < 0.2 in both traits)
# This will be input as the correlation in gwas-pw
# (gwas-pw: -cor [float] if the two GWAS were performed using overlapping cohorts, use this flag to specify the expected correlation in summary statistics under the null (defaults to zero))


### Environment ###

module load gzip/1.13

module load R/4.3.1


### Preamble ###

directory=/path/06_GWAS_pw/MDD_metS/maleMDD_metS/ # working directory

input1_seg=maleMDD.segbfs.gz # fgwas output .segbfs.gz file from 1st set of summary stats 

input1_snp=maleMDD.bfs.gz # fgwas output .bfs file from 1st set of summary stats

input1_sumstats=maleMDD_sumstats_matched_distinct.txt # 1st set of summary stats

input2_seg=metS.segbfs.gz # fgwas output .segbfs.gz file from 2nd set of summary stats

input2_snp=metS.bfs.gz # fgwas output .bfs file from 2nd set of summary stats

input2_sumstats=metS_sumstats_matched_distinct.txt # 2nd set of summary stats

output=maleMDD_metS # name of output file that will print correlation between the two traits



### Run script ###

cd ${directory}

# Unzip fgwas output files
gzip -dk ${input1_seg} ${input1_snp} ${input2_seg} ${input2_snp}


# For each trait, retain segments with ppa < 0.2 (.segbfs file)
awk '
FNR==1 {
  print
  next
}
($10 < 0.2)
' ${input1_seg%.gz} > ${input1_seg%.segbfs.gz}_PPA02.segbfs


awk '
FNR==1 {
  print
  next
}
($10 < 0.2)
' ${input2_seg%.gz} > ${input2_seg%.segbfs.gz}_PPA02.segbfs


# For each trait, retain SNPs with segment ppa < 0.2 (.bfs file)
# merge list of segments with PPA < 0.2 (.segbfs) with SNP file (.bfs) and only keep those SNPS that are in the segments with ppa < 0.2
# keeps chunks in .bfs that are also in .segbfs
awk 'FNR==NR{a[$1]; next} ($11) in a' ${input1_seg%.segbfs.gz}_PPA02.segbfs ${input1_snp%.gz} > ${input1_snp%.bfs.gz}_PPA02.bfs
awk 'FNR==NR{a[$1]; next} ($11) in a' ${input2_seg%.segbfs.gz}_PPA02.segbfs ${input2_snp%.gz} > ${input2_snp%.bfs.gz}_PPA02.bfs


# Merge .bfs files restricted to SNPS with segment ppa < 0.2 for both traits, retaining SNPs from segments with ppa < 0.2 for both traits
# keeps SNPs from input2 that are also in input1
awk 'FNR==NR{a[$1]; next} ($1) in a' ${input1_snp%.bfs.gz}_PPA02.bfs ${input2_snp%.bfs.gz}_PPA02.bfs > ${output}_PPA02.bfs



# Extract effect sizes for this list of SNPs for each trait
# merge list of SNPs from segments with ppa < 0.2 for both traits with each trait sum stats and only retain sum stats for those SNPs
# keeps snps in MDD summary stats (or metS summary stats) that are also in MDD/metS combined bfs file (keeping header)
awk '
FNR == NR {
    a[$1]
    next
}
FNR == 1 {
    print
    next
}
($1) in a
' "${output}_PPA02.bfs" "${input1_sumstats}" > "${input1_sumstats%.txt}_PPA02.txt"


awk '
FNR == NR {
    a[$1]
    next
}
FNR == 1 {
    print
    next
}
($1) in a
' "${output}_PPA02.bfs" "${input2_sumstats}" > "${input2_sumstats%.txt}_PPA02.txt"



# Submit extracted effect sizes of SNPs from segments with ppa < 0.2 in both traits to R script to calculate correlation in effect sizes between the 2 traits
# This is the correlation between the two traits of SNPs that are not associated with either trait
Rscript --vanilla 03_genome_wide_cor.R ${input1_sumstats%.txt}_PPA02.txt ${input2_sumstats%.txt}_PPA02.txt ${output}




