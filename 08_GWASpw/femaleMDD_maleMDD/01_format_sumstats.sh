#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l mem=50GB
#PBS -l ncpus=1

# This script formats the GWAS sum stats so they can be input into fgwas and GWAS-pw

### Environment ###


### Preamble ###

directory=/path/06_GWAS_pw/ # working directory

sumstats1=/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt # first set of summary statistics

sumstats2=/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt # 2nd set of summary statistics

output1=female_MDD_sumstats.txt # name of output based on first set of sum stats

output2=male_MDD_sumstats.txt # name of output based on 2nd set of sum stats

gwaspw_output=female_male_MDD_sumstats_formatted_gwaspw.txt # name of output that will be used for fgwas and gwas-pw

# NOTE: will need to check code is referring to correct column numbers for input files

### Submit script ###

cd ${directory}


## Create Z-score and Variance columns ##

# Z = the signed Z score measuring the evidence for association to phenotype at the SNP = Effect/StdErr
# V = variance = SE^2
# Add "chr" before chromosome number (so matches style in bed file input for fgwas and gwas-pw)

# first file

awk 'BEGIN{OFS="\t"} 
    NR == 1 {
        print $0, "Z", "V";  # Print header and add new columns "Z" and "V"
    } 
    NR > 1 {
        $18 = "chr" $18;  # Prepend "chr" to the CHR column, which is $18
        if ($9 != 0) {
            print $0, $8 / $9, $9 * $9;  # Print original fields and calculated Z and V
        } else {
            print $0, "NaN", "NaN";  # Handle division by zero
        }
    }
' ${sumstats1} > ${output1}


# second file

awk 'BEGIN{OFS="\t"}
    NR == 1 {
        print $0, "Z", "V";  # Print header and add new columns "Z" and "V"
    }
    NR > 1 {
        $18 = "chr" $18;  # Prepend "chr" to the CHR column, which is $18
        if ($9 != 0) {
            print $0, $8 / $9, $9 * $9;  # Print original fields and calculated Z and V
        } else {
            print $0, "NaN", "NaN";  # Handle division by zero
        }
    }
' ${sumstats2} > ${output2}


## Create input files in the format for fgwas: SNPID CHR POS Z SE ##
# Documentation states if include SE 'this column will be taken as a direct estimate of the standard errors of the regression coefficient and will override the sample size and allele frequency columns.'
# However, if include SE and no F or N then get WARNING: detected SE in header, will override F and N and ERROR: cannot find F in header
# So include F = : the allele frequency of one of the alleles of the SNP and N (for case/control studies should provide NCASE and NCONTROL but SE should override this anyway)


awk '
FNR == 1 {
    print "SNPID", "CHR", "POS", "Z", "SE", "F", "N"
}
FNR > 1 {
    print $1, $18, $19, $21, $9, $4, $16
}
' "${output1}" > "${output1%.txt}_formatted_fgwas.txt"

awk '
FNR == 1 {
    print "SNPID", "CHR", "POS", "Z", "SE", "F", "N"
}
FNR > 1 {
    print $1, $18, $19, $21, $9, $4, $16
}
' "${output2}" > "${output2%.txt}_formatted_fgwas.txt"


# Sort rows by chromosomal position
sort -k2,2V -k3,3n ${output1%.txt}_formatted_fgwas.txt > ${output1%.txt}_formatted_fgwas.sorted
sort -k2,2V -k3,3n ${output2%.txt}_formatted_fgwas.txt > ${output2%.txt}_formatted_fgwas.sorted

# Remove chr23 as not in bed file (ideally create bed file including chr23)
awk '($2 != "chr23")' ${output1%.txt}_formatted_fgwas.sorted > ${output1%.txt}_formatted_fgwas_no_chr23.sorted
awk '($2 != "chr23")' ${output2%.txt}_formatted_fgwas.sorted > ${output2%.txt}_formatted_fgwas_no_chr23.sorted


# fgwas throws errors if have multiple rows at the same location (i.e. multiallelic SNPs with one row for each possible SNP at the same location)
# So if there are multiple SNPs at the same location just keep the first SNP (not ideal)
awk '!seen[$2,$3]++' ${output1%.txt}_formatted_fgwas_no_chr23.sorted > ${output1%.txt}_formatted_fgwas_no_chr23_removed_multiallelics.sorted


awk '!seen[$2,$3]++' ${output2%.txt}_formatted_fgwas_no_chr23.sorted > ${output2%.txt}_formatted_fgwas_no_chr23_removed_multiallelics.sorted


## Create input file in the format needed by gwas-pw: SNPID, CHR, POS, Z_[pheno1], V_[pheno1], Z_[pheno2], V_[pheno2] ##

# $1 = MarkerName
# $21 = Z-score
# $22 = Variance
# $18 = CHR
# $19 = BP
# array femalemarker: key = $1 MarkerName, value = = $1 MarkerName
# array femaleZ: key = $1 MarkerName, value = $21 Z-score
# array femaleV: key = $1 MarkerName, value = $22 Variance
# array femaleCHR: key = $1 MarkerName, value = $18 CHR
# array femaleV: key = $1 MarkerName, value = $19 BP


awk '
# First pass: read first file and store data in arrays
NR==FNR {
    femalemarker[$1]=$1;
    femaleZ[$1]=$21;
    femaleV[$1]=$22;
    femaleCHR[$1] = $18;
    femaleBP[$1] = $19;
    seenInFemale[$1]=1;
    next
}

# Second pass: process second file
{
    if (FNR == 1) {
        print "SNPID", "CHR", "POS", "Z_femaleMDD", "V_femaleMDD", "Z_maleMDD", "V_maleMDD";
    } else {
        if ($1 in femalemarker) {
            print $1, $18, $19, femaleZ[$1], femaleV[$1], $21, $22;
            seenInMale[$1] = 1;
        } else {
            print $1, $18, $19, "NA", "NA", $21, $22;
        }
    }
}

# After processing second file (male), print remaining female-only records
END {
    for (marker in femalemarker) {
	if (marker != "MarkerName" && !(marker in seenInMale)) {
            print marker, femaleCHR[marker], femaleBP[marker], femaleZ[marker], femaleV[marker], "NA", "NA";
        }
    }
}
' ${output1} ${output2} > ${gwaspw_output}


# Sort rows by chromosomal position
sort -k2,2V -k3,3n ${gwaspw_output} > ${gwaspw_output%.txt}.sorted


# Remove chr23 as not in bed file (ideally create bed file including chr23)
awk '($2 != "chr23")' ${gwaspw_output%.txt}.sorted > ${gwaspw_output%.txt}_no_chr23.sorted

# gwas-pw throws errors if have multiple rows at the same location (i.e. multiallelic SNPs with one row for each possible SNP at the same location)
# So if there are multiple SNPs at the same location just keep the first SNP (not ideal)
awk '!seen[$2,$3]++' ${gwaspw_output%.txt}_no_chr23.sorted > ${gwaspw_output%.txt}_no_chr23_removed_multiallelics.sorted

