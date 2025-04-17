#!/bin/bash
#PBS -l walltime=01:30:00
#PBS -l mem=150GB
#PBS -l ncpus=1
#PBS -J 0-6


# This script performs QC and formats the meta-analysis output files.


### Preamble ###

directory=/path/03_Metal/Females/

files=(/path/03_Metal/Females/Metaanalysis_MDD_female_*_1.tbl)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/03_Metal/Females/||' | sed 's|_1.tbl||')

dbSNP=/path/Meta-analysis/dbSNP155_GRCh37.p13_vs_dbSNP155_GRCh38.p13/dbSNP_MAF0_0001/NCBI_Hsapiens_dbSNP155_GRCh37_GRCh38.p13_split_multiallelics_MAF0_0001_ALL_matches.txt


### Run script ###

cd ${directory}

## QC ##
# Remove if SNP was not in 3+ studies
awk '$17 >= 3 {print}' ${files[${PBS_ARRAY_INDEX}]} > ${name}_QCed.txt


## Format ##

# Add columns CHR and BP to file - so can make Manhattan plot
awk '{ if (NR==1) 
	print $0,"CHR","BP"; 
	else {split($1, arr, ":"); print $0,arr[1],arr[2] } 
}' ${name}_QCed.txt > ${name}_temp.txt && mv ${name}_temp.txt ${name}_QCed.txt

# Change in CHR column X to 23 - so can make Manhattan plot
awk '{ gsub("X", "23", $18); print }' ${name}_QCed.txt > ${name}_temp.txt && mv ${name}_temp.txt ${name}_QCed.txt


# Change column P-value to P - so R script to make Manhattan plot works
awk '{ if (NR==1) 
	{ sub("P-value", "P", $10); print } 
	else print }' ${name}_QCed.txt > ${name}_temp.txt && mv ${name}_temp.txt ${name}_QCed.txt


# Add rsID from reference dbSNP file (all CHR, BP are in build 37 so match to build 37 columns and output rsID_build37)
awk 'NR==FNR { 
	rsid[$4]=$1;
	next }

     FNR==1 { 
	print $0, "rsID_build37"; 
	next }
{
    if ($1 in rsid)
        print $0, rsid[$1];
    else
        print $0, "NA";
}' ${dbSNP} ${name}_QCed.txt > ${name}_QCed_rsID.txt

## Checks ##

echo "Number SNPs missing an rsID:" >> ${name}_check_QCed_rsID_file.txt
tail -n +2 ${files[${PBS_ARRAY_INDEX}]} | grep -c "NA" >> ${name}_check_QCed_rsID_file.txt
echo "" >> ${name}_check_QCed_rsID_file.txt

echo "Number of rsIDs with > 1 lines = multiallelic SNPS:" >> ${name}_check_QCed_rsID_file.txt
awk '{print $20}' ${name}_QCed_rsID.txt | sort | uniq -c | awk '$1 > 1 {print}' | wc -l >> ${name}_check_QCed_rsID_file.txt
echo "" >> ${name}_check_QCed_rsID_file.txt

echo "Number of SNPs with significant P and significant hetergeneity P (both < 5e-08):" >> ${name}_check_QCed_rsID_file.txt
awk '($10 < 5e-08) && ($15 < 5e-08)' ${name}_QCed_rsID.txt | wc -l >> ${name}_check_QCed_rsID_file.txt
echo "" >> ${name}_check_QCed_rsID_file.txt


# Before QC: Check min and max of each column to see if anything weird
echo "Min and max of each column from METAL output" >> ${name}_check_QCed_rsID_file.txt
awk 'NR==1 {
        for (i=1; i<=NF; i++)
            headers[i] = $i  # Store column headers in an array
     }
     NR > 1 {
        for (i=1; i<=NF; i++) {
            if (NR == 2 || $i < min[i]) {  # Update min value for each column
                min[i] = $i
                minHeader[i] = headers[i]  # Store corresponding header
            }
            if (NR == 2 || $i > max[i]) {  # Update max value for each column
                max[i] = $i
                maxHeader[i] = headers[i]  # Store corresponding header
            }
        }
     }
     END {
        for (i=1; i<=NF; i++) {
            print "Min value for column " headers[i] ": " min[i]
            print "Max value for column " headers[i] ": " max[i]
        }
     }' ${files[${PBS_ARRAY_INDEX}]} >> ${name}_check_QCed_rsID_file.txt
echo "" >> ${name}_check_QCed_rsID_file.txt



# After QC and formatting: Check min and max of each column to see if anything weird

echo "Min and max of each column after QC and formatting" >> ${name}_check_QCed_rsID_file.txt
awk 'NR==1 {
        for (i=1; i<=NF; i++)
            headers[i] = $i  # Store column headers in an array
     }
     NR > 1 {
        for (i=1; i<=NF; i++) {
            if (NR == 2 || $i < min[i]) {  # Update min value for each column
                min[i] = $i
                minHeader[i] = headers[i]  # Store corresponding header
            }
            if (NR == 2 || $i > max[i]) {  # Update max value for each column
                max[i] = $i
                maxHeader[i] = headers[i]  # Store corresponding header
            }
        }
     }
     END {
        for (i=1; i<=NF; i++) {
            print "Min value for column " headers[i] ": " min[i]
            print "Max value for column " headers[i] ": " max[i]
        }
     }' ${name}_QCed_rsID.txt >> ${name}_check_QCed_rsID_file.txt
echo "" >> ${name}_check_QCed_rsID_file.txt


# Check lowest p-value retained after QC
# Check minimum het P-value









