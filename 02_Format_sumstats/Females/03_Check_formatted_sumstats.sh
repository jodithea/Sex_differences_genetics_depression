#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=20GB
#PBS -l ncpus=1
#PBS -J 0-5

# After matching each sumstats file to the reference dbSNP file (build 37 and 38 combined with SNPs only MAF > = 0.0001) to get a consistent marker name
# Check these files


### Preamble ###

directory=/path/01_Format_sumstats/Females/ # This is where the plots will be output

files=(/path/01_Format_sumstats/Females/*_MDD_female_sumstats_formatted.txt)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/01_Format_sumstats/Females/||' | sed 's|_MDD_female_sumstats_formatted.txt||')


### Run script ###

cd ${directory}

echo "Checking file: ${files[${PBS_ARRAY_INDEX}]}" > ${name}_MDD_female_sumstats_check_formatted_file.txt


## Check how many SNPs did not match to marker name ##

echo "Number rows missing (including non-SNPs):" >> ${name}_MDD_female_sumstats_check_formatted_file.txt
grep -c "missing" ${files[${PBS_ARRAY_INDEX}]} >> ${name}_MDD_female_sumstats_check_formatted_file.txt
echo "" >> ${name}_MDD_female_sumstats_check_formatted_file.txt

echo "Number rows missing (SNPs only):" >> ${name}_MDD_female_sumstats_check_formatted_file.txt 
awk '$13 == "missing" && length($3) == 1 && length($4) == 1' ${files[${PBS_ARRAY_INDEX}]} | wc -l >> ${name}_MDD_female_sumstats_check_formatted_file.txt
echo "" >> ${name}_MDD_female_sumstats_check_formatted_file.txt

echo "Number rows with \"no_match\" (SNP present in the reference dbSNP file but not in the particular build using):" >> ${name}_MDD_female_sumstats_check_formatted_file.txt
grep -c "no_match" ${files[${PBS_ARRAY_INDEX}]} >> ${name}_MDD_female_sumstats_check_formatted_file.txt
echo "" >> ${name}_MDD_female_sumstats_check_formatted_file.txt

echo "Number rows with \"no_rsid\" (SNP present in one build, mapped to the other build and coords found but no matching rsid found):" >> ${name}_MDD_female_sumstats_check_formatted_file.txt
grep -c "no_rsid" ${files[${PBS_ARRAY_INDEX}]} >> ${name}_MDD_female_sumstats_check_formatted_file.txt
echo "" >> ${name}_MDD_female_sumstats_check_formatted_file.txt



## Remove lines with no match to dbSNP reference file ##

awk '($13 != "missing" && !($13 ~ /^no_match/) && $14 != "no_rsid") {print}' ${files[${PBS_ARRAY_INDEX}]} > ${name}_MDD_female_sumstats_formatted_formetaanalysis.txt



## Check min and max of each column to see if anything weird ##

# Before removing missing lines

echo "Min and max of each column before removing missing lines" >> ${name}_MDD_female_sumstats_check_formatted_file.txt
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
     }' ${files[${PBS_ARRAY_INDEX}]} >> ${name}_MDD_female_sumstats_check_formatted_file.txt
echo "" >> ${name}_MDD_female_sumstats_check_formatted_file.txt


# After removing missing lines
echo "Min and max of each column after removing missing lines" >> ${name}_MDD_female_sumstats_check_formatted_file.txt
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
     }' ${name}_MDD_female_sumstats_formatted_formetaanalysis.txt >> ${name}_MDD_female_sumstats_check_formatted_file.txt
echo "" >> ${name}_MDD_female_sumstats_check_formatted_file.txt




# Check lowest p-value retained after matching consistent marker
# Check QC on MAF and info done correctly


