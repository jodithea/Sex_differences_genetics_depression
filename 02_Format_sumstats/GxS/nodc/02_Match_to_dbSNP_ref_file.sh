#!/bin/bash
#PBS -l walltime=01:30:00
#PBS -l mem=150GB
#PBS -l ncpus=1
#PBS -J 0-5

# This script matches each sumstats file to the reference dbSNP file (build 37 and 38 combined with SNPs only MAF > = 0.0001 to make this run in a reasonable timeframe but without dropping SNPs)
# As its SNPs only this will remove any indels etc
# It gives each SNP a consistent markername to use in the meta-analysis and converts between builds
# Here most of the cohorts are build 37, with two build 38 so convert all to build 37

### Preamble ###

directory=/path/01_Format_sumstats/GxS/nodc/ # This is where the plots will be output

dbSNP=/path/Meta-analysis/dbSNP155_GRCh37.p13_vs_dbSNP155_GRCh38.p13/dbSNP_MAF0_0001/NCBI_Hsapiens_dbSNP155_GRCh37_GRCh38.p13_split_multiallelics_MAF0_0001_ALL_matches.txt

files=(/path/01_Format_sumstats/GxS/nodc/*_MDD_GxS_nodc_sumstats.txt)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/01_Format_sumstats/GxS/nodc/||' | sed 's|_build3._MDD_GxS_nodc_sumstats.txt||')

### Run script ###

cd ${directory}

# Check if the input file name contains "build37" or "build38" and run the appropriate block of code
# MarkerName is then created based on the dbSNP ref file so it is consistent across all cohort sumstats so corrcet matching occurs in the meta-analysis
# Note that MarkerName is in the format CHR:BP:Allele:Allele but doesn't necessarily match to A1:A2 (alleles are not flipped - but this is not a problem as METAL deals with allele flipping)

if [[ "${files[${PBS_ARRAY_INDEX}]}" == *"build37"* ]]; then
    echo "Processing file for build 37: ${files[${PBS_ARRAY_INDEX}]}"
    awk 'NR==FNR {
            marker[$5,$6,$2,$3]=$4;
            rsid[$5,$6,$2,$3]=$1;
            next
         }
         {
            if (FNR == 1)
                print $0,"MARKER_build37","rsID_build37";
            else if (($1,$2,$3,$4) in marker)
                print $0,marker[$1,$2,$3,$4],rsid[$1,$2,$3,$4];
            else if (($1,$2,$4,$3) in marker)
                print $0,marker[$1,$2,$4,$3],rsid[$1,$2,$4,$3];
            else
                print $0,"missing","missing";
         }' ${dbSNP} ${files[${PBS_ARRAY_INDEX}]} > ${name}_MDD_GxS_nodc_sumstats_formatted.txt


elif [[ "${files[${PBS_ARRAY_INDEX}]}" == *"build38"* ]]; then
    echo "Processing file for build 38: ${files[${PBS_ARRAY_INDEX}]}"
    awk 'NR==FNR {
            marker[$11,$12,$8,$9]=$4;
            rsid[$11,$12,$8,$9]=$1;
            next
         }
         {
            if (FNR == 1)
                print $0,"MARKER_build37","rsID_build37";
            else if (($1,$2,$3,$4) in marker)
                print $0,marker[$1,$2,$3,$4],rsid[$1,$2,$3,$4];
            else if (($1,$2,$4,$3) in marker)
                print $0,marker[$1,$2,$4,$3],rsid[$1,$2,$4,$3];
            else
                print $0,"missing","missing";
         }' ${dbSNP} ${files[${PBS_ARRAY_INDEX}]} > ${name}_MDD_GxS_nodc_sumstats_formatted.txt
else
    echo "Error: Unrecognized build in file name: ${files[${PBS_ARRAY_INDEX}]}. Skipping."
fi
