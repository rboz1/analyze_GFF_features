#!/bin/bash

# function that checks if the user input chromosome number is valid
check_arg(){
  local script_name="$1"
  local chrom_num="$2"

# check if chrom_num variable has zero length (is empty)
  if [ -z $chrom_num ]
    then
      echo "usage: $script_name <chromosomeID>"
  elif ! [[ $chrom_num =~ ^([1-9]|1[0-9]|2[0-3]|MT|X|Y)$ ]]
    then
      echo "$chrom_num is not a valid chromosomeID ( possible values : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 MT X Y )"
  else
    echo "$chrom_num"
  fi
}

# download and unzip specific chromosome .gff3 file taking in chromosome number as input
download_unzip(){
  file_name='Homo_sapiens.GRCh38.110.chromosome.'
  file_extension='.gff3.gz'
  chrom_num="$1"

# check if file already downloaded
# if downloaded then proceed with GFF analysis
# if not downloaded then download and unzip file first
  if [ -f "chromosome_${chrom_num}.gff3" ]
    then
     :
  else
    curl "ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/${file_name}${chrom_num}${file_extension}" --silent --output chromosome_${chrom_num}${file_extension} && gunzip chromosome_${chrom_num}${file$
  fi
}

# counts all features of specific chromosome taking in the .gff3 file name as input
feature_count(){
  local file_name="$1"
  features=$(awk 'NR>=4 {print $3}' $file_name | sort -f -u)

  for feature in $features
    do
      count=$(awk -v f="$feature" '$3 == f' "$file_name" | wc -l)
      echo "$count $feature"
    done
}

# finds transcripts with the top 10 numbers of exon, CDS, five_prime_UTR, or three_prime_UTR and finds associated genes and gene descriptions
# takes the desired feature as an input
top_ten_count(){
  local file_name="$1"
# desired feature is exon, CDS, five_prime_UTR, or three_prime_UTR
  local desired_feature="$2"

# txt files with subset of .gff3 file
 all_exons=$(awk -v feature="$desired_feature" '$3 ~ feature' "$file_name" > all_exons.txt)
 all_mrna=$(awk '$3 ~ /mRNA/' "$file_name" > all_mrna.txt)
 all_genes=$(awk '$3 ~ /gene/' "$file_name" > all_genes.txt)

# count how many times desired feature (i.e. exon, CDS, etc.) appears for each transcript and save only top 10 transcripts with most desired feature into .txt file
  count=$(awk -v feature="$desired_feature" '{print $9}' all_exons.txt | grep -o -E 'ENST[0-9]{11}' | uniq -c | sort -k1,1 -nr | head -n 10 | sed 's/^//' | awk -v feature="$desired_feature" '{print "Tran$

# using transcript id search for matching gene id and save gene id to .txt file
  > top_10_genes.txt
  for id in $(awk '{print$2}' exon_count.txt)
    do
     awk -v id="$id" '$9 ~ id' all_mrna.txt | grep -o -E 'ENSG[0-9]{11}' | sed 's/^/gene:/' >> top_10_genes.txt
    done

# for each gene id find it's description and save to .txt file
  > descriptions.txt
  for gene in $(cat top_10_genes.txt)
   do
      awk -v gene="$gene" '$9 ~ gene' all_genes.txt | cut -f 9-13 | grep -o 'description=[^;]*' | awk -F'=' '{print $2}' >> descriptions.txt
    done

# print the combined .txt files then remove them
  echo "$(paste exon_count.txt top_10_genes.txt descriptions.txt)"
  rm *.txt
}

# one master function to rule them all (calls the rest of the functions and outputs feature counts and top 10 lists)
gff_analysis(){
  local script_name='./analyze_GFF_features.sh'
  local chrom_num="$1"
  local arg_result=$(check_arg "$script_name" "$chrom_num")

# if user input for chromosome number is valid then proceed with .gff3 analysis
  if [[ $arg_result =~ ^([1-9]|1[0-9]|2[0-3]|MT|X|Y)$ ]]
    then
# download and unzip .gff3 file for the individual chromosome
      download_unzip "$chrom_num"
      echo "                                    "
# count the features and print out
      echo "Feature count chromosome $chrom_num:"
      echo "------------------------------------"
      feature_count "chromosome_${chrom_num}.gff3"
      echo "                                    "
# count the top 10 for each transcript of the chromosome for exons, CDS, five prime UTR, and three prime UTR and print out
      echo "Top 10 count chromosome $chrom_num:"
      echo "------------------------------------"
      echo ">>> transcriptIDs with the highest number of exons"
      top_ten_count "chromosome_${chrom_num}.gff3" 'exon'
      echo "------------------------------------"
      echo ">>> transcriptIDs with the highest number of CDS"
      top_ten_count "chromosome_${chrom_num}.gff3" 'CDS'
      echo "------------------------------------"
      echo ">>> transcriptIDs with the highest number of Five Prime UTR"
      top_ten_count "chromosome_${chrom_num}.gff3" 'five_prime_UTR'
      echo "------------------------------------"
      echo ">>> transcriptIDs with the highest number of Three Prime UTR"
      top_ten_count "chromosome_${chrom_num}.gff3" 'three_prime_UTR'
  fi
}

# call master function and use first user argument as the chromosome number input
gff_analysis "$1"
