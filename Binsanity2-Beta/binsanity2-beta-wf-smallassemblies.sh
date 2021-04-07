#! /bin/bash
###This is a placeholder workflow for the beta version of Binsanity2-beta###

fasta=$1 #Assembly/Fasta File name
coverage_file=$3 #Location of coverage file with full path
location_of_fasta_file="$2" #Location of fasta file give full path leading to directory
threads=$4 #number of threads for checkM
#########PASS1###############
Binsanity2-beta -f $location_of_fasta_file -l $fasta -c $coverage_file --checkm_threads $threads  -o PASS1 -m 2000 -v 200  -p -25 --skip-kmeans --refine-preference -50
cd PASS1
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS1-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS1-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna --tab_table -f PASS1_checkmLineageWF_qa.txt -t $threads . PASS1_checkmLineageWF
checkm_analysis -checkM PASS1_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
        for f in low_completion/*.fna; do cat "$f" >> high_redundancy/low_completion.fna; done
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
wait
rm -r low_completion strain_redundancy
###########PASS2##############
mkdir ../PASS2
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS2/PASS1_HR_LC_SH.all.fa; done
cd ../PASS2
Binsanity2-beta -f . -l PASS1_HR_LC_SH.all.fa -c $coverage_file --checkm_threads $threads  -o . -p -15 --skip-kmeans --refine-preference -35
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS2-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS2-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna --tab_table -f PASS2_checkmLineageWF_qa.txt  -t $threads . PASS2_checkmLineageWF
checkm_analysis -checkM PASS2_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
       for f in low_completion/*.fna; do cat "$f" >> high_redundancy/low_completion.fna; done
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r low_completion strain_redundancy
###########PASS3#####################
mkdir ../PASS3
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS3/PASS2_HR_LC_SH.all.fa; done
cd ../PASS3
Binsanity2-beta -f . -l PASS2_HR_LC_SH.all.fa -c $coverage_file --checkm_threads $threads  -o . -p -5 --skip-kmeans --refine-preference -25
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS3-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS3-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna --tab_table -f PASS3_checkmLineageWF_qa.txt -t $threads . PASS3_checkmLineageWF
checkm_analysis -checkM PASS3_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
        for f in low_completion/*.fna; do cat "$f" >> high_redundancy/low_completion.fna; done
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r low_completion strain_redundancy
###############PASS4###############
mkdir ../PASS4
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS4/PASS3_HR_LC_SH.all.fa; done
cd ../PASS4
Binsanity2-beta -f . -l PASS3_HR_LC_SH.all.fa -c $coverage_file --checkm_threads $threads  -o . -p -3 --skip-kmeans --refine-preference -25
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS4-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS4-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna --tab_table -f PASS4_checkmLineageWF_qa.txt -t $threads . PASS4_checkmLineageWF
checkm_analysis -checkM PASS4_checkmLineageWF_qa.txt
################FinalMags##############
mkdir ../Binsanity2-betaWF-Final-Genomes
cd ../
cp PASS*/high_completion/*.fna Binsanity2-betaWF-Final-Genomes
cp PASS4/high_redundancy/*.fna Binsanity2-betaWF-Final-Genomes
cp PASS4/strain_redundancy/*.fna Binsanity2-betaWF-Final-Genomes
cp PASS4/low_completion/*.fna Binsanity2-betaWF-Final-Genomes
mkdir Binsanity-Records
mv PASS* Binsanity-Records
