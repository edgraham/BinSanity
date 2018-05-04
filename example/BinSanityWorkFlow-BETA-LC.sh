#! /bin/bash
fasta="igm.fa" #Assembly/Fasta File name
coverage_file="/media/acclomator/egraham/applications/BinSanity/example/Infant_gut_assembly.cov.x100.lognorm" #Location of coverage file with full path
location_of_fasta_file="/media/acclomator/egraham/applications/BinSanity/example/" #Location of fasta file give full path leading to directory
threads="50" #number of threads for checkM
clust1="100" 
clust2="50"
#########PASS1##############
Binsanity-lc -f $location_of_fasta_file -l $fasta -c $coverage_file --threads $threads -o PASS1 -m 2000 -v 200  -p -25 -C $clust1 --refine-preference -50 
cd PASS1
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS1-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS1-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna -t $threads . PASS1_checkmLineageWF > PASS1_checkmLineageWF_qa.txt
checkm_analysis -checkM PASS1_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
	mv low_completion/*.fna high_completion
fi
if [ "$(ls -A strain_redundancy)" ]; then
	mv strain_redundancy/*.fna high_redundancy
fi
wait
rm -r low_completion strain_redundancy
###########PASS2##############
cd PASS1
mkdir ../PASS2
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS2/PASS1_HR_LC_SH.all.fa; done
cd ../PASS2
Binsanity-lc -f . -l PASS1_HR_LC_SH.all.fa -c $coverage_file --threads $threads -o . -p -15 -C $clust2 --refine-preference -35
#cd PASS2
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS2-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS2-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna -t $threads . PASS2_checkmLineageWF > PASS2_checkmLineageWF_qa.txt
checkm_analysis -checkM PASS2_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
        mv low_completion/*.fna high_completion
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r low_completion strain_redundancy
mkdir ../PASS3
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS3/PASS2_HR_LC_SH.all.fa; done
###########PASS3############
cd ../PASS3
Binsanity-lc -f . -l PASS2_HR_LC_SH.all.fa -c $coverage_file --threads $threads -o . -p -10 -C $clust2 --refine-preference -25
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS3-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS3-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna -t $threads . PASS3_checkmLineageWF > PASS3_checkmLineageWF_qa.txt
checkm_analysis -checkM PASS3_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
        mv low_completion/*.fna high_completion
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r low_completion strain_redundancy
###############PASS4###############
mkdir ../PASS4
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS4/PASS3_HR_LC_SH.all.fa; done
cd ../PASS4
Binsanity-lc -f . -l PASS3_HR_LC_SH.all.fa -c $coverage_file --threads $threads -o . -p -5 -C $clust2 --refine-preference -25
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS4-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS4-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna -t $threads . PASS4_checkmLineageWF > PASS4_checkmLineageWF_qa.txt
checkm_analysis -checkM PASS4_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
        mv low_completion/*.fna high_redundancy
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r low_completion strain_redundancy
mkdir ../PASS5
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS5/PASS4_HR_LC_SH.all.fa; done
##########PASS5#########
cd ../PASS5
Binsanity-lc -f . -l PASS4_HR_LC_SH.all.fa -c $coverage_file --threads $threads -o . -p -3 -C $clust2 --refine-preference -25
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS5-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS5-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna -t $threads . PASS5_checkmLineageWF > PASS5_checkmLineageWF_qa.txt
checkm_analysis -checkM PASS5_checkmLineageWF_qa.txt
if [ "$(ls -A low_completion)" ];then
        mv low_completion/*.fna high_redundancy
fi
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r low_completion strain_redundancy
mkdir ../PASS6
for f in high_redundancy/*.fna;do cat "$f" >> ../PASS6/PASS5_HR_LC_SH.all.fa; done
################PASS6################
cd ../PASS6
Binsanity-lc -f . -l PASS5_HR_LC_SH.all.fa -c $coverage_file --threads $threads -o . -p -3 -C $clust2 --refine-preference -25
mv BinSanity-Final-bins/*.fna .
rm -r BinSanity-Final-bins
num=1
for file in *.fna; do
       echo -e "$file \t $(printf "PASS6-Bin-%u" $num)" >> file_name_record.txt
       mv "$file" "$(printf "PASS6-Bin-%u" $num).fna"
       let num=$num+1
done
checkm lineage_wf -x fna -t $threads . PASS6_checkmLineageWF > PASS6_checkmLineageWF_qa.txt
checkm_analysis -checkM PASS6_checkmLineageWF_qa.txt
if [ "$(ls -A strain_redundancy)" ]; then
        mv strain_redundancy/*.fna high_redundancy
fi
rm -r strain_redundancy
if [ "$(ls -A high_redundancy)" ]; then
#########PASS7########
        cd high_redundancy
        mkdir ../../PASS7
        mkdir ../../PASS7/BinSanity-LogFiles
        for f in *.fna; do Binsanity -f . -l "$f" --log BinSanity-LogFiles/"$f".binsanity.log -p -2 -c $coverage_file -o ../../PASS7 ;done
        cd ../../PASS7
        mkdir Binsanity-records
        mv * Binsanity-records
        cd Binsanity-records
        checkm lineage_wf -x fna -t $threads . CheckmWF_Initial > CheckmWF_Initial_qa.txt
        checkm_analysis -checkM CheckmWF_Initial_qa.txt
        if [ "$(ls -A strain_redundancy)" ]; then
        	mv strain_redundancy/*.fna high_redundancy
	fi
        rm -r strain_redundancy
        cd high_redundancy
        mkdir ../REFINED
        mkdir ../REFINED/LogFiles
        for f in *.fna; do Binsanity-refine -f . -l "$f" --log ../REFINED/LogFiles/"$f".binsanity-refine.log -p -10 -c $coverage_file -o ../REFINED;done
        cd ../../
        cp Binsanity-records/REFINED/*.fna .
        cp Binsanity-records/high_completion/*.fna .
        cp Binsanity-records/low_completion/*.fna .
        num=1
       	for file in *.fna; do
        	echo -e "$file \t $(printf "PASS7-Bin-%u" $num)" >> file_name_record.txt
        	mv "$file" "$(printf "PASS7-Bin-%u" $num).fna"
                let num=$num+1
        done
        checkm lineage_wf -x fna -t $threads . PASS7_CheckmWF > PASS7_CheckmWF_qa.txt
        checkm_analysis -checkM PASS7_CheckmWF_qa.txt
        #########PASS8#############
        if [ "$(ls -A strain_redundancy)" ]; then
        	mv strain_redundancy/*.fna high_redundancy
	fi
        rm -r strain_redundancy
	if [ "$(ls -A high_redundancy)" ];then
        	cd high_redundancy
        	mkdir ../../PASS8
        	mkdir ../../PASS8/BinSanity-LogFiles
        	for f in *.fna; do Binsanity -f . -l "$f" --log BinSanity-LogFiles/"$f".binsanity.log -p -1 -c $coverage_file -o ../../PASS8 ;done
        	cd ../../PASS8
        	mkdir Binsanity-records
        	mv * Binsanity-records
        	cd Binsanity-records
        	checkm lineage_wf -x fna -t $threads . CheckmWF_Initial > CheckmWF_Initial_qa.txt
        	checkm_analysis -checkM CheckmWF_Initial_qa.txt
        	mv strain_redundancy/*.fna high_redundancy
        	rm -r strain_redundancy
        	cd high_redundancy
        	mkdir ../REFINED
        	mkdir ../REFINED/LogFiles
        	for f in *.fna; do Binsanity-refine -f . -l "$f" --log ../REFINED/LogFiles/"$f".binsanity-refine.log -p -3 -c $coverage_file -o ../REFINED;done
        	cd ../../
        	cp Binsanity-records/REFINED/*.fna .
        	cp Binsanity-records/high_completion/*.fna .
        	cp Binsanity-records/low_completion/*.fna .
		num=1
        	for file in *.fna; do
                	echo -e "$file \t $(printf "PASS8-Bin-%u" $num)" >> file_name_record.txt
                	mv "$file" "$(printf "PASS8-Bin-%u" $num).fna"
                	let num=$num+1
        	done
        	checkm lineage_wf -x fna -t $threads . PASS8_CheckmWF > PASS8_CheckmWF_qa.txt
        	checkm_analysis -checkM PASS8_CheckmWF_qa.txt
       		mkdir ../Final-Genomes
        	cp ../PASS*/high_completion/*.fna ../Final-Genomes
		if [ "$(ls -A ../PASS6/low_completion/)" ]; then
                	cp ../PASS6/low_completion/*.fna ../Final-Genomes
                fi
        	if [ "$(ls -A ../PASS7/low_completion/)" ]; then
                        cp ../PASS7/low_completion/*.fna ../Final-Genomes
                fi
		if [ "$(ls -A ../PASS8/low_completion/)" ]; then
                        cp ../PASS8/low_completion/*.fna ../Final-Genomes
                fi
                if [ "$(ls -A ../PASS8/strain_redundancy/)" ]; then
                        cp ../PASS8/strain_redundancy/*.fna ../Final-Genomes
                fi

        	cd ../
        	mkdir Binsanity-records
        	mv PASS* Binsanity-records
	else
		mkdir ../Final-Genomes
		cp ../PASS*/high_completion/*.fna ../Final-Genomes
                if [ "$(ls -A ../PASS6/low_completion/)" ]; then
                        cp ../PASS6/low_completion/*.fna ../Final-Genomes
                fi
                if [ "$(ls -A ../PASS7/low_completion/)" ]; then
                        cp ../PASS7/low_completion/*.fna ../Final-Genomes
                fi
                if [ "$(ls -A ../PASS7/strain_redundancy/)" ]; then
                        cp ../PASS7/strain_redundancy/*.fna ../Final-Genomes
                fi
		mkdir ../Binsanity-records
		mv ../PASS* ../Binsanity-records
	fi
else
        mkdir ../Final-Genomes
        cp ../PASS*/high_completion/*.fna ../Final-Genomes
        if [ "$(ls -A ../PASS6/low_completion/)" ]; then
        	cp ../PASS6/low_completion/*.fna ../Final-Genomes
        fi
	if [ "$(ls -A ../PASS6/strain_redundancy/)" ]; then
        	cp ../PASS6/low_completion/*.fna ../Final-Genomes
        fi
	cd ../
        mkdir Binsanity-Records
        mv PASS* Binsanity-Records
fi

