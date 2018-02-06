#!/usr/bin/bash

project_names=$@

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for project in ${project_names[@]}; do

  #mkdir $DIR/00_project_raw_data/$project/sorted_bams/filtering/
	#run inner_distance.py 
  	for filename in `ls $DIR/00_project_raw_data/$project/sorted_bams/ | grep sorted.bam$`; do
  		echo "$filename exists"
  		(ls $DIR/00_project_raw_data/$project/sorted_bams/$filename.inner_distance.txt && 
  			echo "Processing $filename") || 
  			inner_distance.py -i $DIR/00_project_raw_data/$project/sorted_bams/$filename -o $DIR/00_project_raw_data/$project/sorted_bams/$filename -r $DIR/reference/mm10_RefSeq.bed
	   done



	#filter inserts	
	#mkdir $DIR/00_project_raw_data/$project/sorted_bams/filtering/
	
	for filename in `ls $DIR/00_project_raw_data/$project/sorted_bams/ | grep inner_distance.txt$`; do
		  cat $DIR/00_project_raw_data/$project/sorted_bams/$filename > $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename'_unfiltered.txt';
  		cat $DIR/00_project_raw_data/$project/sorted_bams/$filename | grep '\tsameTranscript=Yes,sameExon=No,dist=mRNA' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename'_STY_SEN_dMRNA.txt';
  		cat $DIR/00_project_raw_data/$project/sorted_bams/$filename | grep '\tsameTranscript=Yes,nonExonic=Yes,dist=genomic' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename'_STY_NEY_dG.txt';
  		cat $DIR/00_project_raw_data/$project/sorted_bams/$filename | grep '\tsameTranscript=Yes,sameExon=Yes,dist=mRNA' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename'_STY_SEY_dMRNA.txt';
  	done
  	echo "Done filtering 'insert_size' files" 

done


for project in ${project_names[@]}; do
	
	#strip strings from filtered files
	mkdir $DIR/00_project_raw_data/$project/sorted_bams/filtering/stripped/

    for filename in `ls $DIR/00_project_raw_data/$project/sorted_bams/filtering/ | grep unfiltered.txt$`; do
  		#echo $filename
  		cat $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename | awk -F '\t' '{print $2}' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/stripped/$filename;
  	done
    echo "Done filtering 'unfiltered.txt' files"

  	for filename in `ls $DIR/00_project_raw_data/$project/sorted_bams/filtering/ | grep STY_SEN_dMRNA.txt$`; do
  		#echo $filename
  		cat $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename | awk -F '\t' '{print $2}' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/stripped/$filename;
  	done
  	echo "Done filtering 'STY_SEN_dMRNA.txt' files"

  	for filename in `ls $DIR/00_project_raw_data/$project/sorted_bams/filtering/ | grep STY_NEY_dG.txt$`; do
  		#echo $filename
  		cat $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename | awk -F '\t' '{print $2}' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/stripped/$filename;
  	done
  	echo "Done filtering 'STY_NEY_dG.txt' files"
	
	for filename in `ls $DIR/00_project_raw_data/$project/sorted_bams/filtering/ | grep STY_SEY_dMRNA.txt$`; do
  		#echo $filename
  		cat $DIR/00_project_raw_data/$project/sorted_bams/filtering/$filename | awk -F '\t' '{print $2}' > $DIR/00_project_raw_data/$project/sorted_bams/filtering/stripped/$filename;
  	done
  	echo "Done filtering 'STY_SEY_dMRNA.txt' files"

done