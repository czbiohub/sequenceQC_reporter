
#!/bin/bash

DIR="$( cd . "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


#Get Novaseq MACA project names
aws s3 ls s3://czbiohub-maca/remux_data/ > $DIR/projects.txt
cat $DIR/projects.txt | awk -F\| 'gsub(/^ */,"",$0){print $0}' | awk -F ' |/' '{print $2}' > $DIR/projectNames.txt

cat $DIR/projectNames.txt | while read project
do
	echo $project

	mkdir $DIR/00_project_raw_data/$project
	mkdir $DIR/00_project_raw_data/$project/processed_fastq_raw_data
	path=$DIR/00_project_raw_data/$project/processed_fastq_raw_data
	aws s3 ls s3://czbiohub-maca/remux_data/$project/rawdata --recursive > $path/$project'_counts.txt'
	cat $path/$project'_counts.txt' | awk -F ' ' '{print $3, $4}' > $path/$project'_processedCounts.txt'
	echo "Done!"
done