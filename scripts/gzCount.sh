
#!/bin/bash

project_names=$@

DIR="$( cd . "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# script_base=`echo $(basename $1) | rev | cut -f 2- -d '.' | rev`
# echo $script_base

for project in ${project_names[@]}; do
	echo "Looking for files to download for $project to $DIR"
	mkdir $DIR/00_project_raw_data/$project/processed_fastq_raw_data
	path=$DIR/00_project_raw_data/$project/processed_fastq_raw_data
	#fn=$path/$project'_counts.txt'
	aws s3 ls s3://czbiohub-seqbot/fastqs/$project/ --recursive > $path/$project'_counts.txt'
	fn=$path/$project'_counts.txt'
  	cat $fn | awk -F ' ' '{print $3, $4}' > $path/$project'_processedCounts.txt'
done




aws s3 ls s3://czbiohub-seqbot/fastqs/$project/ --recursive > $path/$project'_counts.txt'