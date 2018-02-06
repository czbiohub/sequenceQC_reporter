
#!/bin/bash

BASHDIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"





project_names=$@

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for project in ${project_names[@]}; do
	echo "Looking for files to download for $project to $DIR"
	mkdir $DIR/00_project_raw_data/$project/processed_fastq_raw_data
	path=$DIR/00_project_raw_data/$project/processed_fastq_raw_data
	fn=$path/$project'_counts.txt'
	aws s3 ls s3://czbiohub-seqbot/fastqs/$project/ --recursive > $path/$project'_counts.txt'
	fn=$path/$project'_counts.txt'
  	cat $fn | awk -F ' ' '{print $3, $4}' > $path/$project'_processedCounts.txt'
done

