#!/bin/bash

project_names=$@

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


for project in ${project_names[@]}; do
echo "Looking for data to download for $project to $DIR"

	# if [ ! -f $DIR/00_project_raw_data/$project/star_logs/*.log.final.out ] ; then
	# 	echo "Loading log.final.out for $project to $DIR"
	# 	aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ $DIR/00_project_raw_data/$project/star_logs --exclude "*" --include "*.log.final.out"
	# fi

	aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ $DIR/00_project_raw_data/$project/star_logs --exclude "*" --include "*.log.final.out"
	aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ $DIR/00_project_raw_data/$project/htseq-counts --exclude "*" --include "*htseq-count.txt"
	aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ $DIR/00_project_raw_data/$project/sorted_bams/ --exclude "*" --include "*.sorted.bam"
	aws s3 sync s3://czbiohub-seqbot/reports/$project $DIR/00_project_raw_data/$project/reports/ #--exclude "*" --include "laneBarcode.html"
	aws s3 sync s3://czbiohub-seqbot/reports/$project $DIR/00_project_raw_data/$project/reports/ --exclude "*" --include "laneBarcode.html"
	aws s3 sync s3://czbiohub-seqbot/sample-sheets/ $DIR/00_project_raw_data/$project/sample-sheets/ --exclude "*" --include "$project.csv"

done
