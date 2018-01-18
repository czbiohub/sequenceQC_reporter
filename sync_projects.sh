#!/bin/bash

project_names=$@

for project in ${project_names[@]}; do
  aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ ~/data/runs/$project/star_logs --exclude "*" --include "*.log.final.out"
  aws s3 sync s3://czbiohub-seqbot/reports/$project ~/data/runs/$project/reports/
  aws s3 sync s3://czbiohub-seqbot/sample-sheets/ ~/data/runs/$project/sample-sheets/ --exclude "*" --include "$project.csv"
done


 
