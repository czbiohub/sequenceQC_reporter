## README for Sequence QC Reporter Tool

###### Note: you will need an AWS account with CZ Biohub access. Please see the `eng-support` channel for more information.

Create your project directory by using bash to `aws s3 sync` the required data from the S3 bucket for czbiohub-seqbot to a local directory. This is an example bash script:

```bash
project_names=('171215_M05295_0067_000000000-BHWTV' '171219_M05295_0068_000000000-BHW7B'  '171221_M05295_0070_000000000-BHWGN' '171222_M05295_0071_000000000-BJT3K' '180105_NB501961_0050_AH73JJBGX5' '180105_M05295_0072_000000000-BJR4D')

for project in "${project_names[@]}"; do
  aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ ~/data/runs/$project/star_logs --exclude "*" --include "*.log.final.out"
  aws s3 sync s3://czbiohub-seqbot/reports/$project ~/data/runs/$project/reports/
  aws s3 sync s3://czbiohub-seqbot/sample-sheets/ ~/data/runs/$project/sample-sheets/ --exclude "*" --include "$project.csv"
done
```


Using bash in your Terminal, first run `nano sequenceQC_YourReport.sh`. Then, copy the contents of the above chunk and change `project_names` to the correct sequencing run ID (you can add a single or multiple run ID). `Ctrl-X` and `Y` to save the bash script. To view use `cat sequenceQC_YourReport.sh`. If this looks good, then `source sequenceQC_YourReport.sh`...


For general use of the `sequenceQC_reporter` functions, use `source sequenceQC_reporter_functions.R` and execute in your R environment.

For detailed information for each function open the R Markdown file `sequenceQC_reporter.Rmd` in a working RStudio environment.
