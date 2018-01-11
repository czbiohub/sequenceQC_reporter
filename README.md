## README for Sequence QC Reporter Tool
###### Note: you will need an AWS account with CZ Biohub access. Please see the `eng-support` channel for more information.

> Sequence QC Reporter Tool helps visualize quality control (QC) metrics from sequencing data that is either *unmapped* or *mapped* genome.

* *Unmapped* -- right off the sequencer, Illumina's `bcl2fastq.py` generates a report called "laneBarcode.html".
  * These data give early information about well-to-well clusters passing the `bcl2fastq.py` filter

* *Mapped* -- After your (s)pliced (t)ranscripts (a)lign to your (r)eference, STAR generates a standard report called "log.final.out".
  * These data give information ranging from "Number of input reads" to "Number of splices: Annotated (sjdb)"

* This reporter tool aims to quantify these data and faciliates ease-of-use to determine whether your biological sequence is good to go.

#### Getting started
1. From your terminal
  * Create your project directory by using bash to `aws s3 sync` the required data from the S3 bucket for czbiohub-seqbot to a local directory. This is an example bash script:

```bash
project_names=('171215_M05295_0067_000000000-BHWTV' '171219_M05295_0068_000000000-BHW7B'  '171221_M05295_0070_000000000-BHWGN' '171222_M05295_0071_000000000-BJT3K' '180105_NB501961_0050_AH73JJBGX5' '180105_M05295_0072_000000000-BJR4D')

for project in "${project_names[@]}"; do
  aws s3 sync s3://czbiohub-seqbot/fastqs/$project/results/ ~/data/runs/$project/star_logs --exclude "*" --include "*.log.final.out"
  aws s3 sync s3://czbiohub-seqbot/reports/$project ~/data/runs/$project/reports/
  aws s3 sync s3://czbiohub-seqbot/sample-sheets/ ~/data/runs/$project/sample-sheets/ --exclude "*" --include "$project.csv"
done
```


  * Type `nano sequenceQC_YourReport.sh`
  * Copy contents of the above `bash` chunk and change the `project_names` string to some number correct sequencing run ID (Note: this can be a single run ID or multiple run IDs). 
  * `Ctrl-X` and `Y` to save you script. 
  * To view your script `cat sequenceQC_YourReport.sh`
  * If it looks good: `source sequenceQC_YourReport.sh`...


2. For general use of the current `sequenceQC_reporter` functions
  * Within your R environment, you can import the functions with `source sequenceQC_reporter_functions.R`

3. For detailed information for each function checkout [`sequenceQC_reporter.Rmd`](https://github.com/czbiohub/sequenceQC_reporter/blob/master/sequenceQC_reporter.Rmd) in a working [RStudio](https://www.rstudio.com/) environment.
