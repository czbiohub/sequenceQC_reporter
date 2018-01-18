# sequenceQC_reporter
> sequenceQC_reporter helps you visualize early quality control (QC) metrics from your sequencing run.
This reporter tool aims to use standard report output from `bcl2fastq.py` and `aws_star.py` as first-pass metrics to determine whether your biological sequence is ready for downstream analysis.

###### Note: you will need an AWS account with CZ Biohub access. Please see the `eng-support` channel for more information.
# install
1. Clone this repository `git clone https://github.com/czbiohub/sequenceQC_reporter/`
2. From your unix terminal, use the `aws s3 sync` commands to sync your data from `czbiohub-seqbot` S3 bucket to local.
* Make `sync_project.sh` executable (only need to do this the first time): `chmod +x sync_project.sh`
* Run with `./sync_project.sh [YOUR_RUN_IDS]`
  * For a single run: `./sync_project.sh 171215_M05295_0067_000000000-BHWTV`
  * For multiple runs: `./sync_project.sh 171215_M05295_0067_000000000-BHWTV 171221_M05295_0070_000000000-BHWGN`
3. Install [RStudio](https://www.rstudio.com/)
  
4. For general use of the current `sequenceQC_reporter` functions
  * Within your R environment, you can import the functions with `source sequenceQC_reporter_functions.R`

5. For detailed information for each function open [`sequenceQC_reporter.Rmd`](https://github.com/czbiohub/sequenceQC_reporter/blob/master/sequenceQC_reporter.Rmd) in your working RStudio environment.
