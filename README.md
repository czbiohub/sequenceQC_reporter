# sequenceQC_reporter
> sequenceQC_reporter helps you visualize early quality control (QC) metrics from your sequencing run.
This reporter tool aims to use standard report output from `bcl2fastq.py` and `aws_star.py` to determine whether your sequence sample is ready for downstream analysis.

###### Note: you will need an AWS account with CZ Biohub access. Please see the `eng-support` channel for more information.
# install
1. At your terminal, download this repository `git clone https://github.com/czbiohub/sequenceQC_reporter`
2. Change to working directory `cd sequenceQC_reporter`
3. From command-line, sync your project data from `czbiohub-seqbot` S3 bucket to a local directory with `sync_project.sh`
  * Make this script executable with `chmod +x sync_project.sh`
4. Use current functions in `sequenceQC_reporter` in [RStudio](https://www.rstudio.com/)
  * In an `.Rmd` or `.R` file import functions with `source(sequenceQC_reporter_functions.R)`

# usage
1. Setup project directory
*  Example sync `./sync_project.sh [YOUR_RUN_IDS]`
  * Single run: `./sync_project.sh 171215_M05295_0067_000000000-BHWTV`
  * Multiple runs: `./sync_project.sh 171215_M05295_0067_000000000-BHWTV 171221_M05295_0070_000000000-BHWGN`

2. For walkthrough, view [`walkthrough_sequenceQC_reporter.Rmd`](https://github.com/czbiohub/sequenceQC_reporter/blob/master/walkthrough_sequenceQC_reporter.Rmd)
3. For detailed usage of each function, view [`sequenceQC_reporter.Rmd`](https://github.com/czbiohub/sequenceQC_reporter/blob/master/sequenceQC_reporter.Rmd)
