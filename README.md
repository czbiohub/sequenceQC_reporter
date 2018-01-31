# sequenceQC_reporter
sequenceQC_reporter reports QC measurements for single-cell RNAseq data and generates a standard report.

###### Note: AWS account with CZ Biohub access required. Please see the `eng-support` channel for more information.
# install
1. At terminal, clone repository `git clone https://github.com/czbiohub/sequenceQC_reporter`
2. Change to working directory `cd sequenceQC_reporter`

# usage
1. Load project data from `czbiohub-seqbot` S3 bucket to a local directory. 
2. To sync: `bash sync_project.sh yourRunID` 
  * Input `yourRunID` as a single project `bash sync_project.sh 171215_M05295_0067_000000000-BHWTV` 
  * Separate with " " for multiple projects `bash sync_project.sh 171215_M05295_0067_000000000-BHWTV 171221_M05295_0070_000000000-BHWGN`
  * Example project directory
       ```bash
       ├── 171215_M05295_0067_000000000-BHWTV
       │   ├── reports
       │   ├── sample-sheets
       │   └── star_logs
       ├── 171221_M05295_0070_000000000-BHWGN
       │   ├── reports
       │   ├── sample-sheets
       │   └── star_logs
       ```
2. Load `sequenceQC_reporter` functions in R environment with `source("~/sequenceQC_reporter/sequenceQC_reporter_functions.R")`
