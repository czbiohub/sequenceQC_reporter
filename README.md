# sequenceQC_reporter
sequenceQC_reporter reports QC measurements for single-cell RNAseq data and generates a standard report.

###### Note: AWS account with CZ Biohub access required. Please see the `eng-support` channel for more information.
# install
1. At terminal, clone repository `git clone https://github.com/czbiohub/sequenceQC_reporter`
2. Change to working directory `cd sequenceQC_reporter`

# usage
1. Sync by sequencing run ID from `czbiohub-seqbot` S3 bucket: `bash sync_project.sh yourRunID` 
  * Input `yourRunID` as a single project `bash sync_project.sh 171215_M05295_0067_000000000-BHWTV` 
  * Separate with " " for multiple `bash sync_project.sh 171215_M05295_0067_000000000-BHWTV 171221_M05295_0070_000000000-BHWGN`
  * Example project directory
       ```bash
         .
         ├── 00_project_raw_data
         │   ├── 180118_M05295_0075_000000000-D3H9T
         │   │   ├── htseq-counts
         │   │   ├── reports
         │   │   ├── sample-sheets
         │   │   ├── sorted_bams
         │   │   └── star_logs
         │   └── 180126_M05295_0077_000000000-BJNBC
         │       ├── htseq-counts
         │       ├── reports
         │       ├── sample-sheets
         │       ├── sorted_bams
         │       └── star_logs
       ```
2. Load `sequenceQC_reporter` functions in R environment with `source("~/sequenceQC_reporter/sequenceQC_reporter_functions.R")`
