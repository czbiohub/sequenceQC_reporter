# sequenceQC_reporter
> sequenceQC_reporter enables user to apply multiple quality control (QC) metrics to single-cell RNA sequences.
Aim for this tool is to generate a standard report as well as provide the necessary data-types for the user to determine if their is sufficient for downstream analysis.

###### Note: AWS account with CZ Biohub access required. Please see the `eng-support` channel for more information.
# install
1. At terminal, clone repository `git clone https://github.com/czbiohub/sequenceQC_reporter`
2. Change to working directory `cd sequenceQC_reporter`

# usage
1. Load project data from `czbiohub-seqbot` S3 bucket to a local directory. See usage for `sync_project.sh`:
   * Option A. Make executable with `chmod +x` > `./sync_project.sh [YOUR_RUN_IDS]`
   * Option B. Run as bash script > `bash sync_project.sh [YOUR_RUN_IDS]`
         * Single run > ` bash sync_project.sh 171215_M05295_0067_000000000-BHWTV`
         * Multiple runs > `bash sync_project.sh 171215_M05295_0067_000000000-BHWTV 171221_M05295_0070_000000000-BHWGN`
   * Example project structure
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
2. Load `sequenceQC_reporter` function in R environment with `source("~/sequenceQC_reporter/sequenceQC_reporter_functions.R")`
3. For walkthrough, please see [`walkthrough_sequenceQC_reporter.Rmd`](https://github.com/czbiohub/sequenceQC_reporter/blob/master/walkthrough_sequenceQC_reporter.Rmd)
4. For furhter usage, please see [`sequenceQC_reporter.Rmd`](https://github.com/czbiohub/sequenceQC_reporter/blob/master/sequenceQC_reporter.Rmd)
