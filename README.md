# sequenceQC_reporter
>sequenceQC_reporter reports QC measurements for single-cell RNAseq data and generates a standard report.

## Usage
Use with your sequencing run prefix (should have format that looks like this `YYMMDD_A00111`). AWS account with CZ Biohub access required. 

```bash
git clone https://github.com/czbiohub/sequenceQC_reporter
cd sequenceQC_reporter/scripts/
bash scripts/sync_projects.sh YYMMDD_A00111
```

Example project directory
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

Open R file or [RStudio](https://www.rstudio.com/) notebook after sync completes.
```r
source('~/sequenceQC_reporter/sequenceQC_reporter_functions.R')
yourRunID = '180202_NB501961_0059_AHFLYGBGX5'
projectDir = paste0("~/sequenceQC_reporter/00_project_raw_data/", yourRunID)
samplesheet = loadSamplesheet(projectDir)
plotReads = mapSheet(projectDir, samplesheet, gzCount(projectDir), "Processed fastq.gz counts")
```

