# sequenceQC_reporter
>sequenceQC_reporter reports QC measurements for single-cell RNAseq data and generates a standard report.

## Usage
Use with your sequencing run prefix (should have format that looks like this `YYMMDD_A00111`). AWS account with CZ Biohub access required. 

```bash
git clone https://github.com/czbiohub/sequenceQC_reporter
cd sequenceQC_reporter
bash scripts/sync_projects.sh YYMMDD_A00111
```

The script will sync the files needed to run the R script.

Example project directory
```bash
  .
  ├── 00_project_raw_data
  │   ├── YYMMDD_A00111
  │   │   ├── htseq-counts
  │   │   ├── reports
  │   │   ├── sample-sheets
  │   │   ├── sorted_bams
  │   │   └── star_logs
```



Open R file or [RStudio](https://www.rstudio.com/) notebook after sync completes. On the RStudio navigation pane: File > 'Create New project' > 'Existing directory' > Browse and set path to `/sequenceQC_reporter/` > Create project.

Navigate to the parent project directory (this was just created) > Open `platemap_tutorial.Rmd` in RStudio and follow instructions there.

```r
source('~/sequenceQC_reporter/sequenceQC_reporter_functions.R')
yourRunID = 'YYMMDD_A00111'
projectDir = paste0("~/sequenceQC_reporter/00_project_raw_data/", yourRunID)
samplesheet = loadSamplesheet(projectDir)
key_parameters = c('Uniquely mapped reads number', 'Number of input reads')
starlog = loadStarLog(projectDir, key_parameters)
df = sortSheetData(projectDir, samplesheet, whatever = starlog, saveFile = FALSE)
heatmap = map_whatever(projectDir, data = df, column_name = 'clusters', log2 = TRUE, savePlot = FALSE)
```

