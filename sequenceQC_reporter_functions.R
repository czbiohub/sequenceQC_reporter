#Generate alphanumeric index of 384-well plate format. This function can be used alone or called within other functions. Input is `NUMBER_OF_PLATES`. 
generate_wells = function(NUMBER_OF_PLATES){
  w_i = {}
  well_index = {}
  for(i in 1:NUMBER_OF_PLATES){
    for (j in LETTERS[1:16]){
      w_i[[j]] = paste0(j, rep(1:24, each=1))
    }
    well_index[[i]] = unlist(w_i)
  }
  return(unlist(well_index))
}

#Sort data to samplesheet index to "Sample_ID"
sortSheetData = function(samplesheet, starLog, laneBarcode){
  df = left_join(samplesheet, starLog, by = "Sample_ID", sort = FALSE)
  df = left_join(df, laneBarcode, by = "Sample_ID", sort = FALSE)
  return(df)
}

#Load samplesheet for data generation and plotting. This function can be called separately for csv-only (data) and png-only (plot) generation. 
loadSamplesheet = function(projectDir){
  
  runID = list(strsplit(projectDir,"00_project_raw_data/")[[1]][2])
  for(i in 1:length(projectDir)){
    
    temp = list.files(projectDir, pattern = paste0(runID, ".csv"), full.names = TRUE, recursive = TRUE)
    samplesheet = read_csv(temp, col_names = TRUE, col_types = cols(), skip = 1)
    
    NUMBER_OF_PLATES = nrow(samplesheet)/384
    samplesheet$well_index = generate_wells(NUMBER_OF_PLATES)[1:nrow(samplesheet)]
    samplesheet$row = as.factor(samplesheet$well_index %>% str_extract(pattern = "^[A-Z]{1}"))
    samplesheet$col = as.factor(samplesheet$well_index %>% str_extract(pattern = "[0-9]{1,2}"))
    samplesheet$Description = rep(1:NUMBER_OF_PLATES, each=384)
    
    if(is.integer(samplesheet$`Sample_ID`) == TRUE){
      samplesheet$`Sample_ID` = as.character(samplesheet$`Sample_ID`)
    }
    print(paste0("Samplesheet loaded for sequencing run: ", runID))
  }
  return(samplesheet)
}

#Generate your FASTQ stats - this function scrapes Illumina bcl2fastq html report files and sorts data like "passed-filter" clusters, samples, and barcodes.
loadLaneBarcode = function(projectDir){
  
  runID = {} ; for(i in projectDir){runID[[i]]= strsplit(i,"00_project_raw_data/")[[1]][2]}
  laneBarcodeList = {}
  for(i in 1:length(projectDir)){
    barcodeFile = list.files(projectDir, pattern = "laneBarcode.html", full.names = TRUE, recursive = TRUE)
    for(j in barcodeFile){

      sampleList = {} ; barcodeList = {} ; clusterList = {}
      temp = read_lines(j, skip = 42)
      sampleList[[j]] = temp[seq(2, length(temp), 14)]
      barcodeList[[j]] = temp[seq(3, length(temp), 14)]
      clusterList[[j]] = temp[seq(4,length(temp),14)]
    }
    
    samples = {} ; clusters = {} ; barcodes = {}
    for(j in barcodeFile){
      samples[[j]] = unlist(sampleList[[j]])
      clusters[[j]] = unlist(clusterList[[j]])
      barcodes[[j]] = unlist(barcodeList[[j]])
    }
    
    laneBarcodeReport = list(unlist(samples), unlist(clusters), unlist(barcodes))
    
    for (x in 1:length(laneBarcodeReport)){laneBarcodeReport[[x]] = str_extract(laneBarcodeReport[[x]], "(?<=>)[^<]+")}
    
    laneBarcodeReport[[2]] = as.numeric(gsub(",", "", laneBarcodeReport[[2]]))
    df = data_frame(Sample_ID = laneBarcodeReport[[1]], clusters = laneBarcodeReport[[2]], barcodes = laneBarcodeReport[[3]])
    print(paste0("Demultiplexation report loaded for sequencing run: ", runID))
    
  }
  return(df)
}

#Return key-value pair for yourParameter
return.parameterKey = function(projectDir, yourParameter, filePattern){
  
  filenames = list.files(projectDir, pattern = filePattern, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  temp = read_lines(filenames[1])
  allParameters = {}
  
  for(i in 1:length(temp)){
    
    allParameters[[i]] = strsplit(temp[i], " \\|")[[1]][1]
    print("Got keys")
    allParameters = gsub("^\\s+","", allParameters)
    print("Loading parameters")
  }
  
  names(allParameters) = allParameters
  keys = {}
  keys$index = match(yourParameter, allParameters)
  print("Match parameters")
  keys$keynames = allParameters[yourParameter]
  
  checklist = {}; for(j in yourParameter){checklist[[j]] = j %in% keys$keynames}
  
  if(any(checklist == FALSE)){
    print(paste0("Your parameter ", j, " wasn't found. Stopping program."))
  } else {
    for(j in yourParameter){
      print(paste0("Loading your parameter: ", j))
    }
    return(keys)
  }
}

#Generate your STAR stats - This input for his function is any single line or multiple lines as a PARAMETER (i.e. "Mapped Reads", "Percent Mapped Reads", "Insert Average Length") from the "log.final.out" file. This file a standard log generated afer successful STAR alignment.

loadStarLog = function(projectDir, yourParameter){
  
  filenames = list.files(projectDir, pattern = "*log.final.out", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  runID = list(strsplit(projectDir,"00_project_raw_data/")[[1]][2])
  for(i in 1:length(filenames)){runID$sampleID[[i]] = strsplit(filenames[i], "logs/|_S")[[1]][2]}
  
  filePattern = "*log.final.out"
  runID$parameterKey = return.parameterKey(projectDir, yourParameter, filePattern)
  
  for(j in 1:length(runID$parameterKey$index)){
    currentKey = runID$parameterKey$index[j]
    fileData = {}
    for(i in 1:length(filenames)){
      if(is.na(i)){
        print("No log files for this run. Skipping.")
        break
      } else {
        temp = read_lines(filenames[i])
        fileData[[i]] = strsplit(temp[currentKey], "\t")[[1]][2]
        fileData = gsub("%","", fileData)
      }
    }
    runID$parameterData[[j]] = fileData
  }
  
  log.star = {} ; for(i in 1:length(runID$parameterData)){log.star[[i]] = as.numeric(unlist(runID$parameterData[i]))}
  names(log.star) = unlist(runID$parameterKey$keynames)
  log.star$Sample_ID= unlist(runID$sampleID)
  df = as_tibble(log.star[1:length(log.star)])
  return(df)
  
}

#Generate figures from Lane Barcode data
heatmap.laneBarcode = function(projectDir, samplesheet, laneBarcode){
  runID = strsplit(projectDir,"00_project_raw_data/")[[1]][2]
  laneBarcodeReport = left_join(samplesheet, laneBarcode, by = "Sample_ID", sort = FALSE)
  write_csv(x = laneBarcodeReport, path = paste0(projectDir,"/", runID,"_laneBarcodeReport.csv"))
  
  print(paste0("Wrote your Lane Barcode Report for ", runID, "to ", projectDir,"/"))
  
  data = laneBarcodeReport
  data$clusters = as.numeric(data$clusters)
  
  plotTitle = paste0("Lane Barcode Heatmap: ", runID, " Parameter: Unmapped clusters")
  
  log2_heatmap = ggplot() +
    geom_raster(data, 
                mapping = aes(x = col, y = row, fill = log2(data$clusters))) +
    scale_x_discrete(limits = c(1:24)) +
    scale_y_discrete(limits = rev(levels(data$row))) +
    coord_equal() +
    scale_fill_viridis(guide_legend(title = "log2")) +
    theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    facet_wrap(~data$Description) +
    ggtitle(plotTitle)
  ggsave(paste0(projectDir,"/", runID, "_laneBarcode_log2scale_heatmap.png"))
  print("Saved heatmap for log2-scaled clusters")
  
  unscaled_heatmap = ggplot() +
    geom_raster(data, 
                mapping = aes(x = col, y = row, fill = data$clusters)) +
    scale_x_discrete(limits = c(1:24)) +
    scale_y_discrete(limits = rev(levels(data$row))) +
    coord_equal() +
    scale_fill_viridis(guide_legend(title = "unscaled")) +
    theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    facet_wrap(~data$Description) +
    ggtitle(plotTitle)
  ggsave(paste0(projectDir,"/", runID, "_laneBarcode_unscaled_heatmap.png"))
  print("Saved heatmap for unscaled clusters")
  
  heatmaps = list(log2_heatmap, unscaled_heatmap)
  
  return(heatmaps)
  
}

#Generate figures from STAR log data.

heatmap.starLog = function(projectDir, yourParameter, samplesheet, starLog){
  
  runID = strsplit(projectDir,"00_project_raw_data/")[[1]][2]
  starLogReport = left_join(samplesheet, starLog, by = "Sample_ID", sort = FALSE)
  write_csv(x = starLogReport, path = paste0(projectDir,"/", runID,"_starLogReport.csv"))
  print(paste0("Wrote your Star Log Report for ", runID, "to ", projectDir,"/"))
  data = starLogReport
  figures = {}
  
  for(i in 1:length(yourParameter)){
    
    plotTitle = paste0("STAR alignment parameter: ", yourParameter[i])
    
    if(str_detect(yourParameter[i], "%")){
      parameterName = gsub("%", "percent", yourParameter[i])
      parameterName = gsub(" ", "", parameterName)
    } else {parameterName = gsub(" ", "", yourParameter[i])}
    
    print(paste0("Generating platemaps for ", parameterName))
    
    log2_heatmap = ggplot() +
      geom_raster(data,
                  mapping = aes(x = col, y = row, fill = log2(data[[yourParameter[i]]]))) +
      scale_x_discrete(limits = c(1:24)) + scale_y_discrete(limits = rev(levels(data$row))) + coord_equal() +
      scale_fill_viridis(guide_legend(title = paste0("log2"))) +
      theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
            legend.title = element_text(size = 7), legend.text = element_text(size = 7),
            axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
      facet_wrap(~data$Description) +
      ggtitle(plotTitle)
    ggsave(paste0(projectDir,"/", runID, "_", parameterName,"_starLogReport_log2scale_heatmap.png"))
    
    unscaled_heatmap = ggplot() +
      geom_raster(data,
                  mapping = aes(x = col, y = row, fill = data[[yourParameter[i]]])) +
      scale_x_discrete(limits = c(1:24)) + scale_y_discrete(limits = rev(levels(data$row))) + coord_equal() +
      scale_fill_viridis(guide_legend(title = paste0("unscaled"))) +
      theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
            legend.title = element_text(size = 7), legend.text = element_text(size = 7),
            axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
      facet_wrap(~data$Description) +
      ggtitle(plotTitle)
    ggsave(paste0(projectDir,"/", runID, "_", parameterName,"_starLogReport_unscaled_heatmap.png"))
    figures$log2[[i]] = log2_heatmap
    figures$unscaled[[i]] = unscaled_heatmap
  }
  return()
}


#Load counts by geneName (important: geneName must be a character string that exists in the htseq count file)
countsByGene = function(projectDir, geneName = "Rn45s"){
  
  
  filenames = list.files(projectDir, pattern = "*htseq-count", full.names = TRUE, recursive = TRUE)
  indexFile = read_lines(filenames[1])
  temp = {} ; for(i in 1:length(indexFile)){
    temp[[i]] = strsplit(indexFile[i], "\t")[[1]][1]
  }
  
  geneKey = temp
  #geneName = "Rn45s"
  index = as.numeric(grep(geneName, geneKey))
  
  cellNames = {} ; totalCounts = {}
  for(i in 1:length(filenames)){
    temp = read_lines(filenames[i])
    
    countValue = as.numeric(strsplit(temp[index], "\t")[[1]][2])
    if(countValue != 0){
      #print("Nonzero")
      targetGene = strsplit(temp[index], "\t")[[1]][1]
      if(targetGene == geneName){
        totalCounts[[i]] = countValue
        cellNames[[i]] = strsplit(filenames[i], "*counts/|_S")[[1]][2]
      }
    }
  }
  geneCounts = list(totalCounts, cellNames)
  return(geneCounts)
  
}







