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
sortSheetData = function(projectDir, samplesheet, whatever, saveFile = FALSE){
  runID = list(strsplit(projectDir,"00_project_raw_data/")[[1]][2])
  df = left_join(samplesheet, whatever, by = "Sample_ID", sort = FALSE)
  if(saveFile == TRUE){
    dir.create(paste0(projectDir, "/datasheets/"), showWarnings = FALSE, recursive = TRUE)
    filePath = paste0(projectDir, "/datasheets/")
    write_csv(df, path = paste0(filePath, runID,"_datasheet.csv"))}
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
    
    if(is.integer(samplesheet$`Sample_ID`) == TRUE){samplesheet$`Sample_ID` = as.character(samplesheet$`Sample_ID`)}
    
    if(samplesheet$Sample_ID[1] != samplesheet$Sample_Name[1]){samplesheet$Sample_ID[1] = samplesheet$Sample_Name[1]}
    
    print(paste0("Samplesheet loaded for sequencing run: ", runID))
  }
  
  return(samplesheet)
}

#Parse demultiplex html report files and sort data. Load demultiplex html report data
parseLaneBarcode = function(barcodeFile){
  
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
  df$path = barcodeFile
  return(df)
  
}
loadLaneBarcode = function(projectDir){
  
  runID = {} ; for(i in projectDir){runID[[i]]= strsplit(i,"00_project_raw_data/")[[1]][2]}
  laneBarcodeList = {}
  for(i in 1:length(projectDir)){
    barcodeFile = list.files(projectDir, pattern = "laneBarcode.html", full.names = TRUE, recursive = TRUE)
    
    #If no barcode files exist
    if(length(barcodeFile) == 0){
      print(paste0("No `laneBarcode.html` was found for sequencing run: ", runID))
      break
      
    }
    
    #For batch directories of barcode files
    if(length(barcodeFile) >= 1){
      allLanes = data_frame()
      for(j in 1:length(barcodeFile)){
        fn = parseLaneBarcode(barcodeFile[j])
        allLanes = rbind(allLanes, fn)
        
      }
      print(paste0("Binding all reports to a data frame for ", runID))
    }
    
    print(paste0("Demultiplexation report loaded for sequencing run: ", runID))
    
    allLanes = allLanes %>% group_by(Sample_ID) %>% distinct(x, .keep_all = TRUE) %>% ungroup()
    allLanes$clusters = as.numeric(allLanes$clusters)
    return(allLanes)
  }
}






#Generate your STAR stats - This input for his function is any single line or multiple lines as a PARAMETER (i.e. "Mapped Reads", "Percent Mapped Reads", "Insert Average Length") from the "log.final.out" file. This file a standard log generated afer successful STAR alignment.
#Return key-value pair for parameters
return.parameterKey = function(projectDir, key_parameter, filePattern){
  
  filenames = list.files(projectDir, pattern = filePattern, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  temp = read_lines(filenames[1])
  allParameters = {}
  
  for(i in 1:length(temp)){
    allParameters[[i]] = strsplit(temp[i], " \\|")[[1]][1]
    allParameters = gsub("^\\s+","", allParameters)
  }
  
  names(allParameters) = allParameters
  keys = {}
  keys$index = match(key_parameter, allParameters)
  print("Match parameters")
  keys$keynames = allParameters[key_parameter]
  
  checklist = {}; for(j in key_parameter){checklist[[j]] = j %in% keys$keynames}
  
  if(any(checklist == FALSE)){print(paste0("Your parameter ", j, " wasn't found. Stopping program."))
  } else {return(keys)}
  
}

loadStarLog = function(projectDir, key_parameter){
  
  filenames = list.files(projectDir, pattern = "*log.final.out", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  if(length(filenames) == 0){print(paste0("No `log.final.out` files were found for sequencing run: ", runID)) ; break}
  
  runID = list(strsplit(projectDir,"00_project_raw_data/")[[1]][2])
  for(i in 1:length(filenames)){runID$sampleID[[i]] = strsplit(filenames[i], "logs/|_S")[[1]][2]}
  
  filePattern = "*log.final.out"
  runID$parameterKey = return.parameterKey(projectDir, key_parameter, filePattern)
  
  for(j in 1:length(runID$parameterKey$index)){
    currentKey = runID$parameterKey$index[j]
    fileData = {}
    for(i in 1:length(filenames)){
      if(is.na(i)){print("No log files for this run. Skipping.") ; break
      } else {
        temp = read_lines(filenames[i])
        fileData[[i]] = strsplit(temp[currentKey], "\t")[[1]][2]
        fileData = gsub("%","", fileData)
      }
    }
    runID$parameterData[[j]] = fileData
  }
  
  df = {} ; for(i in 1:length(runID$parameterData)){df[[i]] = as.numeric(unlist(runID$parameterData[i]))}
  names(df) = unlist(runID$parameterKey$keynames)
  df$Sample_ID= unlist(runID$sampleID)
  df = as_tibble(df[1:length(df)])
  df[1] = round(df[1]/4, digits = 0)
  return(df)
}









#Parse htseq-count files to count by a given input for geneName (important: geneName must be a character string that exists in the htseq count file). Load a dataframe with counts for a list containing multiple geneNames indexed by 'Sample_ID'
countsByGene = function(projectDir, geneName){
  
  filenames = list.files(projectDir, pattern = "*htseq-count", full.names = TRUE, recursive = TRUE)
  indexFile = read_lines(filenames[1])
  temp = {} ; for(i in 1:length(indexFile)){
    temp[[i]] = strsplit(indexFile[i], "\t")[[1]][1]
  }
  
  geneKey = temp
  index = as.numeric(grep(geneName, geneKey))
  
  cellNames = {} ; totalCounts = {}
  for(i in 1:length(filenames)){
    temp = read_lines(filenames[i])
    
    countValue = as.numeric(strsplit(temp[index], "\t")[[1]][2])
    if(countValue != 0 & !is.na(countValue)){
      #print("Nonzero")
      targetGene = strsplit(temp[index], "\t")[[1]][1]
      if(targetGene == geneName){
        totalCounts[[i]] = countValue
        cellNames[[i]] = strsplit(filenames[i], "*counts/|_S")[[1]][2]
      }
    }
  }
  
  geneCounts = list(totalCounts, cellNames)
  print(paste0("Counts loaded for: ", geneName))
  return(geneCounts)
  
}

loadGeneCounts = function(dashList){
  
  lil_dashGenes = {}
  for(i in 1:length(dashList)){
    geneCounts = countsByGene(projectDir, dashList[i])
    lil_dashGenes[[i]] = geneCounts
  }
  return(lil_dashGenes)
}


lil_gene_todf = function(lil_dashGenes, dashList){
  
  genesData = {}
  for(i in 1:length(geneList)){
    genesData$Sample_ID[[i]] = unlist(geneCounts[[i]][2])
    genesData$geneCounts[[i]] = as.numeric(unlist(geneCounts[[i]][1]))
    genesData$geneName[[i]] = geneList[i]
  }
  df = data_frame()
  for(i in 1:length(geneList)){
    if (!is.null(genesData$Sample_ID[[i]])){
      tmp = data_frame("geneCounts" = as.numeric(unlist(genesData$geneCounts[[i]])), 
                       "Sample_ID" = unlist(genesData$Sample_ID[[i]]), 
                       "geneName" = geneList[i])
      df = rbind(df, tmp)} else {print(paste0('is.null ', geneList[i]))}
  }
  
  return(df)
  
}


map_whatever = function(projectDir, data, column_name, log2 = FALSE, savePlot = FALSE){
  runID = strsplit(projectDir,"00_project_raw_data/")[[1]][2]
  plotTitle = paste0(column_name, ": ", runID)
  
  
  print(column_name[1])
  if(log2 == TRUE){
    fill_value = log2(data[[column_name[1]]])
    column_name = paste0("log2_",column_name)
  } else {fill_value = data[[column_name[1]]]}
  
  heatmap = ggplot() +
    geom_raster(data, mapping = aes(x = col, y = row, fill = fill_value)) +
    scale_x_discrete(limits = c(1:24)) + scale_y_discrete(limits = rev(levels(data$row))) + coord_equal() +
    scale_fill_viridis(guide_legend(title = column_name)) +
    theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank()) + 
    facet_wrap(~data$Description) + 
    ggtitle(plotTitle)
  
  if(savePlot == TRUE){
    dir.create(paste0(projectDir, "/figures/heatmaps/"), showWarnings = FALSE, recursive = TRUE)
    plotsPath = paste0(projectDir, "/figures/heatmaps/")
    ggsave(heatmap, paste0(plotsPath, "/heatmap_", runID, "_", column_name, ".png"))
    print(paste0("Saved plot to ", plotsPath))
  }
  
  return(heatmap)
}









