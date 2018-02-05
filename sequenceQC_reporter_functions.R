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
sortSheetData = function(samplesheet, whatever){
  
  df = left_join(samplesheet, whatever, by = "Sample_ID", sort = FALSE)
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

#Parse demultiplex html report files and sort data
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



#Parse read counts for RL's amplicon fastqs
parseAmpliconCounts = function(countFile){
  countList = {} ; sampleList = {}
  temp = read_lines(countFile)
  for(j in 1:length(temp)){
    countList[[j]] = strsplit(temp[j], "/|\"\t")[[1]][3]
    sampleList[[j]] = strsplit(temp[j], "/|\"\t")[[1]][2]
  }
  ampliconCounts = data_frame("readCounts" = unlist(countList), "Sample_ID" = unlist(sampleList))
  return(ampliconCounts)
}


#Load read counts for RL's amplicon fastqs
loadAmpliconCounts = function(projectDir){
  
  runID = {} ; for(i in projectDir){runID[[i]]= strsplit(i,"00_project_raw_data/")[[1]][2]}
  for(i in 1:length(projectDir)){
    countFile = list.files(projectDir, pattern = "readcount.txt", full.names = TRUE, recursive = TRUE)
    allAmps= data_frame()
    for(j in 1:length(countFile)){
      fn = parseAmpliconCounts(countFile[j])
      allAmps = rbind(allAmps, fn)
    }
  }
  return(allAmps)
}




#Load demultiplex html report data
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
    return(allLanes)
  }
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
  if(length(filenames) == 0){
    
    print(paste0("No `log.final.out` files were found for sequencing run: ", runID))
    break
    
  }
  
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
plot.laneBarcode = function(projectDir, samplesheet, laneBarcode){
  dir.create(paste0(projectDir, "/figures/"), showWarnings = FALSE)
  dir.create(paste0(projectDir, "/figures/laneBarcode/"), showWarnings = FALSE)
  plotsPath = paste0(projectDir, "/figures/laneBarcode/")
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
  ggsave(paste0(plotsPath, runID, "_laneBarcode_log2scale_heatmap.png"))
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
  ggsave(paste0(plotsPath, runID, "_laneBarcode_unscaled_heatmap.png"))
  print("Saved heatmap for unscaled clusters")
  
  heatmaps = list(log2_heatmap, unscaled_heatmap)
  
  return(heatmaps)
  
}

#Generate figures from STAR log data.

plot.starLog = function(projectDir, yourParameter, samplesheet, starLog){
  dir.create(paste0(projectDir, "/figures/"), showWarnings = FALSE)
  dir.create(paste0(projectDir, "/figures/star/"), showWarnings = FALSE)
  plotsPath = paste0(projectDir, "/figures/star/")
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
    ggsave(paste0(plotsPath, runID, "_", parameterName,"_starLogReport_log2scale_heatmap.png"))
    
    unscaled_heatmap = ggplot() +
      geom_raster(data,
                  mapping = aes(x = col, y = row, fill = data[[yourParameter[i]]])) +
      scale_x_discrete(limits = c(1:24)) + scale_y_discrete(limits = rev(levels(data$row))) + coord_equal() +
      scale_fill_viridis(guide_legend(title = paste0("unscaled"))) +
      theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
            legend.title = element_text(size = 7), legend.text = element_text(size = 7),
            axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
      facet_wrap(~data$Description) +
      ggtitle(plotTitle) +
      ggsave(paste0(plotsPath, runID, "_", parameterName,"_starLogReport_unscaled_heatmap.png"))
    
    
    log2_boxplot = ggplot(data) +
      geom_boxplot(mapping = aes(x = as.factor(Description), y = log2(data[[yourParameter[i]]]), fill = as.factor(Description))) +
      labs(x = "Plate number", y = parameterName) + 
      scale_fill_viridis(discrete = TRUE) + guides(fill = FALSE) +
      theme_classic() +
      theme(axis.text.x=element_text(angle = 90)) +
      ggtitle(plotTitle) +
      ggsave(paste0(plotsPath, runID, "_", parameterName, "_starLogReport_log2scale_boxplot.png"))
    
    unscaled_boxplot = ggplot(data) +
      geom_boxplot(mapping = aes(x = as.factor(Description), y = data[[yourParameter[i]]], fill = as.factor(Description))) +
      labs(x = "Plate number", y = parameterName) + 
      scale_fill_viridis(discrete = TRUE) + guides(fill = FALSE) +
      theme_classic() +
      theme(axis.text.x=element_text(angle = 90)) +
      ggtitle(plotTitle) +
      ggsave(paste0(plotsPath, runID, "_", parameterName, "_starLogReport_unscaled_boxplot.png"))
  }
  
  if(is.vector(data$`Uniquely mapped reads number`) == TRUE){
    
    if(is.vector(data$`Number of input reads`) == TRUE){
      ggplot(data) +
        geom_point(mapping = aes(x = log2(`Uniquely mapped reads number`), y = log2(`Number of input reads`)), alpha = 0.6, size = 0.2) +
        labs(x = "Uniquely mapped reads number", y = "Number of input reads") +
        facet_wrap(~Description) + ggtitle("Correlation of Unmapped Input to Mapped Unique") +
        theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
        ggsave(paste0(plotsPath, runID, "_", parameterName, "_starLogReport_log2scale_scatterplot.png"))
      
      if(is.vector(data$`Uniquely mapped reads %`) == TRUE){
        ggplot(data) +
          geom_point(mapping = aes(x = log2(`Uniquely mapped reads number`), y = `Uniquely mapped reads %`), alpha = 0.6, size = 0.2) +
          labs(x = "Uniquely mapped reads number", y = "Uniquely mapped reads %") +
          facet_wrap(~Description) + ggtitle("Correlation of Mapped Reads % to Mapped Unique") +
          theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
          ggsave(paste0(plotsPath, runID, "_", parameterName, "_starLogReport_log2scale_mappedReads2percent_scatterplot.png"))
        
        
        ggplot(data) +
          geom_point(mapping = aes(x = log2(`Number of input reads`), y = `Uniquely mapped reads %`), alpha = 0.6, size = 0.2) +
          labs(x = "Number of input reads", y = "Uniquely mapped reads %") +
          facet_wrap(~Description) + ggtitle("Correlation of Mapped Reads % to Unmapped Input") +
          theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
          ggsave(paste0(plotsPath, runID, "_", parameterName, "_starLogReport_log2scale_unmappedReads2percent_scatterplot.png"))
      }
    }
  }
  return()
}


#Parse htseq-count files to count by a given input for geneName (important: geneName must be a character string that exists in the htseq count file)
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
  print(paste0("Counts loaded for: ", geneName))
  return(geneCounts)
  
}

#Load a dataframe with counts for a list containing multiple geneNames indexed by 'Sample_ID'
loadGeneCounts = function(dashList){
  
  lil_dashGenes = {}
  for(i in 1:length(dashList)){
    geneCounts = countsByGene(projectDir, dashList[i])
    lil_dashGenes[[i]] = geneCounts
  }
  
  df = data_frame()
  for(i in 1:length(dashList)){
    genesData = data_frame("geneCounts" = as.numeric(unlist(lil_dashGenes[[i]][[1]])), "Sample_ID" = unlist(lil_dashGenes[[i]][[2]]), "geneName" = dashList[i])
    df = rbind(df, genesData)
  }
  
  print(paste0("Loaded gene counts indexed by samplesheet 'Sample_ID' for: ", dashList))
  
  return(df)
  
}

#Generate figures from htseq-count data.
plot.geneCounts = function(projectDir, samplesheet, geneCounts, dashList, starLog){
  dir.create(paste0(projectDir, "/figures/"), showWarnings = FALSE)
  dir.create(paste0(projectDir, "/figures/htseq-count/"), showWarnings = FALSE)
  plotsPath = paste0(projectDir, "/figures/htseq-count/")
  
  runID = strsplit(projectDir,"00_project_raw_data/")[[1]][2]
  data = data_frame()
  
  for(i in 1:length(dashList)){
    temp = geneCounts[grep(dashList[i], geneCounts$geneName), ]
    countsReport = left_join(samplesheet, temp, by = "Sample_ID", sort = FALSE)
    data = countsReport
    log2_heatmap = ggplot() +
      geom_raster(data, 
                  mapping = aes(x = col, y = row, fill = log2(data$geneCounts))) +
      scale_x_discrete(limits = c(1:24)) +
      scale_y_discrete(limits = rev(levels(data$row))) +
      coord_equal() +
      scale_fill_viridis(guide_legend(title = "log2")) +
      theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
            legend.title = element_text(size = 7), legend.text = element_text(size = 7),
            axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
      facet_wrap(~data$Description) + 
      ggtitle(paste0("Gene counts (log2) for: ", dashList[[i]])) +
      ggsave(paste0(plotsPath, runID, "_", dashList[[i]],"_log2scale_heatmap.png"))
    
    unscaled_heatmap = ggplot() +
      geom_raster(data, 
                  mapping = aes(x = col, y = row, fill = data$geneCounts)) +
      scale_x_discrete(limits = c(1:24)) +
      scale_y_discrete(limits = rev(levels(data$row))) +
      coord_equal() +
      scale_fill_viridis(guide_legend(title = "log2")) +
      theme(strip.text = element_blank(), plot.title = element_text(size = 9, hjust = 0.5),
            legend.title = element_text(size = 7), legend.text = element_text(size = 7),
            axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
      facet_wrap(~data$Description) + 
      ggtitle(paste0("Gene counts (unscaled) for: ", dashList[[i]])) +
      ggsave(paste0(plotsPath, runID, "_", dashList[[i]],"_unscaled_heatmap.png"))
    
    if(missing(starLog)){
      print("No STAR log data found. Plotting only gene counts data.")
    } else {
      
      print("Found STAR log data. Plotting genes and reads. ")
      genesReads = left_join(countsReport, starLog, by = "Sample_ID", sort = FALSE)
      data = genesReads
      
      if(is.vector(data$`Uniquely mapped reads number`) == TRUE){
        ggplot(data) +
          geom_point(mapping = aes(x = log2(`Uniquely mapped reads number`), y = log2(data$geneCounts)), alpha = 0.6, size = 0.2) +
          labs(x = "Uniquely mapped reads number", y = "Gene Counts") + 
          ggtitle(paste0("Correlation of gene counts (log2) for: ", dashList[[i]], " to Mapped Reads")) +
          facet_wrap(~Description) +
          theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
          ggsave(paste0(plotsPath, runID, "_", dashList[[i]], "_log2scale_Mapped_scatterplot.png"))
        
        
        
        if(is.vector(data$`Number of input reads`) == TRUE){
          data$ratioDashedGene = (data$geneCounts/data$`Number of input reads`)*100
          ggplot(data) +
            geom_point(mapping = aes(x = log2(`Number of input reads`), y = log2(data$geneCounts)), alpha = 0.6, size = 0.2) +
            labs(x = "Number of input reads", y = "Gene Counts") + 
            ggtitle(paste0("Correlation of gene counts (log2) for: ", dashList[[i]], " to Unmapped Reads")) +
            facet_wrap(~Description) +
            theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
            ggsave(paste0(plotsPath, runID, "_", dashList[[i]], "_Unmapped_log2scale_scatterplot.png"))
          
          ggplot(data) +
            geom_boxplot(mapping = aes(x = as.factor(Description), y = ratioDashedGene, fill = as.factor(Description))) +
            labs(x = "Plate number", y = "Gene/Total %") +
            scale_fill_viridis(discrete = TRUE) + guides(fill = FALSE) +
            ggtitle(paste0("Proportion of gene counts for: ", dashList[[i]], " to Unmapped Reads")) +
            theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
            ggsave(paste0(plotsPath, runID, "_", dashList[[i]], "_totalGene2input_unscaled_boxplot.png"))
          
          if(is.vector(data$`Uniquely mapped reads %`) == TRUE){
            ggplot(data) +
              geom_point(mapping = aes(x = `Uniquely mapped reads %`, y = log2(data$geneCounts)), alpha = 0.6, size = 0.2) +
              labs(x = "Uniquely mapped reads %", y = "Gene Counts") + 
              ggtitle(paste0("Correlation of gene counts (log2) for: ", dashList[[i]], " to Mapped Read %")) +
              facet_wrap(~Description) +
              theme_classic() + theme(axis.text.x=element_text(angle = 90)) +
              ggsave(paste0(plotsPath, runID, "_", dashList[[i]], "_MappedPercent_log2scale_scatterplot.png"))
          }
        }
      }
    }
  }
}





