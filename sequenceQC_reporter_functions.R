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

#Load samplesheet for data generation and plotting. This function can be called separately for csv-only (data) and png-only (plot) generation. 
loadSamplesheet = function(projectDir){
  
  dirs = list.dirs(projectDir, recursive = FALSE)
  
  runID = list(strsplit(projectDir,"runs/")[[1]][2])
  for(i in 1:length(projectDir)){
    
    temp = list.files(dirs[2], full.names = TRUE, recursive = TRUE)
    
    samplesheet = read_csv(temp, col_names = TRUE, col_types = cols(), skip = 1)
    
    NUMBER_OF_PLATES = nrow(samplesheet)/384
    
    samplesheet$well_index = generate_wells(NUMBER_OF_PLATES)[1:nrow(samplesheet)]
    samplesheet$row = as.factor(samplesheet$well_index %>% str_extract(pattern = "^[A-Z]{1}"))
    samplesheet$col = as.factor(samplesheet$well_index %>% str_extract(pattern = "[0-9]{1,2}"))
    #samplesheet$Sample_Project = as.vector(runID[i])
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
  
  dirs = list.dirs(projectDir, recursive = FALSE)
  
  runID = {}
  for(i in projectDir){
    
    runID[[i]]= strsplit(i,"runs/")[[1]][2]
    
  }
  
  
  laneBarcodeList = {}
  for(i in 1:length(projectDir)){
    
    barcodeFile = list.files(dirs[1], pattern = "laneBarcode.html", full.names = TRUE, recursive = TRUE)
    
    for(j in barcodeFile){
      
      
      sampleList = {}
      barcodeList = {}
      clusterList = {}
      
      temp = read_lines(j, skip = 42)
      
      sampleList[[j]] = temp[seq(2, length(temp), 14)]
      barcodeList[[j]] = temp[seq(3, length(temp), 14)]
      clusterList[[j]] = temp[seq(4,length(temp),14)]
      
      
    }
    
    samples = {}
    clusters = {}
    barcodes = {}
    
    for(j in barcodeFile){
      
      samples[[j]] = unlist(sampleList[[j]])
      clusters[[j]] = unlist(clusterList[[j]])
      barcodes[[j]] = unlist(barcodeList[[j]])
      
    }
    
    laneBarcodeReport = list(unlist(samples), unlist(clusters), unlist(barcodes))
    
    for (x in 1:length(laneBarcodeReport)){laneBarcodeReport[[x]] = str_extract(laneBarcodeReport[[x]], "(?<=>)[^<]+")}
    
    laneBarcodeReport[[2]] = as.numeric(gsub(",", "", laneBarcodeReport[[2]]))
    laneBarcode = data_frame(Sample_ID = laneBarcodeReport[[1]], clusters = laneBarcodeReport[[2]], barcodes = laneBarcodeReport[[3]])
    
    #laneBarcodeList[[i]] = laneBarcode
    print(paste0("Demultiplexation report loaded for sequencing run: ", runID))
    
  }
  
  return(laneBarcode)
  
}

#Return key-value pair for yourParameter
return.parameterKey = function(projectDir, yourParameter){
  
  filenames = list.files(projectDir, pattern = "*.out", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  temp = read_lines(filenames[1])
  allParameters = {}
  
  for(i in 1:length(temp)){
    
    allParameters[[i]] = strsplit(temp[i], " \\|")[[1]][1]
    allParameters = gsub("^\\s+","", allParameters)
  }
  
  names(allParameters) = allParameters
  keys = {}
  keys$index = match(yourParameter, allParameters)
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

loadLogFile = function(projectDir, yourParameter){
  
  filenames = list.files(projectDir, pattern = "*.out", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  runID = list(strsplit(projectDir,"runs/")[[1]][2])
  for(i in 1:length(filenames)){runID$sampleID[[i]] = strsplit(filenames[i], "logs/|_S")[[1]][2]}
  runID$parameterKey = return.parameterKey(projectDir, yourParameter)
  
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
  
  log = {} ; for(i in 1:length(runID$parameterData)){log[[i]] = as.numeric(unlist(runID$parameterData[i]))}
  names(log) = unlist(runID$parameterKey$keynames)
  log$Sample_ID= unlist(runID$sampleID)
  df = as_tibble(log[1:length(log)])
  
  df = left_join(samplesheet, df, by = "Sample_ID", sort = FALSE)
  #write_csv(df, paste0(projectDir,"/", runID,"_logfinalout.csv"))
  return(df)
}


#This function generates heatmaps (and other plots) from your FASTQ data. This function currently calls the get_FASTQ_csv function, but should also be able to get called separately to accept the csv generated by get_FASTQ_csv.
heatmap.laneBarcode = function(projectDir, samplesheet, laneBarcode){
  runID = strsplit(projectDir,"runs/")[[1]][2]
  laneBarcodeReport = left_join(samplesheet, laneBarcode, by = "Sample_ID", sort = FALSE)
  write_csv(x = laneBarcodeReport, path = paste0(projectDir,"/", runID,"_laneBarcodeReport.csv"))
  
  print(paste0("Wrote ", runID, ".csv to ", projectDir,"/"))
  
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
    theme(strip.text = element_blank(),
          plot.title = element_text(size = 9, hjust = 0.5),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
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
    theme(strip.text = element_blank(),
          plot.title = element_text(size = 9, hjust = 0.5),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    facet_wrap(~data$Description) +
    ggtitle(plotTitle)
  ggsave(paste0(projectDir,"/", runID, "_laneBarcode_unscaled_heatmap.png"))
  print("Saved heatmap for unscaled clusters")
  
  
  heatmaps = list(log2_heatmap, unscaled_heatmap)
  
  return(heatmaps)
  
}

#This function generates heatmaps (and other plots) from your STAR data. It's important to maintain the same YOUR_PARAMETER call. This function currently calls the get_STAR_csv function, but should also be able to get called separately to accept the csv generated by get_STAR_csv.

plot_STAR_output = function(starPARAMS, PROJECT, YOUR_PARAMETER){
  
  starPARAM_plots = {}
  heatmap_list = {}
  grid_list = {}
  
  for(x in 1:length(starPARAMS)){
    print(paste0("project", x))
    
    data = starPARAMS[x][[1]][[1]]
    PROJECT_TITLE = data$Sample_Project[[1]]
    #PROJECT_TITLE = strsplit(PROJECT_INFO[[x]][1][[1]], "runs/[a-zA-Z0-9_-]{1,}/")[[1]][2]
    
    log2_heatmap = {}
    
    for(i in 1:length(YOUR_PARAMETER)){
      print(paste0("parameter", i))
      
      parameter_to_plot = i + 13
      data[[parameter_to_plot]] = as.numeric(data[[parameter_to_plot]])
      
      PLOT_TITLE = paste0("STAR OUTPUT for ", PROJECT_TITLE, " PARAMETER: ", parameter_to_plot)
      
      log2_heatmap[[i]] = ggplot() +
        geom_raster(data,
                    mapping = aes(x = col, y = row, fill = log2(data[[parameter_to_plot]]))) +
        scale_x_discrete(limits = c(1:24)) +
        scale_y_discrete(limits = rev(levels(data$row))) +
        coord_equal() +
        guides(fill = FALSE) +
        scale_fill_viridis(guide_legend(title = paste0("log2 ", parameter_to_plot))) +
        theme(strip.text = element_blank(),
              plot.title = element_text(size = 5, hjust = 0.5),
              legend.title = element_text(size = 5),
              legend.text = element_text(size = 5),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank()) +
        facet_wrap(~data$Description) +
        ggtitle(PLOT_TITLE)
      
      print(log2_heatmap[[i]])
      
      ggsave(paste0(PROJECT[[1]][x],"/",PROJECT[[2]][x], "_", parameter_to_plot,"_STARreport_log2heatmap.png"))
    }
    
    grid = grid.arrange(grobs = log2_heatmap, 
                        ncol = round(0.5*length(YOUR_PARAMETER), digits = 0), 
                        top = textGrob(paste0(PROJECT[[2]][x], "_", "ALL_PARAMETERS","_STARreport_log2heatmap")), 
                        gp=gpar(fontsize=2), padding = unit(1, "cm"))
    #print(grid)
    
    grid_list[[x]] = grid
    
    ggsave(paste0(PROJECT[[1]][x],"/",PROJECT[[2]][x], "_", "ALL_PARAMETERS","_STARreport_log2heatmap.png"), grid_list[[x]])
    
  }
  
  starPARAM_plots = list(heatmaps = heatmap_list, grids = grid_list)
  
  return(starPARAM_plots)
  
}



#Streamline your FASTQ report output. A start-to-finish wrapper function that streamlines your data output.

getyour_FASTQ_data = function(PROJECT_DIR){
  
  fastqPARAMS = get_BCL2FASTQ_csv(
    
    setup_samplesheet(
      load_projects(PROJECT_DIR))
    
  )
  
  
  FASTQ_PLOTS = plot_fastq_output(
    
    get_BCL2FASTQ_csv(setup_samplesheet(load_projects(PROJECT_DIR))), 
    load_projects(PROJECT_DIR)
    
  )
  
}



#Streamline your STAR statistics output


getyour_STAR_data = function(PROJECT_DIR, YOUR_PARAMETER){
  starPARAMS = get_STAR_csv(
    setup_samplesheet(
      load_projects(PROJECT_DIR)), YOUR_PARAMETER
  )
  
  STAR_PLOTS = plot_STAR_output(
    
    get_STAR_csv(setup_samplesheet(load_projects(PROJECT_DIR)), YOUR_PARAMETER), 
    load_projects(PROJECT_DIR),
    YOUR_PARAMETER
    
  )
}

