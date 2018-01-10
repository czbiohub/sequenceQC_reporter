
#Load R libraries

lbry = c("tidyverse", "stringr", "viridis", "gridExtra", "grid", "SparseM")
lapply(lbry, install.packages, character.only=TRUE)
lapply(lbry, require, character.only=TRUE)


#Load projects into environment. This function lists the projects in the directory with S3-synced files called by the bash script.

#Get project names
load_projects = function(PROJECT_DIR){
  
  directory_names = list.dirs(PROJECT_DIR, recursive = FALSE)
  
  projects = {}
  for(i in directory_names){
    projects[[i]]= strsplit(i,"runs/")[[1]][2]
  }
  
  return(list(directory_names, projects))
}



#Function to generate alphanumeric index for 384-well plates. This function can be used alone or within subsequent functions. If used alone, input is NUMBER_OF_PLATES. 

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



#Setup your samplesheet for data generation and plotting. This function can be called separately for csv-only (data) and png-only (plot) generation. 

setup_samplesheet = function(PROJECT){
  
  setup_samplesheet_list = {}
  for(i in 1:length(PROJECT[1][1][[1]])){
    
    path = paste0(PROJECT[[1]][i],"/", PROJECT[[2]][i])
    dir.create(paste0(path, "/output/"), showWarnings = FALSE)
    
    samplesheet = read_csv(paste0(paste0(PROJECT[[1]][i], "/sample-sheets/", PROJECT[[2]][i],".csv")), col_names = TRUE, skip = 1)
    
    NUMBER_OF_PLATES = nrow(samplesheet)/384
    
    samplesheet$well_index = generate_wells(NUMBER_OF_PLATES)[1:nrow(samplesheet)]
    samplesheet$row = as.factor(samplesheet$well_index %>% str_extract(pattern = "^[A-Z]{1}"))
    samplesheet$col = as.factor(samplesheet$well_index %>% str_extract(pattern = "[0-9]{1,2}"))
    samplesheet$Sample_Project = PROJECT[[2]][i]
    samplesheet$Description = rep(1:NUMBER_OF_PLATES, each=384)
    
    if(is.integer(samplesheet$`Sample_ID`) == TRUE){
      samplesheet$`Sample_ID` = as.character(samplesheet$`Sample_ID`)
    }
    
    setup_samplesheet_list[[i]] = list(path, samplesheet)
  }
  
  return(setup_samplesheet_list)
  
}




#Generate your FASTQ stats - this function scrapes Illumina bcl2fastq html report files and sorts data like "passed-filter" clusters, samples, and barcodes.

get_BCL2FASTQ_csv = function(PROJECT_INFO){
  
  fastq_report_list = {}
  
  for(i in 1:length(PROJECT_INFO)){
    
    filenames = paste0(strsplit(PROJECT_INFO[[i]][1][[1]], "/[a-zA-Z0-9_-]{1,}$"),"/reports")
    
    filenames = list.files(filenames, pattern = "laneBarcode.html", full.names = TRUE, recursive = TRUE)
    
    for(j in filenames){
      
      clusterList = {}
      sampleList = {}
      barcodeList = {}
      
      temp = read_lines(j, skip = 42)
      
      sampleList[[j]] = temp[seq(2, length(temp), 14)]
      clusterList[[j]] = temp[seq(4,length(temp),14)]
      barcodeList[[j]] = temp[seq(3, length(temp), 14)]
      
    }
    
    samples = {}
    clusters = {}
    barcodes = {}
    
    for(j in filenames){
      
      samples[[j]] = unlist(sampleList[[j]])
      clusters[[j]] = unlist(clusterList[[j]])
      barcodes[[j]] = unlist(barcodeList[[j]])
      
    }
    
    fastq_report = list(unlist(samples), unlist(clusters), unlist(barcodes))
    
    for (x in 1:length(fastq_report)){fastq_report[[x]] = str_extract(fastq_report[[x]], "(?<=>)[^<]+")}
    
    fastq_report[[2]] = as.numeric(gsub(",", "", fastq_report[[2]]))
    fastq_report = data_frame(Sample_ID = fastq_report[[1]], clusters = fastq_report[[2]], barcodes = fastq_report[[3]])
    fastq_report = left_join(PROJECT_INFO[i][[1]][[2]], fastq_report, by = "Sample_ID", sort = FALSE)
    
    write_csv(fastq_report, paste0(PROJECT_INFO[[i]][1][[1]],"_bcl2fastq_report.csv"))
    
    fastq_report_list[[i]] = list(fastq_report)
    
  }
  
  return(fastq_report_list)
}



#Generate your STAR stats - This input for his function is any single line or multiple lines as a PARAMETER (i.e. "Mapped Reads", "Percent Mapped Reads", "Insert Average Length") from the "log.final.out" file. This file a standard log generated afer successful STAR alignment.

get_STAR_csv = function(PROJECT_INFO, YOUR_PARAMETER){
  
  star_report_list = {}
  
  for(x in 1:length(PROJECT_INFO)){
    
    filenames = paste0(strsplit(PROJECT_INFO[[x]][1][[1]], "/[a-zA-Z0-9_-]{1,}$"),"/star_logs")
    
    
    filenames = list.files(filenames, pattern = "*.out", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    
    sample_id = {}
    
    for(i in filenames){
      sample = strsplit(i, "./results*")[[1]][1]
      sample = strsplit(sample, "logs/")[[1]][2]
      sample = strsplit(sample, "_S")[[1]][1]
      sample_id[[i]] = strsplit(sample, ".mus*")[[1]][1]
    }
    
    logfinalout = {}
    logfinalout_parameters = {}
    logfinalout_data = {}
    
    for(j in 1:length(filenames)){
      
      if(is.na(j)){
        
        print("Skipping missing")
        break
        
      } else {
        temp = read_lines(filenames[j])
        for(i in 1:length(YOUR_PARAMETER)){
          k = YOUR_PARAMETER[i]
          logfinalout_parameters[[i]] = str_split(temp[k], " \\|")[[1]][1]
          logfinalout_parameters = gsub("^\\s+","", logfinalout_parameters)
          logfinalout_data[[i]] = str_split(temp[k], "\t")[[1]][2]
          logfinalout_data = gsub("%","", logfinalout_data)
        }
        logfinalout[[j]] = list(logfinalout_parameters, logfinalout_data)
      }
    }
    
    secondlevel = {}
    firstlevel = {}
    
    for(i in 1:length(YOUR_PARAMETER)){
      for(j in 1:length(logfinalout)){secondlevel[[j]] = logfinalout[j][[1]][[2]][i]}
      firstlevel[[i]] = secondlevel
    }
    
    SAMPLE_ID = unlist(sample_id)
    star_report = data_frame("Sample_ID" = SAMPLE_ID)
    star_report[ ,ncol(star_report) + 1:length(YOUR_PARAMETER)] = NA
    
    for(i in 1:length(YOUR_PARAMETER)){
      star_report[i+1] = firstlevel[i]
      colnames(star_report)[i+1]= logfinalout[[1]][[1]][i]
    }
    
    star_report = left_join(PROJECT_INFO[x][[1]][[2]], star_report, by = "Sample_ID", sort = FALSE)
    write_csv(star_report, paste0(PROJECT_INFO[[x]][1][[1]],"_star_report.csv"))
    star_report_list[[x]] = list(star_report)
  }
  
  return(star_report_list)
}




#This function generates heatmaps (and other plots) from your FASTQ data. This function currently calls the get_FASTQ_csv function, but should also be able to get called separately to accept the csv generated by get_FASTQ_csv.

plot_fastq_output = function(fastqPARAMS, PROJECT){
  
  fastqPARAMS_plots = {}
  heatmap_list = {}
  grid_list = {}
  
  
  for(x in 1:length(fastqPARAMS)){
    print(paste0("project", x))
    
    data = fastqPARAMS[x][[1]][[1]]
    PROJECT_TITLE = data$Sample_Project[[1]]
    data[[14]] = as.numeric(data[[14]])
    
    log2_heatmap = {}
    heatmap = {}
    
    PLOT_TITLE = paste0("BCL2FASTQ OUTPUT for ", PROJECT_TITLE, " PARAMETER: ", "PASSED-FILTER Clusters")
    
    log2_heatmap[[x]] = ggplot() +
      geom_raster(data, 
                  mapping = aes(x = col, y = row, fill = log2(data$clusters))) +
      scale_x_discrete(limits = c(1:24)) +
      scale_y_discrete(limits = rev(levels(data$row))) +
      coord_equal() +
      #guides(fill = FALSE) +
      scale_fill_viridis(guide_legend(title = "log2")) +
      theme(strip.text = element_blank(),
            plot.title = element_text(size = 5, hjust = 0.5),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()) +
      facet_wrap(~data$Description) +
      ggtitle(PLOT_TITLE)
    ggsave(paste0(PROJECT[[1]][x],"/",PROJECT[[2]][x], "_", "PASSED-FILTER Clusters","_miseq_bcl2fastqreport_log2heatmap.png"))
    
    print(log2_heatmap[[x]])
    
    heatmap[[x]] = ggplot() +
      geom_raster(data, 
                  mapping = aes(x = col, y = row, fill = data$clusters)) +
      scale_x_discrete(limits = c(1:24)) +
      scale_y_discrete(limits = rev(levels(data$row))) +
      coord_equal() +
      #guides(fill = FALSE) +
      scale_fill_viridis(guide_legend(title = "unscaled")) +
      theme(strip.text = element_blank(),
            plot.title = element_text(size = 5, hjust = 0.5),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()) +
      facet_wrap(~data$Description) +
      ggtitle(PLOT_TITLE)
    ggsave(paste0(PROJECT[[1]][x],"/",PROJECT[[2]][x], "_", "PASSED-FILTER Clusters","_miseq_bcl2fastqreport_heatmap.png"))
    
    print(heatmap[[x]])
    
  }
  
  heatmap_list = list(log2_heatmap, heatmap)
  fastqPARAMS_plots = list(all_my_heatmaps = heatmap_list)
  
  return(fastqPARAMS_plots)
  
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
