#' Builds all input files and an script to run all analyses from ontobayes in Galax. 
#'
#' @param ids data.frame, matches between characters and anatomical entities. 
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#' @param runs integer, the number of independent runs of mcmc used in the MrBayes commands block. Default is set to 2.
#' @param burnin integer, the absolute number of trees as the burn-in used in the mcmc. Default is set to 251.
#' @param OS character, sets the operational system to produce batch files for. Options are "windows" and "unix". Default is set to "windows".
#' @param filename character, name of the output file.
#'
#' @return All necessary input files, organized in subfolders inside the GALAX/INPUT/all_partitions folder, 
#'   and an script to run all Galax analyses. 
#'
#' @examples
#' \dontrun{
#' # HAO example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("HAO:0000506",5), rep("HAO:0000513",5), rep("HAO:0000453",5), rep("HAO:0000234",5), 
#' rep("HAO:0001003",5), rep("HAO:0000874",5), rep("HAO:0000583",5), rep("HAO:0000630",5)),
#' c(rep("mandible",5), rep("maxilla",5), rep("labium",5), rep("cranium",5)
#' , rep("tentorium",5), rep("prothorax",5), rep("mesothorax",5), rep("metathorax",5))))
#'
#' # Creates bat file to run all analyses on Galax #
#' galaxscript(ids = ID, runs = 2, burnin = 251, OS = "windows", filename = "galaxruns1")
#'
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # Creates bat file to run all analyses on Galax #
#' galaxscript(ids = ID, runs = 2, burnin = 251, OS = "windows", filename = "galaxruns2")
#' }
#' 
#' @export
galaxscript <- function(ids = ids, runs = 2, burnin = 251, OS = "windows", filename = "galaxscript")
{
  
  # Adjust ontology term names #
  ids <- gsub(ids[,3], pattern = " ", replacement = "_")
  ids <- gsub(ids, pattern = "/", replacement = "-")
  
  # Get all unique terms included in the dataset #
  ids <- unique(ids)
  
  # Create input and output folders and subfolders for Galax #
  dir.create(file.path("GALAX", "INPUT",  "all_partitions"), recursive = T)
  dir.create(file.path("GALAX", "OUTPUT", "all_partitions"), recursive = T)
  
  # Create the batch block #
  block <- character()
  
  for(i in 1:length(ids)){
    
    if(OS == "windows"){
      
      z <- character()
      
      for(j in 1:runs){
        
        # Write galax command lines for all runs #
        z <- c(z, paste0("NEXUS\\ALL_PART\\", i, "_", ids[i], "\\", ids[i], ".nex.run", j, ".t"))
        
      }
      
      # Write among-runs files #
      write(paste(z, sep = "\n"), file = paste0("GALAX/INPUT/all_partitions/", ids[i], ".txt"))
      
      # Write all lines of among-runs galax commands for the bat file #
      galaxline <- paste0("galax -l GALAX\\INPUT\\all_partitions\\", ids[i], ".txt -s ", 
                          burnin, " -g 1 -o GALAX\\OUTPUT\\all_partitions\\", "galax_", ids[i])
      
      block <- c(block, galaxline)
      
    }
    
    if(OS == "unix"){
      
      z <- character()
      
      for(j in 1:runs){
        
        # Write galax command lines for all runs #
        z <- c(z, paste0("NEXUS/ALL_PART/", i, "_", ids[i], "/", ids[i], ".nex.run", j, ".t"))
        
      }
      
      # Write among-runs files #
      write(paste(z, sep = "\n"), file = paste0("GALAX/INPUT/all_partitions/", ids[i], ".txt"))
      
      # Write all lines of among-runs galax commands for the sh file #
      galaxline <- paste0("galax --listfile GALAX/INPUT/all_partitions/", ids[i], ".txt --skip ", 
                          burnin, " --outgroup 1 --outfile GALAX/OUTPUT/all_partitions/", "galax_", ids[i])
      
      block <- c(block, galaxline)
      
    }	  
    
  }  
  
  if(OS == "windows"){
    
    # Create the headings for the bat file #
    headline <- paste(
      "@echo off",
      "set PATH=bin;%PATH%",
      "galax.exe %*",
      sep = "\n")
    
    # Write bat file for Windows #
    write(c(headline, block), file = paste0(filename, ".bat"))
    
  }
  
  if(OS == "unix"){
    
    # Create the headings for the sh file #
    headline <- "#!/bin/bash"
    
    # Write sh file for Unix #
    write(c(headline, block), file = paste0(filename, ".sh"))
    
  }   
  
}  
