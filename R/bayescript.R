#' Builds script to run all analyses from ontobayes. 
#'
#' @param ids data.frame, matches between characters and anatomical entities. 
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#' @param OS character, sets the operational system to produce batch files for. Options are "windows" and "unix". Default is set to "windows".
#' @param filename character, name of the output file.
#'
#' @return A script to run all MrBayes batch files inside NEXUS/ALL_PART subfolders. 
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
#' # Create MrBayes batch files to run all analyses #
#' bayescript(ids = ID, OS = "windows", filename = "bayescript1")
#' 
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # Create MrBayes batch files to run all analyses #
#' bayescript(ids = ID, OS = "windows", filename = "bayescript2")
#' }
#'
#' @export
bayescript <- function(ids = ids, OS = "windows", filename = "bayescript")
{

  # Import ontology identifiers and term names to query #
  ids <- ids

  # Adjust ontology term names #
  ids <- gsub(ids[,3], pattern = " ", replacement = "_")
  ids <- gsub(ids, pattern = "/", replacement = "-")
  
  # Get all unique terms included in the dataset #
  ids <- unique(ids)
  
  if(OS == "windows"){
  
	# Create headings for the bat file #
    echline <- "@echo off"
    runline <- "mb batch.nex"
    retline <- "cd ..\\..\\.."
    
	# Create some variables #
    x <- character()
    y <- character()
    z <- character()
    
    for(i in 1:length(ids)){
	  
	  # Set all command lines #
      copyline <- paste0("copy bin\\mb.exe NEXUS\\ALL_PART\\", i, "_", ids[i])
      middline <- paste0("cd NEXUS\\ALL_PART\\", i, "_", ids[i])
      middline <- c(middline, runline, retline)
      deltline <- paste0("del NEXUS\\ALL_PART\\", i, "_", ids[i], "\\mb.exe")
      
	  # Join and organize command lines #
      x <- paste0(x, copyline, sep = "\n")
      y <- c(y, middline, "")
      z <- paste0(z, deltline, sep = "\n")
	  
    }
    
	# Write batch file for Windows #
    write(c(echline, "", x, y, z), file = paste0(filename, ".bat"))
	
  }
  
  if(OS == "unix"){
  
	# Create headings for the batch file #
    binline <- "#!/bin/bash"
    runline <- "mb batch.nex"
    retline <- "cd ../.."
    
    x <- character()
    
    for(i in 1:length(ids)){
	  
	  # Set all command lines #
      middline <- paste0("cd NEXUS/ALL_PART/", i, "_", ids[i])
      middline <- c(middline, runline, retline)
      
	  # Join and organize command lines #
      x <- c(x, middline, "")
	  
    }
    
	# Write batch file for Unix #
    write(c(binline, "", x), file = paste0(filename, ".sh"))
	
  }
  
}
