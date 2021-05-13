#' Builds script to run all analyses from ontobayes_profiles. 
#'
#' @param tags list, profiles with anatomy ontology terms to query. The list of profiles should be organized as indicated in the examples.
#' @param OS character, sets the operational system to produce batch files for. Options are "windows" and "unix". Default is set to "windows".
#' @param filename character, name of the output file.
#'
#' @return A script to run all MrBayes batch files inside NEXUS subfolders.
#'
#' @examples
#' \dontrun{
#' # HAO example #
#' # List of profiles with terms to query #
#' tags <- list(
#'   c("HAO:0000506", "HAO:0000513", "HAO:0000453"),
#'   c("HAO:0000234", "HAO:0001003"),
#'   c("HAO:0000874", "HAO:0000583", "HAO:0000630"))
#' 
#' # Set profile names and terms #
#' names(tags) <- c("mouthparts_profile1", "head_profile2", "mesosoma_profile3")
#' names(tags[[1]]) <- c("mandible", "maxilla", "labium")
#' names(tags[[2]]) <- c("cranium", "tentorium")
#' names(tags[[3]]) <- c("prothorax", "mesothorax", "metathorax")
#'
#' # Create MrBayes batch files to run all analyses #
#' bayescript_profiles(tags = tags, OS = "windows", filename = "bayesprof1")
#'
#' # UBERON example #
#' # List of profiles with terms to query #
#' tags <- list(
#'   c("UBERON_0002244", "UBERON_0002397", "UBERON_0004742"),
#'   c("UBERON_2000658", "UBERON_2000488"),
#'   c("UBERON_0000151", "UBERON_0000152", "UBERON_0003097"))
#' 
#' # Set profile names and terms #
#' names(tags) <- c("jaw_profile1", "pharynx_profile2", "fins_profile3")
#' names(tags[[1]]) <- c("premaxilla", "maxilla", "dentary")
#' names(tags[[2]]) <- c("epibranchial_bone", "ceratobranchial_bone")
#' names(tags[[3]]) <- c("pectoral_fin", "pelvic_fin", "dorsal_fin")
#'
#' # Create MrBayes batch files to run all analyses #
#' bayescript_profiles(tags = tags, OS = "windows", filename = "bayesprof2")
#' }
#'
#' @export
bayescript_profiles <- function(tags = tags, OS = "windows", filename = "bayescript")
{

  if(OS == "windows"){
  
	# Create headings for the bat file #
    echline <- "@echo off"
    runline <- "mb batch.nex"
    retline <- "cd ..\\.."
    
	# Create some variables #
    x <- character()
    y <- character()
    z <- character()
    
    for(i in 1:length(tags)){
	  
	  # Set all command lines #
      copyline <- paste0("copy bin\\mb.exe NEXUS\\", i, "_", names(tags)[i])
      middline <- paste0("cd NEXUS\\", i, "_", names(tags)[i])
      middline <- c(middline, runline, retline)
      deltline <- paste0("del NEXUS\\", i, "_", names(tags)[i], "\\mb.exe")
      
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
    
    for(i in 1:length(tags)){
	  
	  # Set all command lines #
      middline <- paste0("cd NEXUS/", i, "_", names(tags)[i])
      middline <- c(middline, runline, retline)
      
	  # Join and organize command lines #
      x <- c(x, middline, "")
	  
    }
    
	# Write batch file for Unix #
    write(c(binline, "", x), file = paste0(filename, ".sh"))
	
  }
  
}
