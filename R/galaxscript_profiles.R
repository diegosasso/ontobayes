#' Builds all input files and an script to run all analyses from ontobayes_profiles in Galax. 
#'
#' @param tags list, profiles with anatomy ontology terms to query. The list of profiles should be organized as indicated in the examples.
#' @param runs integer, the number of independent runs of mcmc used in the MrBayes commands block. Default is set to 2.
#' @param burnin integer, the absolute number of trees as the burn-in used in the mcmc. Default is set to 251.
#' @param OS character, sets the operational system to produce batch files for. Options are "windows" and "unix". Default is set to "windows".
#' @param filename character, name of the output file.
#'
#' @return All necessary input files, organized in subfolders inside the GALAX/INPUT folder, 
#'   and an script to run all Galax analyses. 
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
#' # Creates bat file to run all analyses on Galax #
#' galaxscript_profiles(tags = tags, runs = 2, burnin = 251, OS = "windows", filename = "galaxprof1")
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
#' # Creates bat file to run all analyses on Galax #
#' galaxscript_profiles(tags = tags, runs = 2, burnin = 251, OS = "windows", filename = "galaxprof2")
#' }
#' 
#' @export
galaxscript_profiles <- function(tags = tags, runs = 2, burnin = 251, OS = "windows", filename = "galaxscript")
{

  # Load all query terms #
  tags <- tags
  
  # Create a block to store some variables #
  block <- character()
  
  for(i in 1:length(tags)){
    
    # Create input and output folders and subfolders for Galax #
    dir.create(file.path("GALAX", "INPUT",  paste0(i, "_", names(tags)[i])), recursive = T)
    dir.create(file.path("GALAX", "OUTPUT", paste0(i, "_", names(tags)[i])), recursive = T)
    
    # Create some variables #
    x <- character()
    y <- character()
    
    if(OS == "windows"){
      
      for(k in 1:runs){
        
        run <- paste0("NEXUS\\", i, "_", names(tags)[i], "\\", 1:length(tags[[i]]), "_", names(tags[[i]]), ".nex.run", k, ".t")  
        
        # Write among-partitions files #
        write(paste(run, sep = "\n"), file = paste0("GALAX/INPUT/", i, "_", names(tags)[i], "/", names(tags)[i], "_run", k,".txt"))
        
        # Write first lines of among-partitions galax commands for the bat file #
        x <- c(x, paste0("galax -l GALAX\\INPUT\\", i, "_", names(tags)[i], "\\", names(tags)[i], "_run", k,".txt -s ", 
                         burnin, " -g 1 -o GALAX\\OUTPUT\\", i, "_", names(tags)[i], "\\", "galax_", names(tags)[i], "_run", k))
        
      }
      
    }
    
    if(OS == "unix"){
      
      for(k in 1:runs){
        
        run <- paste0("NEXUS/", i, "_", names(tags)[i], "/", 1:length(tags[[i]]), "_", names(tags[[i]]), ".nex.run", k, ".t")  
        
        # Write among-partitions files #
        write(paste(run, sep = "\n"), file = paste0("GALAX/INPUT/", i, "_", names(tags)[i], "/", names(tags)[i], "_run", k,".txt"))
        
        # Write first lines of among-partitions galax commands for the sh file #
        x <- c(x, paste0("galax --listfile GALAX/INPUT/", i, "_", names(tags)[i], "/", names(tags)[i], "_run", k,".txt --skip ", 
                         burnin, " --outgroup 1 --outfile GALAX/OUTPUT/", i, "_", names(tags)[i], "/", "galax_", names(tags)[i], "_run", k))
        
      }
      
    }
    
    for(j in 1:length(tags[[i]])){
      
      if(OS == "windows"){
        
        z <- character()
        
        for(u in 1:runs){
          
		  # Write galax command lines for all runs #
          z <- c(z, paste0("NEXUS\\", i, "_", names(tags)[i], "\\", j, "_", names(tags[[i]][j]), ".nex.run", u, ".t"))
          
        }
        
		# Write among-runs files #
        write(paste(z, sep = "\n"), file = paste0("GALAX/INPUT/", i, "_", names(tags)[i], "/", j, "_", names(tags[[i]][j]), ".txt"))
        
        # Write all lines of among-runs galax commans for the bat file #
        galaxline <- paste0("galax -l GALAX\\INPUT\\", i, "_", names(tags)[i], "\\", j, "_", names(tags[[i]][j]), ".txt -s ", 
                            burnin, " -g 1 -o GALAX\\OUTPUT\\", i, "_", names(tags)[i], "\\", "galax_", j, "_", names(tags[[i]][j]))
        
        y <- c(y, galaxline)
        
      }
      
      if(OS == "unix"){
        
        z <- character()
        
        for(u in 1:runs){
          
		  # Write galax command lines for all runs #
          z <- c(z, paste0("NEXUS/", i, "_", names(tags)[i], "/", j, "_", names(tags[[i]][j]), ".nex.run", u, ".t"))
          
        }
        
		# Write among-runs files #
        write(paste(z, sep = "\n"), file = paste0("GALAX/INPUT/", i, "_", names(tags)[i], "/", j, "_", names(tags[[i]][j]), ".txt"))
        
        # Write all lines of among-runs galax commands for the sh file #
        galaxline <- paste0("galax --listfile GALAX/INPUT/", i, "_", names(tags)[i], "/", j, "_", names(tags[[i]][j]), ".txt --skip ", 
                            burnin, " --outgroup 1 --outfile GALAX/OUTPUT/", i, "_", names(tags)[i], "/", "galax_", j, "_", names(tags[[i]][j]))
        
        y <- c(y, galaxline)
        
      }	  
      
    }
    
	# Create the batch block #
    block <- c(block, "", x, y)
    
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
