#' Builds all input files to run all analyses from resamp_bayes in Galax. 
#'
#' @param roots character, a named vector with the ontology IDs of the anatomical data subsets to produce MrBayes NEXUS files. 
#'   IDs and names should match those in the second and third columns of the ids data.frame, respectively.
#' @param nsamp integer, sets the number of samples for the resampling analysis.
#' @param runs integer, sets the number of independent runs of mcmc in the MrBayes commands block. Default is set to 2.
#' @param burnin integer, sets the absolute number of trees as the burn-in for the mcmc. Default is set to 251.
#' @param OS character, sets the operational system to produce batch files for. Options are "windows" and "unix". Default is set to "windows".
#' @param foldername character name of the output folder with all MrBayes NEXUS files.
#' @param filename character, name of the output file.
#'
#' @return All necessary input files, organized in the GALAX/INPUT/foldername folder, 
#'
#' @export
resamp_galax <- function(roots = roots, nsamp = 100, runs = 4, burnin = 251, OS = "windows", 
foldername = "resamp", filename = "galaxscript")
{

  # Create input and output folders and subfolders for Galax #
  suppressWarnings(dir.create(file.path("GALAX", "INPUT",  foldername), recursive = T))
  suppressWarnings(dir.create(file.path("GALAX", "OUTPUT", foldername), recursive = T))

  # Create the batch block #
  block <- character()

  for(i in 1:nsamp){
  
    for(j in 1:length(roots)){
	
	  if(OS == "windows"){
	  
	    z <- character()
		
		for(k in 1:runs){
		
		  # Write galax command lines for all runs #
		  z <- c(z, paste0("NEXUS\\RESAMP\\", foldername, "\\", names(roots)[j], "_resamp_", i, ".nex.run", k, ".t"))
		
		}
		
		# Write among-runs files #
		write(paste(z, sep = "\n"), file = paste0("./GALAX/INPUT/", foldername, "/", names(roots)[j], "_resamp_", i, ".txt"))
		
		# Write all lines of among-runs galax commands for the bat file #
		galaxline <- paste0("galax -l GALAX\\INPUT\\", foldername, "\\", names(roots)[j], "_resamp_", i, ".txt -s ", 
		burnin, " -g 1 -o GALAX\\OUTPUT\\", foldername, "\\", "galax_", names(roots)[j], "_resamp_", i)
		
		block <- c(block, galaxline)
	  
	  }
	  
	  if(OS == "unix"){
	  
	    z <- character()
		
		for(k in 1:runs){
		
		  # Write galax command lines for all runs #
		  z <- c(z, paste0("NEXUS/RESAMP/", foldername, "/", names(roots)[j], "_resamp_", i, ".nex.run", k, ".t"))
		
		}
		
		# Write among-runs files #
		write(paste(z, sep = "\n"), file = paste0("./GALAX/INPUT/", foldername, "/", names(roots)[j], "_resamp_", i, ".txt"))
		
		# Write all lines of among-runs galax commands for the sh file #
		galaxline <- paste0("galax --listfile GALAX/INPUT/", foldername, "/", names(roots)[j], "_resamp_", i, ".txt --skip ", 
		burnin, " --outgroup 1 --outfile GALAX/OUTPUT/", foldername, "/", "galax_", names(roots)[j], "_resamp_", i)
		
		block <- c(block, galaxline)
	  
	  }	  
	  
	}
  
  }

  if(length(roots) > 1){
    
	if(OS == "windows"){
	
	  # Set path names #
	  x <- paste0("NEXUS\\RESAMP\\", foldername, "\\")
	  y <- paste0("_resamp_", 1:nsamp, ".nex.run")
	  
	  for(v in 1:nsamp){
	  
	    for(u in 1:runs){
		
		  # Write among-partition files #
		  write(paste0(x, names(roots), y[v], u, ".t"), 
		  file = paste0("./GALAX/INPUT/", foldername, "/", paste0(names(roots), collapse = "_"), 
		  "_resamp_", v, "_run_", u,".txt"))
		  
		  # Write all lines of among-partition galax commands for the bat file #
		  block <- c(block, paste0("galax -l GALAX\\INPUT\\", foldername, "\\", paste0(names(roots), collapse = "_"), "_resamp_", v, "_run_", u, ".txt -s ", 
		  burnin, " -g 1 -o GALAX\\OUTPUT\\", foldername, "\\", "galax_", paste0(names(roots), collapse = "_"), "_resamp_", v, "_run_", u))
		
		}
	  
	  }
	
	}
	
	if(OS == "unix"){
	
	  # Set path names #
	  x <- paste0("NEXUS/RESAMP/", foldername, "/")
	  y <- paste0("_resamp_", 1:nsamp, ".nex.run")
	  
	  for(v in 1:nsamp){
	  
	    for(u in 1:runs){
		
		  # Write among-partition files #
		  write(paste0(x, names(roots), y[v], u, ".t"), 
		  file = paste0("./GALAX/INPUT/", foldername, "/", paste0(names(roots), collapse = "_"), "_resamp_", v, "_run_", u,".txt"))
		  
		  # Write all lines of among-partition galax commands for the bat file #
		  block <- c(block, paste0("galax --listfile GALAX/INPUT/", foldername, "/", paste0(names(roots), collapse = "_"), "_resamp_", v, "_run_", u, ".txt --skip ", 
		  burnin, " --outgroup 1 --outfile GALAX/OUTPUT/", foldername, "/", "galax_", paste0(names(roots), collapse = "_"), "_resamp_", v, "_run_", u))
		
		}
	  
	  }
	
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
