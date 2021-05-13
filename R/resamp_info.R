#' Import and summarize results from analyses set by the resamp_galax function.
#'
#' @param roots character, a named vector with the ontology IDs of the anatomical data subsets to produce MrBayes NEXUS files. 
#'   IDs and names should match those in the second and third columns of the ids data.frame, respectively.
#' @param nsamp integer, sets the number of samples for the resampling analysis.
#' @param runs integer, sets the number of independent runs of mcmc in the MrBayes commands block. Default is set to 2.
#' @param foldername character name of the output folder with all Galax output files.
#'
#' @return A list of matrices with values posterior coverage, Bayesian phylogenetic information (BPI) content, and phylogenetic dissonance
#'   among mcmc runs and among data subsets (if roots include more than one data subset) for each data subset specified in roots. 
#' 
#' @export
resamp_info <- function(roots = roots, nsamp = 100, runs = 2, foldername = "resamp")
{

  # Create a vector for the subfolder paths #
  x <- paste0("./GALAX/OUTPUT/", foldername, "/", "galax_")
  
  # Create a list to store results #
  res <- list()
  
  # Among-runs summary #
  for(i in 1:length(roots)){
  
      # Create a vector to store coverage, information and dissonance values #
	  w <- numeric()
	  
	  for(j in 1:nsamp){
	  
	    # Read all .txt files and extract coverage, information and dissonance values #
		z <- readLines(paste0(x, names(roots)[i], "_resamp_", j, ".txt"))
		z <- strsplit(z[grep(z, pattern = "merged")[3]], split = " ")
		suppressWarnings(z <- as.numeric(z[[1]]))
		z <- z[!is.na(z)][c(2,6,8)]
		
		# Build a raw matrix with values from all simulations #
		w <- cbind(w, z)
	  
	  }
	  
	  # Name all rows of each matrix #
	  rownames(w) <- c("Cover", "Info", "Disso")
	  colnames(w) <- paste0("Simu_", 1:nsamp)
	  
	  # Store each matrix in results list #
	  res[[i]] <- w
  
  }
  
  # Name list elements #
  names(res) <- names(roots)
  
  # Among-partitions summary #
  if(length(roots) > 1){
  
    # Create a list to store results #
	res2 <- list()
	
	for(k in 1:nsamp){
	
	  # Create a vector to store coverage, information and dissonance values #
	  s <- numeric()
	  
	  for(u in 1:runs){
	  
	    # Read all .txt files and extract coverage, information and dissonance values #
		v <- readLines(paste0(x, paste0(names(roots), collapse = "_"), "_resamp_", k, "_run_", u, ".txt"))
		v <- strsplit(v[grep(v, pattern = "merged")[3]], split = " ")
		suppressWarnings(v <- as.numeric(v[[1]]))
		v <- v[!is.na(v)][c(2,6,8)]
		
		# Build a raw matrix with values from all simulations #
		s <- cbind(s, v)
	  
	  }
	  
	  # Summarize results from a given simulation across all runs #
	  s <- cbind(s, apply(s[,1:runs], 1, mean))
	  
	  # Name all rows of each matrix #
	  rownames(s) <- c("Cover", "Info", "Disso")
	  colnames(s) <- c(paste0("R", 1:runs), "Mean")
	  
	  # Store summary for each simulation #
	  res2[[k]] <- s
	
	}
	
	# Name list elements #
	names(res2) <- paste0("Simu_", 1:nsamp)
  
  }
  
  # Return results #
  if(length(roots) > 1){ return(list(RUNS = res, PART = res2)) }else{ return(RUNS = res) }
  
}
