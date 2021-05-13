#' Import and summarize results from analyses set by the galaxscrip_profiles function.
#'
#' @param tags list, profiles with anatomy ontology terms to query. The list of profiles should be organized as indicated in the examples.
#' @param runs integer, the number of independent runs of mcmc used in the MrBayes commands block. Default is set to 2.
#'
#' @return A list with matrices. Each matrix shows values of posterior coverage, Bayesian phylogenetic information, and phylogenetic dissonance
#'   among mcmc runs for each subset of data defined by an anatomy ontology term and across all subsets for each profile defined. 
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
#' # Import and summarize results from profiles #
#' profs1 <- bayesinfo_profiles(tags = tags, runs = 2)
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
#' # Import and summarize results from profiles #
#' profs2 <- bayesinfo_profiles(tags = tags, runs = 2)
#' }
#'
#' @export
bayesinfo_profiles <- function(tags = tags, runs = 2)
{

  # Load all query terms #
  tags <- tags
    
  # Create a vector for all paths #
  x <- paste0("GALAX/OUTPUT/", 1:length(names(tags)), "_", names(tags))
    
  # Create a list to store results #
  res <- list()
    
  for(i in 1:length(tags)){
    
    # Create a vector for the subfolder .txt file names #
    y <- list.files(x[i])
    y <- y[grep(y, pattern = ".txt")]
      
    # Create a vector to store coverage, information and dissonance values #
    w <- numeric()
      
    for(j in 1:length(y)){
      
      # Read all .txt files and extract coverage, information and dissonance values #
      z <- readLines(paste0(x[i], "/", y[j]))
      z <- strsplit(z[grep(z, pattern = "merged")[3]], split = " ")
      suppressWarnings(z <- as.numeric(z[[1]]))
      z <- z[!is.na(z)][c(2,6,8)]  
        
      w <- rbind(w, z)
	  
    }
    
    # Build matrix with all information and dissonance values #
    w <- rbind(w, apply(w[(dim(w)[1]-1):dim(w)[1],], 2, mean))
    
    # Name rows and columns of each matrix #
    colnames(w) <- c("Cover", "Info", "Disso")
    rownames(w) <- c(names(tags[[i]]), paste0("Overall_R", 1:(length(y) - length(tags[[i]]))), "Mean")
    
    # Store each matrix in results list #
    res[[i]] <- w
    
  }
  
  # Name list elements using tags provided #
  names(res) <- names(tags)
  
  # Return results #
  return(res)  

}
