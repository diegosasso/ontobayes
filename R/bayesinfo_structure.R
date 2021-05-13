#' Import and summarize results from analyses set by the galaxscrip_structure function.
#'
#' @param struc list, the output from the ontostructure function.
#' @param runs integer, the number of independent runs of mcmc used in the MrBayes commands block. Default is set to 2.
#'
#' @return A list of matrices with values of posterior coverage, Bayesian phylogenetic information, and phylogenetic dissonance
#'   for each internal node of the semantic similarity or dissonance dendrogram.
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
#' # Import and visualize results from pairwise comparisons of terms #
#' pp <- bayesinfo_pairwise(ids = ID, runs = 2, profile = F, cluster = T)
#'
#' # Get hierarchical structure based on dissonance dendrogram #
#' hstruc1 <- ontostructure(dsm = pp, manual = T, plot = F)
#'
#' # Import results from comparisons among terms according to dissonance dendrogram #
#' struc1 <- bayesinfo_structure(struc = hstruc1, runs = 2)
#'
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # Get hierarchical structure based on semantic similarity dendrogram #
#' hstruc2 <- ontostructure(ids = ID, manual = F, plot = F, similarity = "jaccard")
#'
#' # Import results from comparisons among terms according to ontology structure #
#' struc2 <- bayesinfo_structure(struc = hstruc2, runs = 2)
#' }
#'
#' @export
bayesinfo_structure <- function(struc = struc, runs = 2)
{

  # Create a list to store results #
  res <- list()
  
  for(i in 1:length(struc)){
  
    # Create a vector to store coverage, information and dissonance values #
	w <- numeric()
	
	for(j in 1:runs){
	
	  # Read all .txt files and extract coverage, information and dissonance values #
	  z <- suppressWarnings(readLines(paste0("./GALAX/OUTPUT/structure/", "galax_", names(struc)[i], "_run", j,".txt")))
	  z <- strsplit(z[grep(z, pattern = "merged")[3]], split = " ")
	  suppressWarnings(z <- as.numeric(z[[1]]))
	  z <- z[!is.na(z)][c(2,6,8)]
	  
	  # Store values #
	  w <- rbind(w, z)
	
	}
	
	# Summarize values among runs for a given group of terms #
	w <- rbind(w, apply(w[1:runs,], 2, mean))
	
	# Name all columns #
	colnames(w) <- c("Cover", "Info", "Disso")
	rownames(w) <- c(paste0("R", 1:runs), "Mean")
	
	# Store summaries #
	res[[i]] <- w
  
  }
  
  # Name all node referring to subtrees #
  names(res) <- names(struc)
  
  # Return results #
  return(res)

}
