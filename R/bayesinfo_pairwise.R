#' Import and summarize results from analyses set by the galaxscrip_pairwise function.
#'
#' @param ids data.frame or list, if profile = F, ids should be a data.frame providing matches between characters and anatomical entities. 
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#'   Otherwise, ids should be a list of profiles organized as indicated in the examples.
#' @param profile logical, if TRUE specify that the input in ids is a list of profiles to query, organized as indicated in the examples.
#' @param cluster logical. If TRUE plots the resultant dissonance clustering dendrogram.
#'
#' @return A dissimilarity matrix based on values of phylogenetic dissonance obtained from all pairwise comparisons 
#'   between anatomical data subsets.
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
#' # Import and visualize results from all pairwise comparisons of terms #
#' pp1 <- bayesinfo_pairwise(ids = ID, runs = 2, profile = F, cluster = T)
#' 
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
#' # Import and visualize results from all pairwise comparisons of terms in user defined profiles #
#' tt1 <- bayesinfo_pairwise(ids = tags, runs = 2, profile = T, cluster = T)
#'
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # Import and visualize results from all pairwise comparisons of terms #
#' pp2 <- bayesinfo_pairwise(ids = ID, runs = 2, profile = F, cluster = T)
#' 
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
#' # Import and visualize results from all pairwise comparisons of terms in user defined profiles #
#' tt2 <- bayesinfo_pairwise(ids = tags, runs = 2, profile = T, cluster = T)
#' }
#'
#' @export
bayesinfo_pairwise <- function(ids = ids, runs = 2, profile = F, cluster = T)
{
  
  if(profile == T){
	
    # Load all query terms #
    tags <- ids
	
	# Create an empty vector #
    x <- character()
    
    for(i in 1:length(tags)){
      
      # Create vector for the file names #
      x <- c(x, paste0(1:length(tags[[i]]), "_", names(tags[[i]])))
      
    }
    
  }else{
    
    # Create vector for the file names #
    x <- gsub(ids[,3], pattern = " ", replacement = "_")
    x <- gsub(x, pattern = "/", replacement = "-")
    x <- unique(x)
	
  }    
  
  # Create a vector for all possible pairwise comparisons between files #
  y <- combn(x, m = 2, simplify = F)
  y <- lapply(y, paste0, sep = "", collapse = "_")
  
  # Create matrix to store dissonance values #
  mat <- numeric()
  
  for(k in 1:runs){
    
    # Create vectors to store dissonance values #
    r <- numeric()
    
    for(j in 1:length(y)){
      
      # Read all .txt files and extract dissonance values #
      txt <- readLines(paste0("GALAX/OUTPUT/pairwise/galax_", y[[j]], "_run", k,".txt"))
      txt <- strsplit(txt[grep(txt, pattern = "merged")[3]], split = " ")
      
      # Extract dissonance values and convert to numeric #
      txt <- suppressWarnings(as.numeric(txt[[1]]))
      txt <- txt[!is.na(txt)][8]
      
      # Store dissonance values into vectors #
      r[j] <- txt
      
    }
    
    # Combine dissonance values from all runs #
    mat <- rbind(mat, r)
    
  }
  
  # Create a matrix to store results #
  z <- matrix(NA, length(x), length(x))
  
  # Store mean dissonance values among runs into the matrix #
  z[upper.tri(z, diag = F)] <- apply(mat, 2, mean)
  
  if(profile == T){
  
    # Name rows and columns of the matrix #
	colnames(z) <- unlist(lapply(tags, names))
	rownames(z) <- unlist(lapply(tags, names))
  
  }else{
  
    # Name rows and columns of the matrix #
	colnames(z) <- x
	rownames(z) <- x
  
  }
  
  # Standardize dissonance values #
  clst <- z/100
  
  # Set diagonal values to 0 #
  diag(clst) <- 0
  
  # Transpose matrix #
  clst <- t(clst)
  
  if(cluster == T){
    
    # Plot hierarquical cluster #
    plot(hclust(as.dist(clst)))
    
  }
  
  # Return results #
  return(clst)
  
}
