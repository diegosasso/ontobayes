#' A wraper function to visualize the posterior distribution of tree topologies for data subsets based on selected anatomy ontology terms. 
#' 
#' @description This wraper uses the functions makeplot.topology and makeplot.treespace from package rwty (Warren et al. 2017) to provide an easy visualization 
#'   of the tree space for the posterior distributions of data subsets defined by selected anatomy ontology terms.
#'
#' @param ids data.frame, matches between characters and anatomical entities. 
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#' @param tags character, a vector with the names of the anatomical data subsets to produce the tree spaces. 
#'   Names should match those in the third column of the ids data.frame.
#' @param burnin integer, the absolute number of trees as the burn-in used in the mcmc. Default is set to 251.
#' @param res integer, if mode = "densitymap", sets the number of trees to sample from the posterior distribution.
#' @param tdist character, sets the type of tree distance metric to use. The options are "PD" for path distance or "RB" for Robinson Foulds distance.
#'   Default is set to "PD". See the help of makeplot.topology function in the rwty package for more details.
#' @param mode character, sets the type of visualization to use. The options are "distribution" to use the tree distance trace plot 
#'   from the makeplot.treespace function of rwty or "densitymap" to use the density map plot from the makeplot.treespace function of rwty. 
#'   See the help of makeplot.topology  and makeplot.treespace functions in the rwty package for more details.
#'
#' @return A tree trace plot or a density map plot for data subsets based on selected anatomy ontology terms. 
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
#' # Choose some terms annotated to characters in the matrix #
#' tg1 <- c("mandible", "maxilla")
#'
#' # Plot tree trace plot from RWTY package #
#' plot_treedistr(ids = ID, tags = tg1, burnin = 251, res = 500, tdist = "PD", mode = "distribution")
#'
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # Choose some terms annotated to characters in the matrix #
#' tg2 <- c("premaxilla", "maxilla")
#'
#' # Plot tree trace plot from RWTY package #
#' plot_treedistr(ids = ID, tags = tg2, burnin = 251, res = 500, tdist = "PD", mode = "distribution")
#' }
#'
#' @export
plot_treedistr <- function(ids = ids, tags = tags, burnin = 251, res = 100, tdist = "PD", mode = "distribution")
{
  
  # Import all terms #
  x <- unique(ids[,3])
  
  # Reorganize term names #
  x <- gsub(x, pattern = " ", replacement = "_")
  x <- gsub(x, pattern = "/", replacement = "-")
  
  # Reorganize tag names #
  tags <- gsub(tags, pattern = " ", replacement = "_")
  tags <- gsub(tags, pattern = "/", replacement = "-")
  
  # Get path to all tree folders #
  folders <- paste0("./NEXUS/ALL_PART/", 1:length(x), "_", x)
  
  # Name all paths #
  names(folders) <- x
  
  # Get correct paths according to tags #
  y <- folders[names(folders) %in% tags]
  
  # Create a list to store all trees #
  trees <- list()
  
  for(i in 1:length(y)){
    
    # Import all trees from MrBayes analyses #
    trees <- c(trees, rwty:::load.multi(path = y[i]))
    
  }
  
  # Reorganize tree names #
  names(trees) <- gsub(names(trees), pattern = "\\.nex\\.", replacement = "_")
  names(trees) <- gsub(names(trees), pattern = "\\.t", replacement = "")
  
  if(mode == "distribution"){
    
    # Make a distribution of tree topologies based on tree distances #
    topo <- suppressWarnings(rwty:::makeplot.topology(trees, burnin = burnin, independent.chains = F, 
                                                      treedist = tdist, approx.ess = F))
    
    # Plot distribution #
    return(topo$density.plot)
    
  }
  
  if(mode == "densitymap"){
    
    # Make a densitymap of tree topologies based on tree distances #
    space <- suppressWarnings(rwty:::makeplot.treespace(trees, burnin = burnin, 
                                                        n.points = res, fill.color = "LnL"))
    
    # Plot densitymaps #
    return(space$treespace.heatmap)
    
  }
  
}
