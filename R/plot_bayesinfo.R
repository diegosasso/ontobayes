#' Plot Bayesian phylogenetic information or dissonance values at tips or internal nodes of a reference semantic similarity or dissonance dendrogram.
#'
#' @param phy phylo or hclust, a reference semantic similarity or dissonance dendrogram.
#' @param bpi matrix or list, information teory metrics to be ploted onto the reference dendrogram. If mode = "terminal", 
#'   a matrix with values of posterior coverage, Bayesian phylogenetic information, and phylogenetic dissonance as the output of the bayesinfo function.
#'   If mode = "internal", a list of matrices with values of posterior coverage, Bayesian phylogenetic information, and phylogenetic dissonance for each 
#'   internal node of the semantic similarity or dissonance dendrogram as the output of the bayesinfo_structure function.
#' @param struc list, the output from the ontostructure function.
#' @param stdz numeric, amount to standardize values Bayesian phylogenetic information and dissonance. If mode = "terminal", stdz should be a single number. 
#'   If mode = "internal", stdz should be a pair of numbers. Default is set to c(1,1).
#' @param plot logical, if TRUE plots values of Bayesian phylogenetic information or dissonance at internal nodes when mode = "internal".
#'   Otherwise, only returns a matrix with mean values between mcmc runs of Bayesian phylogenetic information and dissonance for all internal nodes of the reference dendrogram.
#' @param tcex numeric, sets the font size for tip labels in the dendrogram. Default is set to 0.5.
#' @param type character, sets the type of information theory metrics to plot. If type = "information", plots Bayesian phylogenetic information.
#'   If type = "dissonance", plots phylogenetic dissonance. Defaut is set to "information".
#' @param mode character, sets whether information theory metrics should be plotted for terminals or internal nodes of the reference dendrogram.
#'   If mode = "terminal", plots the chosen IT metric for individual anatomical data subsets defined by the ontology terms at tips.
#'   If mode = "internal", plots the chosen IT metric for comparisons among anatomical data subsets subtended by the internal nodes indicated in the reference dendrogram. 
#'
#' @return A matrix with mean values between mcmc runs of Bayesian phylogenetic information and dissonance for all internal nodes of the reference dendrogram if mode = "internal", 
#'   otherwise only plots the dendrogram with the chosen IT metrics (i.e., information or dissonance) and mode (i.e., terminal or internal). 
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
#' # Get a reference dendrogram based on phylogenetic dissonance #
#' phy1 <- hclust(as.dist(pp))
#'
#' # Get hierarchical structure based on dissonance dendrogram #
#' hstruc1 <- ontostructure(dsm = pp, manual = T, plot = F)
#'
#' # Import results from comparisons among terms according to dissonance dendrogram #
#' struc1 <- bayesinfo_structure(struc = hstruc1, runs = 2)
#'
#' # Import and summarize results from comparisons among runs for all terms #
#' runs1 <- bayesinfo(ids = ID) 
#'
#' # Plot dendrogram with information of individual anatomical data subsets at tips # 
#' plot_bayesinfo(phy = phy1, bpi = runs1, stdz = 1, mode = "terminal")
#'
#' # Plot dendrogram with information of among anatomical data subsets at internal nodes # 
#' plot_bayesinfo(phy = phy1, struc = hstruc1, bpi = struc1, stdz = c(10, 10), mode = "internal", type = "information")
#'
#' # Plot dendrogram with dissonance of among anatomical data subsets at internal nodes #
#' plot_bayesinfo(phy = phy1, struc = hstruc1, bpi = struc1, stdz = c(10, 10), mode = "internal", type = "dissonance")
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
#'
#' # Build a semantic similarity dendrogram"] #
#' sm <- paste0("http://purl.obolibrary.org/obo/", gsub(unique(ID[,2]), pattern = ":", replacement = "_"))
#' sm <- rphenoscape:::jaccard_similarity(terms = sm, .colnames = c("label"))
#' phy2 <- hclust(as.dist(1 - sm))
#'
#' # Import and summarize results from comparisons among runs for all terms #
#' runs2 <- bayesinfo(ids = ID) 
#'
#' # Plot dendrogram with information of individual anatomical data subsets at tips # 
#' plot_bayesinfo(phy = phy2, bpi = runs2, stdz = 1, mode = "terminal")
#'
#' # Plot dendrogram with information of among anatomical data subsets at internal nodes # 
#' plot_bayesinfo(phy = phy2, struc = hstruc2, bpi = struc2, stdz = c(10, 10), mode = "internal", type = "information")
#'
#' # Plot dendrogram with dissonance of among anatomical data subsets at internal nodes #
#' plot_bayesinfo(phy = phy2, struc = hstruc2, bpi = struc2, stdz = c(10, 10), mode = "internal", type = "dissonance")
#' }
#'
#' @export
plot_bayesinfo <- function(phy = phy, struc = struc, bpi = bpi, stdz = c(1, 1), 
                           plot = T, tcex = 0.5, type = "information", mode = "internal")
{
  
  # Import term tree #
  if(class(phy) == "hclust"){ phy <- ape:::as.phylo(phy) }
  
  # Reorganize tip labels #
  phy$tip.label <- gsub(phy$tip.label, pattern = " ", replacement = "_")
  phy$tip.label <- gsub(phy$tip.label, pattern = "/", replacement = "-")
  
  if(mode == "internal"){
  
  # Create some vectors to store values #
  nodes <- numeric()
  info <- numeric()
  diss <- numeric()
  
  # Extract all nodes corresponding to term clusters #
  for(i in 1:length(struc)){
    
	# Get nodes #
	nodes <- c(nodes, ape:::getMRCA(phy, struc[[i]]))
	
    # Get all clusters at a given level #
    groups <- struc[[i]]
    
	# Get information #
	info <- c(info, bpi[[i]]["Mean", 2])
	
	# Get dissonance #
	diss <- c(diss, bpi[[i]]["Mean", 3])
    
  }
    
  # Concatenate and organize the vectors #
  res <- rbind(nodes, info, diss)
  res <- res[,order(res[1,])]
  colnames(res) <- res[1,]
  res <- res[-1,]
  
  if(plot == T){
    
    if(type == "information"){
      
      # Plot term tree #
      ape:::plot.phylo(ape:::as.phylo(phy), cex = tcex, label.offset = 0.001, main = "Information")
      
      #ape:::nodelabels(res[1,]), frame = "none", col = "blue", cex = 0.8)
      ape:::nodelabels(pch = 21, lwd = 0.75, cex = round((res[1,]/stdz[1]), 2), 
	  bg = scales:::alpha("blue", alpha = 0.3), col = NULL)
      
    }
    
    if(type == "dissonance"){
      
      # Plot term tree #
      ape:::plot.phylo(ape:::as.phylo(phy), cex = tcex, label.offset = 0.001, main = "Dissonance")    
      
      #ape:::nodelabels(res[2,]), frame = "none", col = "red", cex = 0.8)
	  ape:::nodelabels(pch = 21, lwd = 0.75, cex = round((res[2,]/stdz[2]), 2), 
	  bg = scales:::alpha("red", alpha = 0.3), col = NULL)
      
    }
    
  }
  
  # Return results #
  return(res)
  
  }
  
  if(mode == "terminal"){

    # Ensure tree is in an apropriate format #
    phy <- phytools::untangle(phy, method = "read.tree")
    
    # Reorganize tip labels #
    phy$tip.label <- gsub(phy$tip.label, pattern = " ", replacement = "_")
    phy$tip.label <- gsub(phy$tip.label, pattern = "/", replacement = "-")
	
	# Reorganize some labels #
	bpi <- bpi[match(phy$tip.label, rownames(bpi)), 2]
  
    # Plot term tree with barplot of informational content #
	phytools:::plotTree.barplot(phy, bpi/stdz, 
                            args.plotTree = list(fsize = tcex),
                            args.barplot = list(col = scales:::alpha("blue", alpha = 0.3), 
                                                border = NA, main = "Information"))
												
  }
  
}
