#' Get node specific componenets of Bayesian phylogenetic information, posterior probability, and internode confidence for all majority rule trees 
#'   obtained from the analyses of data subsets defined by anatomy ontology terms.
#'
#' @param phylo, a reference phylogenetic species tree.
#' @param foldername character, the name of the Galax output folder with the results from ontobayes analysis.
#'
#' @return A list of matrices with values of Bayesian phylogenetic information, posterior probability, and internode confidence  
#'   for each node in the reference phylogenetic species tree.
#'
#' @examples
#' \dontrun{
#' # Example: HAO or UBERON #
#' # Import a reference phylogenetic species tree #
#' reftree <- ape:::read.nexus(file = "./data/reftree.tre")
#'
#' # Get node-specific components from each tree inferred from all different data subsets defined in ontobayes # 
#' M <- get_bayesinfoHM(phy = reftree, foldername = "all_partitions")
#' }
#'
#' @export
get_bayesinfoHM <- function(phy = phy, foldername = "all_partitions")
{

	# Organize tip labels of species tree #
	phy$tip.label <- gsub(phy$tip.label, pattern = " ", replacement = "_")
	phy$tip.label <- gsub(phy$tip.label, pattern = "/", replacement = "-")
	
	# Get path to MJ trees inside the Galax folder #
	path <- paste0("./GALAX/OUTPUT/", foldername,"/")
	
	# Get all file names #
	files <- grep(list.files(path), pattern = "merged.tre", value = T)
	
	# Organize tags for partition names #
	tags <- gsub(files, pattern = "galax_", replacement = "")
	tags <- gsub(tags,  pattern = "-merged.tre", replacement = "")
	
	# Create list to store trees #
	trees <- list()
	
	# Import MJ trees from all partitions #
	for(i in 1:length(files)){
	
	  # Read tree nexus files (as lines due to some issues with Galax tree nexus files) #
	  nexus <- readLines(paste0(path, files[i]))
	  
	  # Get taxon labels #
	  taxa <- nexus[(grep(nexus, pattern = "translate") + 2):(grep(nexus, pattern = "tree majrule = ") - 1)]
	  
	  # Organize taxon labels #
	  taxa <- gsub(taxa, pattern = "\\s", replacement = "")
	  taxa <- gsub(taxa, pattern = "\\d{1,}", replacement = "")
	  taxa <- gsub(taxa, pattern = "'", replacement = "")
	  taxa <- gsub(taxa, pattern = ",", replacement = "")
	  taxa <- gsub(taxa, pattern = ";", replacement = "")
	  
	  # Extract the actual MJ tree #
	  tree <- grep(nexus, pattern = "tree majrule", value = T)
	  tree <- gsub(tree, pattern = "tree majrule = ", replacement = "")
	  
	  # Replace taxon labels #
	  for(j in 1:length(taxa)){
	  
	    tree <- gsub(tree, pattern = paste0(j, ":"), replacement = paste0(taxa[j], ":"))
	  
	  }
	  
	  # Store trees #
	  trees[[i]] <- ape:::read.tree(text = tree)
	
	}
	
	# Name all trees #
	names(trees) <- tags
	
	# Create some variables to store values #
	PP <- numeric()
	IF <- numeric()
	IC <- numeric()
	
	# Extract values from all nodes of all partition trees #
	x <- lapply(lapply(trees, function(x) x$node.label), gsub, pattern = "\"", replacement = "")
	x[lapply(x, length) < 1] <- NA
	
	for(k in 1:length(x)){
	
	  # Get all values from all nodes of a given partition tree #
	  y <- x[[k]]
	  
	  # Extract desired values from nodes of a given partition tree #
	  pp <- stringr::str_extract(y, pattern = "P=\\d{1,}\\.\\d{1,}")
	  pp <- as.numeric(gsub(pp, pattern = "P=", replacement = ""))
	  pp <- round(pp, 2)
	  
	  ii <- stringr::str_extract(y, pattern = "I=\\d{1,}\\.\\d{1,}")
	  ii <- as.numeric(gsub(ii, pattern = "I=", replacement = ""))
	  ii <- round(ii, 2)
	  
	  ic <- stringr::str_extract(y, pattern = "IC=\\d{1,}\\.\\d{1,}")
	  ic <- as.numeric(gsub(ic, pattern = "IC=", replacement = ""))
	  ic <- round(ic, 2)
	  
	  # Match nodes from reference tree to a given partition tree #
	  nodes.pp <- phytools:::matchNodes(phy, trees[[k]])[,2]
	  nodes.if <- phytools:::matchNodes(phy, trees[[k]])[,2]
	  nodes.ic <- phytools:::matchNodes(phy, trees[[k]])[,2]
	  
	  # Reorder and organize values #
	  nodes.pp[order(nodes.pp)][1:length(pp)] <- pp
	  nodes.if[order(nodes.if)][1:length(ii)] <- ii
	  nodes.ic[order(nodes.ic)][1:length(ic)] <- ic
	  
	  # Build a matrix with all values from all partitions #
	  PP <- rbind(PP, nodes.pp)
	  IF <- rbind(IF, nodes.if)
	  IC <- rbind(IC, nodes.ic)
	
	}
	
	# Name rows and columns #
	rownames(PP) <- tags
	rownames(IF) <- tags
	rownames(IC) <- tags
	colnames(PP) <- paste0("N", 1:dim(PP)[2])
	colnames(IF) <- paste0("N", 1:dim(IF)[2])
	colnames(IC) <- paste0("N", 1:dim(IC)[2])
	
	# Return results #
	return(list(PP = PP, IF = IF, IC = IC))

}
