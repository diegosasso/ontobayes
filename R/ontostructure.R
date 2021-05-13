#' Gets the ontology hierarchical structure represented as a semantic similarity or dissonance dendrogram #
#'
#' @param ids data.frame, if manual = F, ids should be a data.frame providing matches between characters and anatomical entities.
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#'   Otherwise, ids can be ommited.
#' @param dsm matrix, if manual = T, dsm should be a dissimilarity matrix obtained from the function bayesinfo_pairwise.
#'   Otherwise, dsm can be ommited.
#' @param similarity character, the metric chosen to produce a semantic similarity dendrogram. Default is set to "jaccard".
#'   This option is the only one currently available. Future developments will allow to use other metrics such as "resnik".
#'   This options is available only if Uberon ontology terms are employed since it uses dependencies from rphenoscape to build
#'   a semantic similarity dendrogram as a representation of the ontology structure.
#' @param manual logical, if TRUE specify that a dissimilarity matrix based on phylogenetic dissonance will be provide to produce a dendrogram.
#'   Otherwise, a dendrogram will be produced using a dissimilarity matrix based on the semantic similarity metric chosen.
#' @param plot logical, If TRUE plots the resultant semantic similarity or dissonance dendrogram.
#'
#' @return A list with all terms subtended by each internal node of the semantic similarity or dissonance dendrogram.
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
#' }
#'
#' @export
ontostructure <- function(ids = ids, dsm = dsm, similarity = "jaccard", manual = F, plot = F)
{
    
   if(manual == T){
     
     # Get hierarchical cluster using a provided dissimilarity matrix #
     htree <- hclust(as.dist(dsm))
     
   }else{
     
     # Organize vector of IRIs #
     x <- paste0("http://purl.obolibrary.org/obo/", gsub(unique(ids[,2]), pattern = ":", replacement = "_"))
     
     if(similarity == "jaccard"){
       
       # Get Jaccard semantic similarity matrix for all terms #
       sm <- rphenoscape:::jaccard_similarity(terms = x, .colnames = c("label"))
       
     }
     
     # Get hierarchical cluster using ss #
     htree <- hclust(as.dist(1 - sm))
     
   }
   
   if(plot == T){
     
     # Plot hierarchical cluster dendrogram #
     plot(htree, cex = 0.4)
     
   }
   
   # Convert dendrogram to a phylo object #
   phy <- ape:::as.phylo(htree)
   
   # Organize tip labels #
   phy$tip.label <- gsub(phy$tip.label, pattern = " ", replacement = "_")
   phy$tip.label <- gsub(phy$tip.label, pattern = "/", replacement = "-")
   
   # Extract all subtrees #
   phy <- ape:::subtrees(phy)
   
   # Get all terms in each subtree #
   phy <- lapply(phy, function(x) x$tip.label)
   
   # Name all nodes #
   names(phy) <- paste0("N", (length(phy) + 2): (length(phy) * 2 + 1))
   
   # Return results #
   return(phy)
   
}
