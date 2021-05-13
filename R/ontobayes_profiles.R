#' Builds MrBayes NEXUS files for subsets of characters annotated with ontology anatomy terms based on a chosen profile of terms. 
#'
#' @param tags list, profiles with anatomy ontology terms to query. The list of profiles should be organized as indicated in the examples.
#' @param ids data.frame, matches between characters and anatomical entities. 
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#' @param manual logical, if TRUE specify that an ontology file (OBO format) is provided by the user. 
#'   This option is the only one currently available. Future developments will allow to access Uberon ontology directly using rphenoscape dependencies.
#' @param ontology character, name of the OBO file with the ontology.
#' @param relations character, set the types of ontological relations to use. Default is set to c("part_of", "is_a"). 
#' @param dataset character, name of the NEXUS file with the character matrix to be read. 
#' @param outgroup character, name of the outgroup taxon.
#' @param rates character, sets the model of among-character rate variation in the MrBayes commands block. Default is set to "equal". 
#' @param coding character, sets the coding bias in the MrBayes commands block. Default is set to "all".
#' @param ratepr character, sets the character specific or partition specific rates model in the MrBayes commands clock. Default is set to "fixed".
#' @param brlenspr character, sets the branch lengths prior in the MrBayes commands block. Deafault is set to "Unconstrained:Exp(10)".
#' @param gen integer, sets the number of generations of mcmc in the MrBayes commands block. Default is set to 1000000.
#' @param runs integer, sets the number of independent runs of mcmc in the MrBayes commands block. Default is set to 2.
#' @param chains integer, sets the number of chains for each run in the MrBayes commands block. Default is set to 4.
#' @param freq integer, sets the sampling frequency for the mcmc in the MrBayes commands block. Default is set to 1000.
#'
#' @return Individual NEXUS files for each subset of characters annotated with ontology anatomy terms.
#'   The files will be organized in individual subfolders named according to each profile inside the folder NEXUS.
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
#' # Create MrBayes files #
#' ontobayes_profiles(tags = tags, ids = ID, manual = T, relations = c("BFO:0000050", "is_a"), ontology = "data/HAO.obo", 
#' dataset = "data/matrix1.nex", outgroup = "sp1", rates = "gamma", coding = "variable", 
#' gen = 1000000, runs = 2, chains = 4)
#' 
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
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
#' # Create MrBayes files #
#' ontobayes_profiles(tags = tags, ids = ID, manual = T, relations = c("part_of", "is_a"), ontology = "data/UBERON.obo", 
#' dataset = "data/matrix2.nex", outgroup = "sp1", rates = "gamma", coding = "variable", 
#' gen = 1000000, runs = 2, chains = 4)
#' }
#' 
#' @export
ontobayes_profiles <- function(tags = tags, ids = ids, manual = T, relations = c("part_of", "is_a"), 
ontology = "data/ontology.obo", dataset = "data/matrix.nex", outgroup = "outroup_taxon", rates = "gamma", coding = "all", 
ratepr = "fixed", brlenspr = "Unconstrained:Exp(10)", gen = 1000000, runs = 2, chains = 4, freq = 1000)
{

  # Import ontology identifiers and profiles of terms to query #
  tags <- tags
  ids <- ids
  
  # Change character IDs (characters must be ordered in the dataset!) #
  ids[,1] <- 1:length(ids[,1])

  if(manual == T){
    
    # Load ontology #
    ont <- ontologyIndex::get_OBO(file = ontology, extract_tags = "everything", propagate_relationships = relations)
    
  }
  
  # Adjust names of ontology identifiers and terms to query #
  if(manual == T){
    
	ids[,2] <- gsub(ids[,2], pattern = "_", replacement = ":")
    tags <- lapply(tags, gsub, pattern = "_", replacement = ":")
	
  }else{
    
    ids[,2] <- gsub(ids[,2], pattern = ":", replacement = "_")
    ids[,2] <- paste0("http://purl.obolibrary.org/obo/", ids[,2])  
    
  }  
  
  # Import matrix #
  suppressWarnings(nex <- readLines(con = dataset))

  # Extract total number of characters in the matrix #
  z <- stringr::str_extract(stringr::str_extract(grep(x = nex, pattern = "NCHAR", value = T), pattern = "(NCHAR)\\s*(=)\\s*(\\d{1,})(;)"), pattern = "\\d{1,}")
  z <- paste(1:as.numeric(z))
  
  # Create MrBayes text blocks #
  # MrBayes settings block #
  block <- paste(
    "",
    "BEGIN MRBAYES;",
    "",
    paste0("outgroup ", outgroup, ";"),
    "log start filename = &&logfile.txt;",
    "",
    paste0("lset rates = ", rates, " coding = ", coding, ";"),
    "",
    paste0("charset IN = ", paste0(z, sep = " ", collapse = ""), ";"),
    "exclude &&exclude;",
    "partition FULL = 1: IN;",
    "set partition = FULL;",
    "",
	paste0("prset ratepr = ", ratepr," brlenspr = ", brlenspr, ";"),
    "",
    paste0("mcmc ngen = ", format(gen, scientific = F), " nruns = ", runs, 
       " nchains = ", chains, " printfreq = ", freq, " samplefreq = ", freq, 
       " diagnfreq = 1000 temp = 0.025;"),
    "",
    "sumt filename = &&sumt.nex;",
    "sump filename = &&sump.nex;",
    "",
    "END;",
    sep = "\n"
  )
  
  # MrBayes batch file #
  batch <- paste(
    "#NEXUS",
    "",
    "BEGIN MRBAYES;",
    "",
    "set autoclose = yes nowarn = yes;",
	"",
    sep = "\n")

  for(j in 1:length(tags)){
  
    # Create folders #
    dir.create(file.path("NEXUS", paste0(j, "_", names(tags)[j])), recursive = T)
    
	# Get all descendent terms from the query list of profiles #
    if(manual == T){
      
      # Query terms #
      query <- list()
      
      for(u in 1:length(tags[[j]])){
        
		query[[u]] <- ids[,1][ids[,2] %in% ontologyIndex::get_descendants(ontology = ont, roots = tags[[j]][u])]
		
      }
      
      names(query) <- names(tags[[j]])	  
      
    }else{
      
	  # Query terms #
      vv <- gsub(tags[[j]], pattern = ":", replacement = "_")
      vv <- paste0("http://purl.obolibrary.org/obo/", vv)
      
      query <- list()
      
      for(u in 1:length(tags[[j]])){
        
        query[[u]] <- ids[,1][rphenoscape:::pk_is_descendant(term = vv[u], candidates = ids[,2], includeRels = "part_of")]
        
      }
      
      names(query) <- names(tags[[j]])
      
    }
    
    # Copy original MrBayes settings block #
    block2 <- block
    
    # Build MrBayes settings block for each partition #
    tag <- 0
    for (v in 1:length(query)){
	
      block3 <- gsub(pattern = "&&logfile.txt", replace = paste(v, "_", names(query)[v], ".txt", sep = ""), x = block2)
      block3 <- gsub(pattern = "&&sumt.nex", replace = paste(v, "_", names(query)[v], ".nex", sep = ""), x = block3)
      block3 <- gsub(pattern = "&&sump.nex", replace = paste(v, "_", names(query)[v], ".nex", sep = ""), x = block3)
      bayes <- gsub(pattern = "&&exclude", replace = paste0(setdiff(z, query[[v]]), sep = " ", collapse = ""), x = block3)
      
      # Write partition files #
      write(c(nex, bayes), file = paste0(getwd(), "/NEXUS/", j, "_", names(tags)[j], "/", v, "_", names(query)[v], ".nex"))
      
      # Create tag for partitions #
      tag[v] <- paste("execute", " ", v, "_" , names(query)[v], ".nex;", sep = "")
	  
    }
    
    # Write batch file #
    write(c(batch, paste(tag, sep = "\n"), " ", "END;"), file = paste0(getwd(), "/NEXUS/", j, "_", names(tags)[j], "/", "batch.nex"))
    
  }

}
