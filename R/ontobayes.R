#' Builds MrBayes NEXUS files for subsets of characters annotated with ontology anatomy terms. 
#'
#' @param ids data.frame, matches between characters and anatomical entities.
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
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
#'   The files will be organized in individual subfolders inside the folder NEXUS/ALL_PART.
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
#' # Create MrBayes partition files #
#' ontobayes(ids = ID, dataset = "data/matrix1.nex", outgroup = "sp1", rates = "gamma", 
#' coding = "all", gen = 1000000, runs = 2, chains = 4)
#'
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # Create MrBayes partition files #
#' ontobayes(ids = ID, dataset = "data/matrix2.nex", outgroup = "sp1", rates = "gamma", 
#' coding = "all", gen = 1000000, runs = 2, chains = 4)
#' }
#'
#' @export
ontobayes <- function(ids = ids, dataset = "data/matrix.nex", 
                      outgroup = "outroup_taxon", rates = "gamma", coding = "all", ratepr = "fixed", brlenspr = "Unconstrained:Exp(10)", 
                      gen = 1000000, runs = 2, chains = 4, freq = 1000)
{
  
  # Import ontology identifiers and term names to query #
  ids <- ids
  
  # Change character IDs (characters must be ordered in the dataset!) #
  ids[,1] <- 1:length(ids[,1])
  
  # Adjust ontology term names #
  ids[,3] <- gsub(ids[,3], pattern = " ", replacement = "_")
  ids[,3] <- gsub(ids[,3], pattern = "/", replacement = "-")
  
  # Get all unique terms included in the dataset #
  tags <- unique(ids[,2])
  
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
  
  for(i in 1:length(tags)){
    
    # Create folders #
    dir.create(file.path("NEXUS", "ALL_PART", paste0(i, "_", unique(ids[,3])[i])), recursive = T)
    
    # Get all descendent terms from the query #  	  
    query <- ids[,1][ids[,2] %in% tags[i]]
    
    # Get partition name #
    part <- unique(ids[,3])[i]
    
    # Copy original MrBayes settings block #
    block2 <- block
    
    # Build MrBayes settings block for each partition #
    tag <- 0
    
    block2 <- gsub(pattern = "&&logfile.txt", replace = paste(part, ".txt", sep = ""), x = block2)
    block2 <- gsub(pattern = "&&sumt.nex", replace = paste(part, ".nex", sep = ""), x = block2)
    block2 <- gsub(pattern = "&&sump.nex", replace = paste(part, ".nex", sep = ""), x = block2)
    bayes <- gsub(pattern = "&&exclude", replace = paste0(setdiff(z, query), sep = " ", collapse = ""), x = block2)
    
    # Write partition files #
    write(c(nex, bayes), file = paste0(getwd(), "/NEXUS/ALL_PART/", i, "_", part, "/", part, ".nex"))
    
    # Create tag for partitions #
    tag <- paste("execute", " ", part, ".nex;", sep = "")
    
    # Write batch files #
    write(c(batch, paste(tag, sep = "\n"), " ", "END;"), file = paste0(getwd(), "/NEXUS/ALL_PART/", i, "_", part, "/", "batch.nex"))
    
  }
  
}
