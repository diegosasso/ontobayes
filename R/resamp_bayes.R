#' Builds MrBayes NEXUS files for different types of resampling analyses.
#'
#' @param ids data.frame, matches between characters and anatomical entities. The first column should provide character ids in the same order presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#' @param type character, sets the type of resampling analysis. Options are "none" to resample the same dataset defined in dataset; 
#'   "char" to resample characters from the dataset defined in dataset; "tax" to resample taxa from the dataset defined in dataset;
#'   "both" to resample characters and taxa from the dataset defined in dataset.
#' @param nsamp integer, sets the number of samples for the resampling analysis.
#' @param nchars integer, sets the number of characters to resample for each data subset if mode = "char" or "both". If more than one data subset
#'   is selected then should be a vector with more than one element.
#' @param ntax integer, sets the number of taxa to resample if mode = "tax" or "both".
#' @param roots character, a named vector with the ontology IDs of the anatomical data subsets to produce MrBayes NEXUS files. 
#'   IDs and names should match those in the second and third columns of the ids data.frame, respectively.
#' @param manual logical, if TRUE specify that an ontology file (OBO format) is provided by the user. This option is the only one currently available. Future developments will allow to access Uberon ontology directly using rphenoscape dependencies.
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
#' @param foldername character name of the output folder to save all MrBayes NEXUS files.
#'
#' @return nsamp NEXUS files for each data subset specified in roots.
#'   The files will be organized in the foldername subfolder inside the folder NEXUS/RESAMP.
#'
#' @export
resamp_bayes <- function(ids = ids, type = "none", nsamp = 100, nchars = c(10, 10), ntax = 10, roots = roots,
manual = T, ontology = "data/ontology.obo", relations = c("part_of", "is_a"), 
dataset = "data/matrix.nex", outgroup = "outgroup_taxon", rates = "gamma", coding = "all", ratepr = "fixed", brlenspr = "Unconstrained:Exp(10)", 
gen = 1000000, runs = 2, chains = 4, freq = 1000, foldername = "resamp")
{

  # Create MrBayes folders #
  suppressWarnings(dir.create(file.path("NEXUS", "RESAMP", foldername), recursive = T))

  # Import matrix #
  suppressWarnings(nex <- readLines(con = dataset))
  
  # Copy the original matrix #
  mat <- nex

  # Extract total number of characters in the matrix #
  z <- stringr::str_extract(stringr::str_extract(grep(x = nex, pattern = "NCHAR", value = T), pattern = "(NCHAR)\\s*(=)\\s*(\\d{1,})(;)"), pattern = "\\d{1,}")
  z <- paste(1:as.numeric(z))

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

  if(manual == T){
    
    # Load ontology #
    ont <- ontologyIndex::get_OBO(file = ontology, extract_tags = "everything", propagate_relationships = relations)
    
	# Adjust names of ontology identifiers #
	ids[,2] <- gsub(ids[,2], pattern = "_", replacement = ":")
	
	# Import terms to query #
	roots <- gsub(roots, pattern = "_", replacement = ":")
	roots <- as.list(roots)
	
  }else{
  
	# Adjust names of ontology identifiers #
	ids[,2] <- gsub(ids[,2], pattern = ":", replacement = "_")
	
	# Import terms to query #
	roots <- gsub(roots, pattern = ":", replacement = "_")
	roots <- paste0("http://purl.obolibrary.org/obo/", roots)
	ids[,2] <- paste0("http://purl.obolibrary.org/obo/", ids[,2])    
    roots <- as.list(roots)
	
  }

  # Change character IDs (characters must be ordered in the dataset!) #
  ids[,1] <- 1:length(ids[,1])
  
  if(type == "none"| type == "tax"){
    
	if(manual == T){
	
	  for(i in 1:length(roots)){
	  
	    # Store characters #
		roots[[i]] <- ids[,1][ids[,2] %in% ontologyIndex::get_descendants(ontology = ont, roots = roots[[i]])]
	  
	  }
	
	}else{
	
	  for(i in 1:length(roots)){
	  
	    # Store characters #
		roots[[i]] <- ids[,1][rphenoscape:::pk_is_descendant(term = roots[[i]], candidates = ids[,2], includeRels = "part_of")]
	  
	  }	
	
	}
    
  }
  
  for(k in 1:nsamp){
    
	if(type == "tax" | type == "both"){
	
	  # Adjust number of taxa in the matrix #
	  nex[grep(x = nex, pattern = "NTAX")] <- gsub(grep(x = nex, pattern = "NTAX", value = T), pattern = "(NTAX)\\s*(=)\\s*(\\d{1,})", replacement = paste0("NTAX = ", ntax))
	  
	  # Copy the original matrix #
	  mat <- nex
	  
	  # Get a random set of taxa #
	  w <- (grep(x = mat, pattern = "MATRIX", value = F) + 1):(grep(x = mat, pattern = "END", value = F) - 2)
	  w <- w[grep(x = mat[w], pattern = ".")]
	  w <- w[!w %in% grep(x = mat, pattern = outgroup)]
	  w <- w[!w %in% sample(w, (ntax - 1))]
	  mat <- mat[-w]
	
	}
	
	# Get a copy of all characters #
	x <- z
	
	for(j in 1:length(roots)){
	  
	  if(type == "char" | type == "both"){
	  
	    # Get mutually exclusive subsamples of characters #
		y <- sample(x, size = nchars[j])
		chars <- paste0(setdiff(z, y), sep = " ", collapse = "")
		x <- x[!x %in% y]
	  
	  }else{
	  
	    chars <- paste0(setdiff(z, roots[[j]]), sep = " ", collapse = "")
	  
	  }
	  
	  # Adjust file names #
	  bayes <- gsub(pattern = "&&exclude", replace = chars, x = block)
	  bayes <- gsub(bayes, pattern = "&&logfile", replacement = paste0(names(roots)[j], "_resamp_", k))
	  bayes <- gsub(bayes, pattern = "&&sumt", replacement = paste0(names(roots)[j], "_resamp_", k))
	  bayes <- gsub(bayes, pattern = "&&sump", replacement = paste0(names(roots)[j], "_resamp_", k))
	  
      # Write MrBayes files #
	  write(c(mat, bayes), file = paste0("./NEXUS/RESAMP/", foldername, "/", names(roots)[j], "_resamp_", k, ".nex"))
	  
	}
  
  }

  # MrBayes batch file #
  batch <- paste(
  "#NEXUS",
  "",
  "BEGIN MRBAYES;",
  "",
  "set autoclose = yes nowarn = yes;",
  "",
  sep = "\n")
  
  # Get file names #
  fnames <- character()

  for(u in 1:length(roots)){
  
    fnames <- c(fnames, paste0("execute ", names(roots)[u], "_resamp_", 1:nsamp, ".nex;"))
  
  }
  
  # Write MrBayes batch file #
  write(c(batch, paste(fnames, sep = "\n"), " ", "END;"), file = paste0("./NEXUS/RESAMP/", foldername, "/", "batch.nex"))

}
