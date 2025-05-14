#Load libraries
library(Biostrings)

#Set correct working directory (file output)
setwd("~/Desktop/...")

#Function to translate sequences that are (5'-->3')
translate_folder_5to3 <- function(input_dir, output_dir) {
  files <- list.files(input_dir, pattern = "profile", full.names = TRUE)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (file in files) {
    dna <- readDNAStringSet(file)
    protein <- translate(dna)
    
    output_name <- sub("profile_", "protein_", basename(file))
    
    output_file <- file.path(output_dir, output_name)
    writeXStringSet(protein, output_file)
  }
}

#Function to create the reverse complement and translate (for 3' -> 5')
translate_folder_3to5 <- function(input_dir, output_dir) {
  files <- list.files(input_dir, pattern = "profile", full.names = TRUE)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (file in files) {
    dna <- readDNAStringSet(file)
    dna_rc <- reverseComplement(dna)
    protein <- translate(dna_rc)
    
    output_name <- sub("profile_", "protein_", basename(file))
    
    output_file <- file.path(output_dir, output_name)
    writeXStringSet(protein, output_file)
  }
}


#Creating figures with saccharis 
#Package needed, one time installation. 
install.packages("ggnewscale", lib = .libPaths()[1])
BiocManager::install("ggtree", lib = .libPaths()[1], type = "source")
install.packages("rsaccharis_1.0.5.tar.gz", repos = NULL, type = "source", lib = .libPaths()[1])

library(rsaccharis)

#Set working directory(file output)
setwd("~/Desktop/...")

A_load_data()
B_plots_all()


## Script for translating protein sequence to nucletide sequence
library(Biostrings)

aa_to_codon <- list(
  A = "GCC", R = "CGC", N = "AAC", D = "GAC", C = "TGC",
  Q = "CAG", E = "GAG", G = "GGC", H = "CAC", I = "ATC",
  L = "CTG", K = "AAG", M = "ATG", F = "TTC", P = "CCC",
  S = "AGC", T = "ACC", W = "TGG", Y = "TAC", V = "GTC",
  "*" = "TAA", X = "NNN", B = "NNN", Z = "NNN", U = "TGA"  # ambiguous/selenocysteine
)

# Function to back-translate protein to nucleotide sequence
back_translate <- function(protein_seq, codon_map = aa_to_codon) {
  # Clean input: remove non-letter characters and uppercase
  clean_seq <- gsub("[^A-Za-z*]", "", protein_seq)
  clean_seq <- toupper(clean_seq)
  
  aa_vec <- strsplit(clean_seq, "")[[1]]
  
  codons <- sapply(aa_vec, function(aa) {
    codon <- codon_map[[aa]]
    if (is.null(codon)) {
      warning(paste("Unknown amino acid:", aa, "- using NNN"))
      codon <- "NNN"
    }
    return(codon)
  })
  
  dna_seq <- paste(codons, collapse = "")
  return(DNAString(dna_seq))
}

protein <- "WNQKWRGPLRDAMNKLKSIADEVLLREFPKISNVGPWEARNQYIQVLVKPEDVNRKNIFLQSILINPSEEKDRATAIRLLEIQKFCLFTFTSCGWFFNDIEGLEPVQNMRYAERAIELLKPFLPPSLPFKDQVLSILARARSNEHKLNGAEIFVQQAIPEVPVALKLMAEQAVIY"
dna <- back_translate(protein)
print(dna)
#as.character(dna)