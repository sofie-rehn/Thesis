##Split "main" fasta-file into smaller files (10mb) that can be uploaded to dbCAN
#Load libraries
library(Biostrings)

#Reading the large FASTA file
fasta <- readDNAStringSet("PSL-8W-2.fasta") #Update depending on file-name

#Setting threshold in number of bases (10 MB ~ 10e6 bases)
threshold <- 10e6

#Initialising variables to accumulate sequences and track sizes
current_group <- list()
current_size <- 0
group_counter <- 1

#Function to write out a group of contigs as a FASTA file
write_group <- function(group, counter) {
  group_set <- do.call(c, group)
  outfile <- paste0("split_", counter, ".fasta")
  writeXStringSet(group_set, filepath = outfile)
  message("Written file: ", outfile)
}

#Looping through each contig
for (i in seq_along(fasta)) {
  seq_length <- width(fasta[i])
  
  #Checking if adding this contig would exceed the threshold
  if (current_size + seq_length > threshold && length(current_group) > 0) {
    write_group(current_group, group_counter)
    group_counter <- group_counter + 1
    current_group <- list()
    current_size <- 0
  }
  
  #Adding the current contig
  current_group[[length(current_group) + 1]] <- fasta[i]
  current_size <- current_size + seq_length
}

#Writing any remaining contigs
if (length(current_group) > 0) {
  write_group(current_group, group_counter)
}


##Create new files for the interesting contigs
#splitting according to defined positions 

#Defining file paths
main_fasta <- "PSL-8W-1.fasta"  #Replace with  actual path
output_file <- file.path("contigs","all.fasta")  #Replace with  actual output directory

#Ensuring the output directory exists
if (!dir.exists(output_file)) {
  dir.create(output_file, recursive = TRUE)
}

#Reading the FASTA file
sequences <- readDNAStringSet(main_fasta)

#Defining the list of interesting contigs
interesting_contigs <- c("contig_7079",
                         "contig_790",
                         "contig_3933",
                         "contig_14795",
                         "contig_6989",
                         "contig_2711",
                         "contig_1721")  #Replace with actual names to extract

selected_sequences <- sequences[interesting_contigs]

#Checking if any contigs were found before writing
if (length(selected_sequences) > 0) {
  writeXStringSet(selected_sequences, output_file)
  cat("Saved combined FASTA file:", output_file, "\n")
} else {
  cat("Warning: No matching contigs found!\n")
}
cat("Process completed!\n")


##Combining all interesting contigs in one fasta file 
fasta_file <- "PSL-8W-1.fasta"  #Update according to file path
sequences <- readDNAStringSet(fasta_file)

contigs_list <- readLines("interesting_contigs.txt")

#Filtering the sequences
filtered_sequences <- sequences[names(sequences) %in% contigs_list]

#Writing the filtered sequences to a new FASTA file:
writeXStringSet(filtered_sequences, "filtered_contigs.fasta")

#Special section for contig 1154
fasta_file <- "contig_1154_W.fasta"
sequences <- readDNAStringSet(fasta_file)

selected_contig <- unique(sequences[names(sequences) == "contig_1154"])

writeXStringSet(selected_contig, "contig_1154.fasta")





#Splitting the longer contigs into specific genes to synthesise
#Defining the file path
file_path <- file.path("untrimmed_contigs", "contig_1251.fasta")

contig <- readDNAStringSet(file_path)

#Extracting the name
contig_name <- names(contig)

#Define the cut-positions
start_pos <- 495454            #Update according to correct positions
end_pos <- 496533

extracted_sequence <- subseq(contig[[1]], start = start_pos, end = end_pos)

output_file <- file.path("profiles", paste0("profile_", contig_name, "_1.fasta"))

writeXStringSet(DNAStringSet(extracted_sequence), filepath = output_file)

#Count length of cut sequence
profile_file.path <- file.path("profiles","profile_contig_1251_1.fasta")
profile_contig <- readDNAStringSet(profile_file.path)
contig_length <- nchar(profile_contig[[1]])

#Printing the length for validation 
cat("The contig contains", contig_length, "base pairs.\n")

