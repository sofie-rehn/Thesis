#Load libraries
library(vegan)
library(dada2)
library(Biostrings)
library(ShortRead)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(edgeR)

#Manually creates a table of sample names
SampleName = c("PSL-0W-1","PSL-0W-2","PSL-2W-1","PSL-2W-2",
                 "PSL-4W-1","PSL-4W-2","PSL-8W-1","PSL-8W-2",
                 "PSL-6M-1","PSL-6M-2")

#Create a list of the raw data files
path.raw <- file.path("data", "raw") 
list.files(path.raw) 
barcode_files <- list.files(path.raw, pattern="barcode",full.names = TRUE) 

#Match the raw data files and sample names
seqID <- setNames(barcode_files, SampleName)

#Check if the data files was correctly matched to the sample names
print(seqID) 

#Check the quality profile of the raw data for a selected sample (random) 
plotQualityProfile(seqID[3])

#Count number of reads and length of reads before doing any trimming, cut or filtering
for (file.path in barcode_files){
  fastq_data <- readFastq(file.path) #Read the fastq file for each file path
  
  num_reads <- length(fastq_data) #Count the number of reads
  print(paste("Number of Reads",file.path,":", num_reads)) 
  
  seq_lengths <- width(fastq_data) #Count the length of Reads
  print(summary(seq_lengths) ) 
}

##Identify and count primers (validation)
#Define the primer sequences
FWD <- "AGRGTTYGATYMTGGCTCAG"
REV <- "GACGGGCGGTGWGTRCA"

#Create all orientations of the primer sequences
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#Pre-filtering by removing undetermined bases (Ns)
#Create new file paths and save the trimmed sequences in new files
barcode_files.filtN <- file.path("data", "filtN", paste0(SampleName,"_filtN.fastq.gz"))
path.filtN <- file.path("data","filtN")

#Remove the Ns
filterAndTrim(barcode_files, barcode_files.filtN, maxN = 0)

#Function for counting the number of times the primers appears
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#Actually counting the number of primers in selected (random) sample files
rbind(FWD.Reads = sapply(FWD.orients, primerHits, fn = barcode_files.filtN[[9]]), 
      REV.Reads = sapply(REV.orients, primerHits,fn = barcode_files.filtN[[9]]))

rbind(FWD.Reads = sapply(FWD.orients, primerHits, fn = barcode_files.filtN[[7]]), 
      REV.Reads = sapply(REV.orients, primerHits,fn = barcode_files.filtN[[7]]))

#Removing primers using cutadapt
#Let R access cutadapt
cutadapt <- "/opt/anaconda3/envs/cutadapt/bin/cutadapt" #Modify file path according to device
system2(cutadapt, args = "--version")

#Creating file-names for cutadapt output and defining trimming parameters
path.cut <- file.path("data", "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
barcode_files.cutP <- file.path(path.cut, paste0(SampleName,"_cutP.fastq.gz"))

#Define the reverse complement of the primer sequences
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

#Trim the forward primer and the reverse-complement of reverse primer
reads.flags <- c("-g", FWD, "-a", REV.RC,"-g", REV, "-a", FWD.RC)

# Run Cutadapt
for(i in seq_along(barcode_files.filtN)) {
  system2(cutadapt, args = c("-u","24",#remove custom 24 nucleotides (barcodes)
                             reads.flags,
                             "-n", 3, #required to remove FWD and REV from reads, and not only best match
                             "-e", 0.2,
                             "-o", barcode_files.cutP[i], #output file
                             barcode_files.filtN[i])) #input file
}

#Count the number of primers after removal (compare manually to the before-count)
rbind(FWD.Reads = sapply(FWD.orients, primerHits, fn = barcode_files.cutP[[9]]), #Update what files are counted depending on purpose.
      REV.Reads = sapply(REV.orients, primerHits,fn = barcode_files.cutP[[9]]))

rbind(FWD.Reads = sapply(FWD.orients, primerHits, fn = barcode_files.cutP[[7]]), #Update what files are counted depending on purpose.
      REV.Reads = sapply(REV.orients, primerHits,fn = barcode_files.cutP[[7]]))


# Trimming the raw data according to determined quality parameters
#Set a new file path and save the output in new files
path.filt <- file.path("data", "filtered")
barcode_files.filtered <- file.path(path.filt, paste0(SampleName,"_filtered.fastq.gz"))

#Filter and trim according to quality parameters
out <- filterAndTrim(barcode_files.cutP, barcode_files.filtered, maxN = 0, maxEE = 3, truncQ = 2, 
                     minLen = 50, rm.phix = TRUE)

#Examine the quality plots for two samples
plotQualityProfile(barcode_files.filtered[9])
plotQualityProfile(barcode_files.filtered[7])

#Count the number of reads and length of sequences again
for (file.path in barcode_files.filtered){
  fastq_data <- readFastq(file.path) #Read the fastq file for each file path
  
  num_reads <- length(fastq_data) #Count the number of reads
  print(paste("Number of Reads",file.path,":", num_reads))
  
  seq_lengths <- width(fastq_data) #Count the length of Reads
  print(summary(seq_lengths) )
}


# Learn and plot the error rates
#NOTE: the processing time for this line of code is LONG
errF <- learnErrors(barcode_files.filtered, multithread = TRUE)
plotErrors(errF, nominalQ=TRUE)

#Check number for unique sequences for the filtered samples... 
#...to determine if dada is an appropriate method to use
derepFastq(barcode_files.filtered[1])
derepFastq(barcode_files.filteredN[2])
derepFastq(barcode_files.filtered[3])

#Sample inference, applied to the dereplicated data
dadaBarcode_files.filtered <- dada(barcode_files.filtered, err = errF, multithread = TRUE)

#Construct sequence table
seqtab <- makeSequenceTable(dadaBarcode_files.filtered)

#Remove chimeras, sequencing artifacts not in the original sample
seqtab.nochim = removeBimeraDenovo(
  seqtab, method="consensus", multithread=TRUE, verbose=TRUE
)

#Check the fraction of chimeras (validation)
sum(seqtab.nochim)/sum(seqtab)
dim(seqtab)

#Define variables for the taxonomic profiling
seqtab_copy <- seqtab
seqtab = t(seqtab)
colnames(seqtab)
rownames(seqtab)
asv = rownames(seqtab)
rownames(seqtab) = paste("ASV", 1:length(asv), sep = "")
rownames(seqtab)

#Inspect distribution of sequence length 
table(nchar(getSequences(seqtab.nochim)))

#Tracking reads through the pipeline (validation)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaBarcode_files.filtered, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- SampleName
head(track)

#Assigning taxonomy 
set.seed(123) 
taxa = assignTaxonomy(
  asv, "data/silva/silva_nr99_v138.1_train_set.fa.gz",
  multithread=TRUE,
  taxLevels = c(
    "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )
)

#naming the rows of "taxa" as the rows of seqtab, ie with ASV IDs
rownames(taxa) = rownames(seqtab)

#Counting the number of undetermined outputs for "genus"
taxa.NA <- is.na(taxa[,"Genus"])
taxa.NAtrue <- which(taxa.NA == TRUE)
taxa.NAfalse <- which(taxa.NA == FALSE)
print(length(taxa.NAtrue))
print(length(taxa.NAfalse))
print(length(taxa[,"Genus"]))

#Normalising ASV count
norm_seqtab = seqtab
for (i in 1:ncol(seqtab)) {
  norm_seqtab[,i] = seqtab[,i]/sum(seqtab[,i])
}

clade_counts = list()
norm_clade_counts = list()

for (i in 1:ncol(taxa)) {
  matr = norm_matr = NULL
  clade = unique(taxa[,i])
  clade = clade[!is.na(clade)]
  for (j in 1:length(clade)) {
    ix = which(clade[j]==taxa[,i])
    if (length(ix) > 1) {
      matr = rbind(matr, apply(seqtab[ix,], 2, sum, na.rm=TRUE))
      norm_matr = rbind(norm_matr, apply(norm_seqtab[ix,], 2, sum, na.rm=TRUE))
    } else {
      matr = rbind(matr, seqtab[ix,])
      norm_matr = rbind(norm_matr, norm_seqtab[ix,])
    }
  }
  rownames(matr) = rownames(norm_matr) = clade
  colnames(matr) = colnames(norm_matr) = SampleName
  clade_counts[[i]] = matr
  norm_clade_counts[[i]] = norm_matr
}

#Validation of clade counts 
clade_counts[[2]][,1]
clade_counts[[6]][3,1:10]


#Creating the taxonomy barplots 
# set what taxonomic level to plot (1 - 6, corresponding to domain - genus)
tax_level = 6

# to select those clades with a relative abundance over a threshold (here 0.01)
ok = which(apply(norm_clade_counts[[tax_level]], 1, mean) > 0.01)

# make the barplot
barplot(
  norm_clade_counts[[tax_level]][ok,],
  col = mycols(length(ok)),
  las = 2,
  names.arg = paste(SampleName)
)

# add a color legend
legend(
  "bottomleft", bty = "n", pch = 19,
  col = mycols(length(ok))[1:length(ok)],
  cex = 1, inset = c(1,0),
  legend = rownames(clade_counts[[tax_level]])[ok]
)

# to make a color palette
mycols = colorRampPalette(c("#a6cee3",
                            "#1f78b4",
                            "#b2df8a",
                            "#33a02c",
                            "#fb9a99",
                            "#e31a1c",
                            "#fdbf6f",
                            "#ff7f00",
                            "#cab2d6",
                            "#6a3d9a",
                            "#ffff99",
                            "#b15928"))

# define the plotting area
par(mfrow=c(1,1), mar=c(7,3,2,10), xpd = TRUE)

# Manually define which samples belong to each incubation time group
zeroW <- which(colnames(clade_counts[[1]]) %in% c("PSL-0W-1", "PSL-0W-2"))
twoW <- which(colnames(clade_counts[[1]]) %in% c("PSL-2W-1", "PSL-2W-2"))
fourW <- which(colnames(clade_counts[[1]]) %in% c("PSL-4W-1", "PSL-4W-2"))
eightW <- which(colnames(clade_counts[[1]]) %in% c("PSL-8W-1", "PSL-8W-2"))
sixM <- which(colnames(clade_counts[[1]]) %in% c("PSL-6M-1", "PSL-6M-2"))

# Final custom sample order (like the tutorial's 'all' vector)
all <- c(zeroW, twoW, fourW, eightW, sixM)
sample_type <- c("zeroW", "zeroW", "twoW", "twoW", "fourW", "fourW", "eightW", "eightW", "sixM", "sixM")

#Redoing the barplot but ordering according to incubation time
# define the plotting area
par(mfrow=c(1,1), mar=c(7,5,3,13.5), xpd = TRUE)

# set what taxonomic level to plot (1 - 6)
tax_level = 6

# make the barplot
ok = which(apply(norm_clade_counts[[tax_level]], 1, mean) > 0.005)
barplot(
  norm_clade_counts[[tax_level]][ok,all],
  col = mycols(length(ok)),
  las = 2,
  names.arg = paste(SampleName[all]),
  ylab = "Relative abundance",
  xlab = "",
  main = "Taxonomic level: Genus"
)

mtext("Sample Name", side = 1, line = 5.5)

# add a color legend
legend(
  "bottomleft", bty = "n", pch = 19,
  col = mycols(length(ok))[1:length(ok)],
  cex = 1, inset = c(1,0),
  legend = rownames(clade_counts[[tax_level]])[ok]
)

#Counting the number of undetermined outputs for each taxa level
taxa.NA <- is.na(taxa[,"Genus"]) #Update according to taxa level of interest
taxa.NAtrue <- which(taxa.NA == TRUE)
taxa.NAfalse <- which(taxa.NA == FALSE)
print(length(taxa.NAtrue))
print(length(taxa.NAfalse))
print(length(taxa[,"Phylum"]))


#Calculating the alpha-diversity (within-sample diversity)
shannon = diversity(seqtab.nochim, MARGIN = 2)

#Making a barplot of the Shannon diversities.
barplot(
  shannon[all], las = 2,
  names.arg = SampleName[all]
)

#Making boxplots, one per sample incubation time
boxplot(
  shannon[zeroW], shannon[twoW],
  shannon[fourW], shannon[eightW],
  shannon[sixM],
  names = c(
    "0 weeks", "2 weeks",
    "4 weeks", "8 weeks",
    "6 months"
  ),
  las = 2,
  ylab = "Shannon diversity index",
  xlab = "",
  main = "Alpha diversity over time"
)

mtext("Incubation time", side = 1, line = 5)

group_names <- unique(sample_type)
group_positions <- as.numeric(factor(sample_type, levels = group_names))
sample_label <- c("0 week", "2 week", "4 week", "8 week", "6 month")

par(mar = c(7, 6, 4, 6))

# Basic scatterplot
plot(
  group_positions, shannon,
  xaxt = "n",
  xlab = "Incubation time",
  ylab = "Shannon diversity index",
  main = "Alpha diversity by sample",
  pch = 19
)
axis(1, at = 1:length(group_names), labels = sample_label)


#checking if the shannon diversity between samples are statistically significant 
#Wilcoxon signed-rank test
wilcox.test(shannon[zeroW], shannon[twoW], paired = TRUE)
wilcox.test(shannon[zeroW], shannon[fourW], paired = TRUE)
wilcox.test(shannon[zeroW], shannon[eightW], paired = TRUE)
wilcox.test(shannon[zeroW], shannon[sixM], paired = TRUE)
wilcox.test(shannon[twoW], shannon[fourW], paired = TRUE)
wilcox.test(shannon[twoW], shannon[eightW], paired = TRUE)
wilcox.test(shannon[twoW], shannon[sixM], paired = TRUE)
wilcox.test(shannon[fourW], shannon[eightW], paired = TRUE)
wilcox.test(shannon[fourW], shannon[sixM], paired = TRUE)
wilcox.test(shannon[eightW], shannon[sixM], paired = TRUE)

#Investigating beta-diversity (between-sample community differences) - Bray Curtis
bray_dist = as.matrix(vegdist(t(norm_seqtab.nochim), method = "bray"))

pheatmap(
  bray_dist,
  clustering_distance_rows = as.dist(bray_dist),
  clustering_distance_cols = as.dist(bray_dist),
  labels_row = paste(SampleName),
  labels_col = paste(SampleName),
  main = "\n\nBeta diversity (Bray Curtis) "
)

#Creating Multi-dimensional clustering graphs 
color_type = rep("white", ncol(norm_seqtab.nochim))
color_type[1:2] = "lightgreen"
color_type[3:4] = "purple"
color_type[5:6] = "blue"
color_type[7:8] = "orange"
color_type[9:10] = "hotpink"

pch_type = rep(21, length(all))

#running the mds algorithm
mds = metaMDS(bray_dist)

mds_data <- data.frame(MDS1 = mds$points[,1], MDS2 = mds$points[,2], 
                       color = color_type, labels = SampleName)
#Plotting the results
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = color)) +
  geom_point(aes(fill = color), shape = 21, size = 3) +  #Circles
  scale_fill_identity() +  #Using defined colours
  scale_color_identity() +
  geom_text_repel(aes(label = labels), box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  #Avoiding label overlap
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2")
  theme(legend.position = "none") #Hiding legend

#Checking if the clustering of samples in NMDS is statisically significant
#r-value close to 1 suggests dissimilarity between groups are larger than within groups
#0 suggests an even distribution of dissimilarities within and between groups

#Defining the sample indices for each time point (this assumes the samples are in specific columns)
time1 <- c(1, 2)      # For time 0 weeks (samples 1 and 2)
time2 <- c(3, 4)       # For time 2 weeks (samples 3 and 4)
time3 <- c(5, 6)      # For time 4 weeks (samples 5 and 6)
time4<- c(7, 8)     # For time 8 weeks (samples 7 and 8)
time5<- c(9, 10)      # For time 6 months (samples 9 and 10)

# Create a vector of group labels corresponding to each time point
timestamp <- c(rep("time1", length(time1)),
               rep("time2", length(time2)),
               rep("time3", length(time3)),
               rep("time4", length(time4)),
               rep("time5", length(time5)))

all_samples <- c("time1", "time2", "time3", "time4", "time5")

anosim_results <- anosim(bray_dist[all_samples,all_samples], timestamp, permutations = 9999)
print(anosim_results)

#Checking if any specific taxa differs between the sample groups 
#Differential abundance analysis 
x <- clade_counts[[6]][, all_samples]

group <- factor(c(
  rep("0W", 2),
  rep("2W", 2),
  rep("4W", 2),
  rep("8W", 2),
  rep("6M", 2)
))

#running the statistical analysis (on genus level)
x = clade_counts[[6]][, all_samples]
y = DGEList(counts=x, group=group)
y = calcNormFactors(y)
design = model.matrix(~group)
y = estimateDisp(y, design)
fit = glmQLFit(y, design)

#Perfoming the abundance analysis
# 0W 2W
contrast_0W2W <- makeContrasts(Intercept - group2W, levels = design)
qlf_0W2W <- glmQLFTest(fit, contrast = contrast_0W2W)

top_0W2W <- topTags(qlf_0W2W, n=1000, p.value=0.05)
sign_clades_0W2W = rownames(top_0W2W$table)
# 0W 4W
contrast_0W4W <- makeContrasts(Intercept - group4W, levels = design)
qlf_0W4W <- glmQLFTest(fit, contrast = contrast_0W4W)

top_0W4W <- topTags(qlf_0W4W, n=1000, p.value=0.05)
sign_clades_0W4W = rownames(top_0W4W$table)
# 0W 8W
contrast_0W8W <- makeContrasts(Intercept - group8W, levels = design)
qlf_0W8W <- glmQLFTest(fit, contrast = contrast_0W8W)

top_0W8W <- topTags(qlf_0W8W, n=1000, p.value=0.05)
sign_clades_0W8W = rownames(top_0W8W$table)
# 0W 6M
contrast_0W6M <- makeContrasts(Intercept - group6M, levels = design)
qlf_0W6M <- glmQLFTest(fit, contrast = contrast_0W6M)

top_0W6M <- topTags(qlf_0W6M, n=1000, p.value=0.05)
sign_clades_0W6M = rownames(top_0W6M$table)
# 2W 4W
contrast_2W4W <- makeContrasts(group2W - group4W, levels = design)
qlf_2W4W <- glmQLFTest(fit, contrast = contrast_2W4W)

top_2W4W <- topTags(qlf_2W4W, n=1000, p.value=0.05)
sign_clades_2W4W = rownames(top_2W4W$table)
# 2W 8W
contrast_2W8W <- makeContrasts(group2W - group8W, levels = design)
qlf_2W8W <- glmQLFTest(fit, contrast = contrast_2W8W)

top_2W8W <- topTags(qlf_2W8W, n=1000, p.value=0.05)
sign_clades_2W8W = rownames(top_2W8W$table)
# 2W 6M
contrast_2W6M <- makeContrasts(group2W - group6M, levels = design)
qlf_2W6M <- glmQLFTest(fit, contrast = contrast_2W6M)

top_2W6M <- topTags(qlf_2W6M, n=1000, p.value=0.05)
sign_clades_2W6M = rownames(top_2W6M$table)
# 4W 8W
contrast_4W8W <- makeContrasts(group4W - group8W, levels = design)
qlf_4W8W <- glmQLFTest(fit, contrast = contrast_4W8W)

top_4W8W <- topTags(qlf_4W8W, n=1000, p.value=0.05)
sign_clades_4W8W = rownames(top_4W8W$table)
# 4W 6M
contrast_4W6M <- makeContrasts(group4W - group6M, levels = design)
qlf_4W6M <- glmQLFTest(fit, contrast = contrast_4W6M)

top_4W6M <- topTags(qlf_4W6M, n=1000, p.value=0.05)
sign_clades_4W6M = rownames(top_4W6M$table)
# 8W 6M
contrast_8W6M <- makeContrasts(group8W - group6M, levels = design)
qlf_8W6M <- glmQLFTest(fit, contrast = contrast_8W6M)

top_8W6M <- topTags(qlf_8W6M, n=1000, p.value=0.05)
sign_clades_8W6M = rownames(top_8W6M$table)

#Summarising the results
all_sign_taxa <- unique(c(
  sign_clades_0W2W, sign_clades_0W4W, sign_clades_0W8W, sign_clades_0W6M,
  sign_clades_2W4W, sign_clades_2W8W, sign_clades_2W6M,
  sign_clades_4W8W, sign_clades_4W6M, sign_clades_8W6M
))

#Plotting the heatmap
heatmap_matrix <- norm_clade_counts[[6]][all_sign_taxa, all_samples]


pheatmap(
  heatmap_matrix,
  cluster_cols=TRUE, cluster_rows=FALSE,
  labels_row = rownames(heatmap_matrix),
  cex = 0.8,
  labels_col = SampleName[all_samples],
  clustering_distance_rows = "correlation",
  main = "Differentially Abundant Genera Across Time"
)