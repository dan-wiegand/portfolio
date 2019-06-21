#Midterm Project for MATH E156 - Harvard University - Daniel Wiegand - Matthew Smith

# The purpose of this script is to import DNA sequences from Illumina Next Generation Gene Sequencing data files
# in order to evaluate the quality of an oligonucleotide library synthesized from a CustomArray 12K
# Oligonucleotide Library Synthesis Array.

# Download Rsamtools and VariantAnootation to read sequence alignment files from UseGalaxy.org
source("https://bioconductor.org/biocLite.R")
biocLite(c("Rsamtools","VariantAnnotation","Biostrings"))
install.packages("gplots", dependencies = TRUE)
install.packages("RColorBrewer", dependencies = TRUE)

# The Rsamtools package provides an interface to BAM files. BAM files are produced by samtools and other software, and
# represent a flexible format for storing ‘short’ reads aligned to reference genomes. BAM files typically contain sequence and
# base qualities, and alignment coordinates and quality measures.

# Call up libraries for 
library(Rsamtools)
library(VariantAnnotation)
library(Biostrings)
library(gplots)
library(RColorBrewer)

# Define function for creating frequencies of nucleotides from aligned sequencing
# data from BAM data file and BAI indexing file generated with the UseGalaxy.org tool
# box. This global function uses the pileup function from the Rsamtools and VariantAnnotation libraries

pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}

gcContent <- function(x){
  alf <- alphabetFrequency(x, as.prob=TRUE)
  sum(alf[c("G", "C")])
}

calcGC <- function(referenceSeq){
  #Calculate GC content
  refGC <- numeric(length(referenceSeq))
  for (i in 1:length(refGC)){
    dString <- DNAString(referenceSeq[i])
    refGC[i] <- gcContent(dString)
  }
  return(refGC)
}

importReference <- function(x){
  referenceLib <- readDNAStringSet(x) #Import reference library
  referenceName <- names(referenceLib) #Call names
  referenceSeq <- paste(referenceLib) #Call sequences
  refFrame <- data.frame(referenceName,referenceSeq)
  return(refFrame)
}

importBAM <- function(bamfile){
  bamfile <- bamfile
  bf <- BamFile(bamfile)
  what <- c("pos", "strand", "rname", "seq", "qual", "qname")
  param <- ScanBamParam(what=what)
  p_param <- PileupParam(max_depth = 1000, ignore_query_Ns = FALSE)
  res <- pileup(bf,scanBamParam = param, pileupParam = p_param)
  return(res)
}

libraryErrorAnalysis <- function(referenceName){
  
  #Declare needed variables
  libraryNames <- unique(dataRef$referenceName)
  readTotal <- numeric(length(libraryNames))
  
  errorTotal <- numeric(length(libraryNames))
  errorMean <- numeric(length(libraryNames))
  errorSD <- numeric(length(libraryNames))
  
  normError <- numeric(length(libraryNames))
  normSD <- numeric(length(libraryNames))
  
  mutTotal <- numeric(length(libraryNames))
  insTotal <- numeric(length(libraryNames))
  delTotal <- numeric(length(libraryNames))
  
  #for large datasets, this script takes a while - about 7 hours for the large dataset
  for (k in 1:length(libraryNames)){
    #call name of sequence being analyzed
    name <- libraryNames[k] 
    #extract the relevant information for this sequence from the BAM file
    subError <- subset(res, seqnames == name, select=c(seqnames,pos,strand,nucleotide,count))
    #count the number of total rows
    rows <- nrow(subError)
    #some sequences were not present in the sequencing data, make a if statement for handling failures
    if (rows > 1) {
      #pass extracted data to pileupFreq global function
      error.pileup <- pileupFreq(subError)
      #extract only the integers the sequence nucleotide frequences
      freqError <-subset(error.pileup , seqnames == name, select=c("A","C","G","T","N","=","-","+"))
      #begin summation of total errors
      N<-numeric(nrow(freqError));m<-numeric(nrow(freqError));ins<-numeric(nrow(freqError));del<-numeric(nrow(freqError));nmean<-numeric(nrow(freqError))
      for (i in 1:length(N)){
        # the highest frequency of nucleotide at each position is assumed to be the correct nucleotide
        # subtract the sum of all nucleotide counts, which the the number of total reads for a given 
        # sequence from this number to yield the total number of errors per sequence
        N[i] <- sum(c(sum(c(freqError[i,1],freqError[i,2],freqError[i,3],freqError[i,4]))-max(freqError[i,c(1,2,3,4)])),freqError[i,8],freqError[i,7])
        m[i] <- sum(c(freqError[i,1],freqError[i,2],freqError[i,3],freqError[i,4]))-max(freqError[i,c(1,2,3,4)])
        ins[i] <- freqError[i,8]
        del[i] <- freqError[i,7]
        nmean[i] <- (sum(c(sum(c(freqError[i,1],freqError[i,2],freqError[i,3],freqError[i,4]))-max(freqError[i,c(1,2,3,4)])),freqError[i,8],freqError[i,7]))/(sum(freqError[i,])*170)
      }
      
      #get the statistics for later; sample from this data extraction
      
      # number of errors in total
      errorTotal[k] <- sum(N) #extract the total number of errors for each position and sum
      errorMean[k] <- mean(N[N>0]) #mean of the total number of errors leaving out zeros
      errorSD[k] <- sd(N[N>0]) #calculate the standard deviation leaving out the zeros
      
      # normalized error (error per number of bases sequenced = reads*170)
      normError[k] <- mean(nmean[nmean>0])
      normSD[k] <- sd(nmean[nmean>0])
      
      #get the statistics for the individual errors needed for barplots
      mutTotal[k] <- sum(m)
      insTotal[k] <- sum(ins)
      delTotal[k] <- sum(del)
      
      #print sequence name when finished, so to indicate that the analysis for a particular sequence was completed
      print(name)
    }
  }
  #turn extracted data into data frame and export as a csv file
  output.ngsdata <- data.frame(referenceName,refGC,errorTotal,errorMean,errorSD,normError,mutTotal,insTotal,delTotal)
  write.table(output.ngsdata, file="error_extraction_NGS.csv", sep=",")
}

# Import necessary files, calculate GC content, and perform error analysis.
# Output is a csv file containing all the information
refData <- importBAM("ngsdata.bam")
dataRef <- importReference("reference_lib.fasta")
refGC <- calcGC(dataRef$referenceSeq)
# libraryErrorAnalysis(referenceName) #creates a data file in main directory containing the information pulled by the function


### Begin statistical analysis here ----------------------------------------------------------------------------

# Variable Key from Extracted Sequencing Data:
# referenceName: name of the sequence
# refGC: GC content of sequence
# errorTotal: the total number of errors, not normalized for number of reads
# errorMean: mean of the total error, does not take missing sequences (0) into account
# errorSD: standard deviation of Total error, does not take missing sequences into account
# normError: the normalized error: totalErrors/number of reads * length of sequence (170)
# mutTotal: the number of mutations (like G to A, etc)
# insTotal: the total number of insertions (adding addition sequence)
# delTotal: the total number of deletions (bases missing)
# gcBin: the bins for high, medium and low -> high is above 65% and low is below 45%
# arrayQuad: the quadrant in which the DNA sequence was synthesized in
# Note the expected error rate is 1/200 or 0.005 for a given sequence

###
### Data import and setup (DATASET POINTS: 1-4, BONUS POINTS: 1,2,14)
###

# Import data so you don't have to rerun the long script again for stats analysis, change zeros to NA's so zeros don't interfere
extractionAll <- read.csv("error_extraction_NGS.csv"); extData<-extractionAll[1:5152,]; extData[extData==0] <- NA

# Add bins for gc distribution for contigency table and analysis
gcBin <- numeric(length(extData$refGC))
for (i in 1:length(gcBin)){
  gc <- extData$refGC[i]
  if (gc < 0.45){ 
    gcBin[i] = "Low"
  }else if (gc > 0.65){
    gcBin[i] = "High"
  }else
    gcBin[i]= "Medium"
}
extData <- cbind(extData, gcBin)

# Set up array quadrants for contigency table and analysis 
arrayQuad <- c(rep(1,1288),rep(2,1288),rep(3,1288),rep(4,1288))
extData <- cbind(extData, arrayQuad)

# include number of reads per sequence
reads = extData$errorTotal / (extData$normError * 170)
extData = cbind(extData, reads)

###
### Function declarations (BONUS POINT: 8)
###

# chisquare functions for test statistics based on observations as in textbook section 3.4.1.
chisq = function(Obs){
  expected = outer(rowSums(Obs), colSums(Obs)) / sum(Obs)
  sum((Obs - expected)^2/expected)
}

# chisquare permutation for analysis of contingency table later
chi.perm = function(input.1, input.2) {
  N = 10^4 - 1
  result = numeric(N)
  for (i in 1:N){
    permutate = sample(input.1)
    new.table = table(input.2, permutate)
    new.table = new.table[,c(2,3,1)]
    result[i] = chisq(new.table)
  }
  return(result)
}

# inner permutation abstraction for x = data to choose from, s = size to choose, foo = data sample
permutation = function(x, s, foo){
  N = 10^4 - 1
  result = numeric(N)
  for (i in 1:N){
    index = sample(x, size = s, replace = FALSE)
    result[i] = mean(foo[index]) - mean(foo[-index])
  }
  return(result)
}

# outer abstracted permutation function for analysis of multiple variables by permutation test 
perm = function(i){
  a = t[t$arrayQuad == val[1,i],]
  b = t[t$arrayQuad == val[2,i],]
  c = rbind(a,b)
  obs = mean(a$normError) - mean(b$normError)
  r = permutation(NROW(c), NROW(a), c$normError)
  s = ((sum(r < obs) + 1) / 10^4)
  hist(r, col = "light blue", main = c("Histogram of Results: ", val[1,i], " ", val[2,i]),
       xlab = "Sample Results") 
  abline(v = obs, col = "red")
  return(s)
}

###
### Statistical Analysis (ANALYSIS POINT: 1,3-4, GRAPHICAL DISPLAY POINT: 2-4)
###

### Permutation tests for independence between normalized error rates and 4 array locations
t = extData[,c("arrayQuad", "normError")]
t[is.na(t)] = 0
val = combn(4,2)
p = numeric(6)
par(mfrow = c(2, 3))
for (j in 1:6) {
  p[j] = perm(j)
}
names(p) = c("1,2", "1,3", "1,4", "2,3", "2,4", "3,4"); p
par(mfrow = c(1, 1))
# COMMENT ON SIGNIFICANCE OF RESULTS?

### Quick clean-up of variables
rm(t, val, i, j, gc, reads, extractionAll, arrayQuad)

### P-value based on a distribution function


### Analysis of the contingency table
array = extData$arrayQuad
bins = extData$gcBin
CT = table(array, bins) 
CT = CT[,c(2,3,1)]; CT
observed = chisq(CT); observed
r = chi.perm(array, bins)
hist(r, breaks = 40, col = "navy", main = "Contingency Table Analysis", col.main = "maroon", probability = TRUE)
abline(v = observed, col = "orange", lwd = 2)

# fits nicely with chisquare density function
curve(dchisq(x,4), col = "red", add = TRUE, lwd = 2)

# simulated p-value
p = (sum(r >= observed) + 1) / 10^4; p
# compared to:
chisq.test(CT)$p.value

# Observed statistic is similar; however the p-value in the simulated version does not come out
#   with the same result as does the built-in function, but at the 10% confidence interval, both
#   have the same conclusion. As time allows, increasing N value in chi.perm() function improves p.

# calculation of kurtosis (should be positive) and skewness (should be positive) by the chisquare density function
integ = integrate(function(x) x*dchisq(x,4), -Inf, Inf)
mu = integ$value
MC2 = integrate(function(x) (x-mu)^2*dchisq(x,4), -Inf, Inf)$value
MC3 = integrate(function(x) (x-mu)^3*dchisq(x,4), -Inf, Inf)$value
sigma = sqrt(MC2); sigma
skewness.1 = MC3/sigma^3; skewness.1
MC4 = integrate(function(x) (x-mu)^4*dchisq(x,4), -Inf, Inf)$value
kurtosis.1 = MC4/sigma^4-3; kurtosis.1

# same calculations for r, the sampled data in our chisquare permutation of the contingency table analysis
mu = mean(r)
MC2 = mean((r-mu)^2)
MC3 = mean((r-mu)^3)
sigma = sqrt(MC2)
skewness.2 = MC3/sigma^3; skewness.2
MC4 = mean((r-mu)^4)
kurtosis.2 = MC4/sigma^4-3; kurtosis.2

# difference between the two:
kurtosis.1 - kurtosis.2
skewness.1 - skewness.2
# this again shows that the data aligns with the chisquare density function fairly well!

### Quick clean-up of variables
rm(array, bins, CT, observed, r, p, mu, MC2, MC3, sigma, skewness, MC4, kurtosis)

###
### Graphical Displays (GRAPHICAL DISPLAY POINTS: 1, BONUS POINTS: 3, 9)
###

### barplot of frequency of GC content bins
gcBin = table(gcBin)
gcBin = gcBin[c(2,3,1)]
barplot(gcBin, main = "Frequency of GC Content Bins", col = c("green", "yellow", "red"),
        col.main = "dark green", family = "serif", font = 4)

### boxplot of normalized error rates by the GC content bin levels
l = extData$gcBin
l = factor(l, c("Low", "Medium", "High"))
boxplot(extData$normError ~ l, main = "Normalized Error Rates by GC Content Bin",
        col = c("green", "yellow", "red"), col.main = "maroon")

### QQ Normal Quantile Plots
gc = extData$refGC
qqnorm(gc, col = c("blue", "red"), main = "Normal Quantile Plot of GC Content Values")
qqline(gc, lwd = 2)
# nearly a normal quantile plot as seen by the fit of the line, at least near the [-2,2] range

### heatmap display of errors

# Begin spatial analysis of 12k Array & generate heatmaps
quadAll <-matrix(extData$normError, nrow=92, ncol=56, byrow = TRUE)

# Split main array into 4 quadrants to compare sequence errors based on position
k <- kronecker(matrix(1:4, 2, byrow = TRUE), matrix(1, 46, 28))
splitQuad <- lapply(split(quadAll, k), matrix, nr = 46) 
quad1 <- as.vector(t(splitQuad$`1`))
quad2 <- as.vector(t(splitQuad$`2`))
quad3 <- as.vector(t(splitQuad$`3`))
quad4 <- as.vector(t(splitQuad$`4`))

my_palette <- colorRampPalette(c("blue", "green", "red"))(n = 5152)

heatmap.2(splitQuad$`4`,
          main="Heatmap of Errors on CustomArray 12K OLS Arrays", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          dendrogram="none" ,   # turns off dendrogram
          trace="none",         # turns off trace lines inside the heat map
          margins=c(12,9),      # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          Colv=FALSE,           # turn off column clustering
          Rowv=FALSE,           # turn off row clustering
          key.xlab = "Array Placement") # label Color Key x-axis
