# Author: Daniel Wiegand
# Institution: Harvard Medical School, Wyss Institute
# Purpose: Analysis of Kinetic Data of In Vitro Protein Synthesis

# Load data from file
experiment.data <- read.csv("invitro_protein_syn_10302017.csv")

initialVelocity <- function(experiment.data,species,cofactor){
  # Initialize the data frame, selecting the appropriate information based on the user input
  x <- data.frame(t(subset(experiment.data, Species == species & Cofactor == cofactor, select = c(X1:X73))))
  #This function requires two replicates 
  colnames(x) <- c("S1","S2","S3"); rownames(x) <- c(1:73); newx <- c(1:73)
  samples <- subset(x, select = c(S1, S2, S3))
  
  # Determine the first derv of fitted function to determine slope for replicates
  spl.s1 <- smooth.spline(samples$S1 ~ newx); spl.s2 <- smooth.spline(samples$S2 ~ newx); spl.s3 <- smooth.spline(samples$S3 ~ newx)
  pred1 <- predict(spl.s1, x=newx, deriv=1); pred2 <- predict(spl.s2, x=newx, deriv=1); pred3 <- predict(spl.s3, x=newx, deriv=1)
  
  # Make plot space 2x1 and plot the raw kinetic data as well as the slopes at each point for the first replicates
  par(mfrow = c(1,2))
  plot(newx,samples$S1, xlab = "Time (min)", ylab = "RFU", main = paste(species,": Raw Kinetic Data for Mg @", cofactor, "mM"))
  plot(newx,pred1$y, xlab = "Time (min)", ylab = "δRFU/δT", main = paste(species,": Raw Kinetic Data for Mg @", cofactor, "mM"))
  
  #Determine the mean and SEM of the replicates
  output1 <- max(pred1$y); output2 <- max(pred2$y); output3<- max(pred3$y)
  pred.mean <- mean(c(output1,output2,output3))
  pred.sem <- sd(c(output1,output2,output3))/sqrt(2)
  
  # Function output is a vector containing these
  calculation <- c(pred.mean,pred.sem)
  return(calculation) # No code to handle if reaction is poor just yet
}

MaxRFU <- function(experiment.data,species,cofactor){
  # Initialize the data frame, selecting the appropriate information based on the user input
  x <- data.frame(t(subset(experiment.data, Species == species & Cofactor == cofactor, select = c(X1:X73))))
  #This function requires two replicates 
  colnames(x) <- c("S1","S2","S3"); rownames(x) <- c(1:73); newx <- c(1:73)
  samples <- subset(x, select = c(S1, S2, S3))

  #Determine the mean and SEM of the replicates
  output1 <- max(samples$S1); output2 <- max(samples$S2); output3<- max(samples$S3)
  pred.mean <- mean(c(output1,output2,output3))
  pred.sem <- sd(c(output1,output2,output3))/sqrt(2)
  
  # Function output is a vector containing these
  calculation <- c(pred.mean,pred.sem)
  return(calculation) # No code to handle if reaction is poor just yet
}

# Set up names for the data analysis
species.names <- c("S.lividans","S.albus","S.coelicolor","S.vioruber","S.albovinae","V.nat_OD1","V.nat_OD3")
cofactor.names <- c(0,2,3,4,5)

# Set up matrix for storage of the data
average.stored <- numeric(0)
sem.stored <- numeric(0)

# Loop through each cofactor and then their respective concentration
# Note: the enzyme is hard-coded here
for (r in 1:length(species.names)){
  for (t in 1:length(cofactor.names)){
    value <- MaxRFU(experiment.data,species.names[r],cofactor.names[t])
    average.stored[length(average.stored)+1] <- value[1]
    sem.stored[length(sem.stored)+1] <- value[2]
  } 
}

# Yield two matrices, first with the average reaction velocity and the second with the SEM
average.matrix <- data.frame(matrix(average.stored, nrow = 7, ncol = 5, byrow = TRUE, 
                                    dimnames = list(c(species.names), c(cofactor.names)))); average.matrix
sem.matrix <- data.frame(matrix(sem.stored, nrow = 7, ncol = 5, byrow = FALSE, 
                                dimnames = list(c(species.names), c(cofactor.names)))); sem.matrix





