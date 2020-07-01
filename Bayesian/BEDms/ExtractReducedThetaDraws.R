
################################# EXTRACT REDUCED DRAWS #################################

# Script to extract all the draws for a determined stanfit inference result for one of the 3 models used in the study and
# from these, extract randomly 400 draws to be used for the optimisation procedure (more can be used, even all of them, but
# this come with a notable increase of computational time)


# MODEL 1

# The name will deppend on the users results name
exp2 <- "ALL_Model1.stan"

draws <- matrix(data = NaN, nrow = 8000, ncol = 16)

# Extract all the draws for the parameters
x <- as.array(readRDS(paste("fit_", exp2, ".rds", sep="")))

draws[,1] <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
draws[,2] <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
draws[,3] <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
draws[,4] <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
draws[,5] <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
draws[,6] <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
draws[,7] <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
draws[,8] <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
draws[,9] <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
draws[,10] <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
draws[,11] <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
draws[,12] <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])
draws[,13] <- c(x[,1,29], x[,2,29], x[,3,29], x[,4,29])
draws[,14] <- c(x[,1,30], x[,2,30], x[,3,30], x[,4,30])
draws[,15] <- c(x[,1,31], x[,2,31], x[,3,31], x[,4,31])
draws[,16] <- c(x[,1,32], x[,2,32], x[,3,32], x[,4,32])

# Save results
fn <- paste("draws_", exp2, ".csv", sep = "")
write.csv(draws, file=fn, row.names=FALSE)

# Extract 400 randomly independent draws
g = sample(1:8000, 400, replace=F)

drawsRed <- matrix(data = NaN, nrow = 400, ncol = 16)

for(i in 1:400){
  drawsRed[i,] = draws[g[i],]
}

# Save results
fn <- paste("drawsRedT_", exp2, ".csv", sep = "")
write.csv(drawsRed, file=fn, row.names=FALSE)



# MODEL 2

exp <- "ALL_Model2.stan"

draws <- matrix(data = NaN, nrow = 8000, ncol = 14)

x <- as.array(readRDS(paste("fit_", exp, ".rds", sep="")))

draws[,1] <- c(x[,1,15], x[,2,15], x[,3,15], x[,4,15])
draws[,2] <- c(x[,1,16], x[,2,16], x[,3,16], x[,4,16])
draws[,3] <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
draws[,4] <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
draws[,5] <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
draws[,6] <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
draws[,7] <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
draws[,8] <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
draws[,9] <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
draws[,10] <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
draws[,11] <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
draws[,12] <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
draws[,13] <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
draws[,14] <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])

fn <- paste("draws_", exp, ".csv", sep = "")
write.csv(draws, file=fn, row.names=FALSE)

g = sample(1:8000, 400, replace=F)

drawsRed <- matrix(data = NaN, nrow = 400, ncol = 14)

for(i in 1:400){
  drawsRed[i,] = draws[g[i],]
}


fn <- paste("drawsRedT_", exp, ".csv", sep = "")
write.csv(drawsRed, file=fn, row.names=FALSE)


# MODEL 2


exp <- "ALL_Long_Model3.stan"

draws <- matrix(data = NaN, nrow = 8000, ncol = 14)

x <- as.array(readRDS(paste("fit_", exp, ".rds", sep="")))

draws[,1] <- c(x[,1,15], x[,2,15], x[,3,15], x[,4,15])
draws[,2] <- c(x[,1,16], x[,2,16], x[,3,16], x[,4,16])
draws[,3] <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
draws[,4] <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
draws[,5] <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
draws[,6] <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
draws[,7] <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
draws[,8] <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
draws[,9] <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
draws[,10] <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
draws[,11] <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
draws[,12] <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
draws[,13] <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
draws[,14] <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])

fn <- paste("draws_", exp, ".csv", sep = "")
write.csv(draws, file=fn, row.names=FALSE)

g = sample(1:8000, 400, replace=F)

drawsRed <- matrix(data = NaN, nrow = 400, ncol = 14)

for(i in 1:400){
  drawsRed[i,] = draws[g[i],]
}


fn <- paste("drawsRedT_", exp, ".csv", sep = "")
write.csv(drawsRed, file=fn, row.names=FALSE)
