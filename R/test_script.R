#' Create PCA Plot
#'
#' @param x csv file
#'
#' @return PCA plot
#'
#' @examples
#' pcaplot("csvfile.csv")
#'
#' @export

pcaplot <- function(csvfile) {

  #Input csv file
  #test <- read.csv("Assignment-2_gene_data-.csv")
  test <- read.csv(csvfile)

  #Removing the column which has the serial numbers (data specific)
  test <- test[, -1]
  test2 <- test[, -1]

  #"unique = TRUE" accounts for genes with same names but different values
  rownames(test2) <- make.names(test[, 1], unique = TRUE)

  dimensions.test2 <- dim(test2)

  #Converts non-numeric elements to NA
  test3 <- as.numeric(unlist(test2))

  #Restores dimensions as the original matrix
  dim(test3) <- dim(test2)

  #Replace the NA values with column mean
  for (i in 1:ncol(test3)) {
    test3[is.na(test3[, i]), i] <- mean(test3[, i], na.rm = TRUE)
  }

  #Assigns row names and column names same as the original matrix
  rownames(test3) <- test[, 1]
  colnames(test3) <- paste("S", 1:30, sep = "")

  #Computes singular value decomposition of a matrix to get the PCA plot
  test3.svd <- svd(scale(t(test3), center = TRUE))
  test3.pca <- data.frame(Sample = colnames(test3), X = (test3.svd$u[, 1] * test3.svd$d[1]), Y = (test3.svd$u[, 2] * test3.svd$d[2]))

  #Calculates percentage variances for a scree plot
  test3.svd.var <- test3.svd$d^2 / ncol(test3)
  test3.svd.var.perc <- round(test3.svd.var/sum(test3.svd.var)*100, 2)

  barplot(test3.svd.var.perc, main = "Scree Plot", xlab = "Principal Components", ylab = "Percentage Variation")

  library (ggplot2)

  #Plot of PC1 vs PC2
  ggplot(data = test3.pca, aes(x = X, y = Y, label = Sample)) +
    geom_text() +
    xlab(paste("PC1")) + ylab(paste("PC2")) +
    ggtitle("PCA")

}

#main()
