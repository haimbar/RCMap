
#' Create adjacency matrices from pile-sorting data
#'
#' N stakeholder sort S cards into piles, and this function converts this
#' data into n SxS adjacency matrices (one per sorter) in which cell i,j
#' contains 1 if user k put cards i and j in the same pile, and 0 otherwise.
#' The first column in the input file contains the sorters IDs, the
#' second column contains a name for the pile, given by the sorter
#' (this column may contain blanks), and the column 3,4,... contain
#' card numbers. If a user x created M piles, the file will contain
#' M rows with x in the first column.
#' Note that a card cannot be in more than one pile.
#' Sorters cannot put all the cards in one pile, but they can leave
#' cards unsorted.
#' @param piledat The pile sorting data.
#' @param showWarnings Print any potential problems in the pile-sorting data
#'  (default=TRUE).
#' @return A list of n SxS 0/1-matrices - one for each sorter.
#' @export
getAdjMatrices <- function(piledat, cardNames, sorters_dict, showWarnings=TRUE) {
  issues <- ""
  if (ncol(piledat) < 3)
    issues <- issues %+% red("Pile sorting file must have at least 3 columns!\n")
  sorters <- unique(piledat[ ,1])
  cardsused <- na.exclude(unique(as.numeric(stack(
    piledat[ ,3:ncol(piledat)])[ ,1])))
  cardsused <- sort(as.numeric(setdiff(cardsused, "")))
  nCards <- length(cardsused)
  mat.list <- list()
  #cat(sorters_dict$orig_ID,"\n")
  # Making a matrix for each individual and storing them in a list
  for (i in 1:length(sorters)) {
    sorterID <- which(sorters_dict$seq_ID == i)
    #cat(i, sorterID, sorters_dict$orig_ID[sorterID], "\n")
    #sorters_dict <- data.frame(orig_ID=sorters, seq_ID=unique(tmpids))
    sorter.dat <- piledat[which(piledat[ ,1] == sorters_dict$orig_ID[sorterID]), ]
    #                                sorters[i]), ]
    #cat(dim(sorter.dat),"\n")
    adj.mat <- matrix(0, nrow=nCards, ncol=nCards)
    for (j in 1:nrow(sorter.dat)) {
      cards <- sorter.dat[j, which(!is.na(sorter.dat[j, ]))]
      cards <- na.omit(as.numeric(sort(unlist(cards[-(1:2)]))))
      if (length(cards) == 0)
        next
      vct   <- rep(0, nCards)
      vct[which(cardNames[ ,1] %in% cards)] <- 1
      adj.mat <- adj.mat + vct%*%t(vct)
    }
    rsum <- rowSums(adj.mat)
    if(any(rsum == 0) & showWarnings) {
      didnotsort <- paste(cardNames[which(rsum == 0),1], collapse=",")
      issues <- issues %+% yellow(sprintf("Sorter  %s did not sort card(s) %s\n",
                                          sorters[i], didnotsort))
    }
    diag(adj.mat) <- 0
    if ((sum(adj.mat) == nCards^2) & showWarnings)
      issues <- issues %+% red("Sorter ",sorters[i], ": All cards in one pile!\n")
    if ((max(rowSums(adj.mat)) >= nCards/3) & showWarnings)
      issues <- issues %+% red("Sorter ",sorters[i], " put more than a third of the cards in one pile.\n")
    if ((sum(adj.mat) == nCards) & showWarnings)
      issues <- issues %+% red("Sorter ",sorters[i], ": Each card in its own pile!\n")
    # notused <- which(apply(adj.mat, 1, max) == 0)
    # if ((length(notused) > 0) & showWarnings)
    #   issues <- issues %+% yellow("Sorter ",i, " did not sort card(s) ", notused, "\n")
    gt1 <- which(apply(adj.mat, 1, max) > 1)
    #if(length(gt1) > 0) { cat(sorters[i], "\n", gt1,"\n\n", issues,"\n") }
    if (length(gt1) > 0) # this should not be allowed:
      issues <- issues %+% red("Sorter ",sorters[i], " put cards ", gt1, "in multiple piles\n")
    mat.list[[i]] <- adj.mat
  }
  cat(issues)
  return(list(mat.list=mat.list, issues=issues))
}


#' Add adjacency matrices from pile-sorting data and return a distance matrix
#'
#' n SxS adjacency matrices are added to create a similarity matrix, which is
#' then converted into a distance matrix in the S-dimensional Euclidean space.
#' @return An SxS distance matrix.
#' @export
distanceMatrix <- function(adjMats) {
  M <- Reduce("+", adjMats)
  n <- length(adjMats)
  S <- nrow(M)
  Imat <- diag(1, S, S)
  Jmat <- matrix(1, S, S)
  set.seed(154204)
  P <- n*(Jmat - Imat) - M
  P <- abs(jitter(P, 0.001)) #in case there are cards with identical pile distributions
  diag(P) <- 0
  return(P)
}


#' Calculate the sum of squares of a vector
#'
#' @param x a numeric vector.
#' @return The sum of squares.
sumsq <- function(x) {
  return(sum(x^2))
}


#' Convert Euclidean distances in the unit disk to a hyperbolic space distance.
#'
#' Converts a Euclidean distance matrix between points in a unit 2-D disk into a
#' (Poincare) distance in a 2-D hyperbolic space.
#' @param x distance matrix obtained for points in 2-D.
#' @return The hyperbolic distances.
#' @export
diskDist <- function(x) {
  sqnorms <- apply(x, 1, sumsq)
  P <- nrow(x)
  dd <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in 1:P) {
      dd[i, j] <- dd[j, i] <- acosh(1 + 2*(sumsq(x[i, ] - x[j, ]))/
                                      ((1 - sqnorms[i]) * (1 - sqnorms[j])))
    }
  }
  return(dd)
}

get_cmat <- function(dists, M, disttype, max_nc, clustMethod="ward.D2") {
  groups <- Cmat <- list()
  for (k in 2:max_nc) {
    fit.Clust <- hclust(as.dist(dists), method=clustMethod)
    Cmat[[k-1]] <- matrix(0, M, M)

    groups[[k-1]] <- cutree(fit.Clust, k=k)
    for (i in 1:k) {
      jj <- which(groups[[k-1]] == i)
      Cmat[[k-1]][jj, jj] <- 1
    }
    diag(Cmat[[k-1]]) <- 0
  }
  return(Cmat)
}

clusterings <- function(adjmat, j=0, max_nc=3, clustMethod="ward.D2") {
  if (j > length(adjmat)) {
    #cat("Parameter (j) is greater the the length of list (adjmat). Using 0\n")
    j <- 0
  }
  A <- adjmat
  if (j > 0) {
    A <- adjmat[-j]
  }
  D1 <- distanceMatrix(A)
  M <- nrow(D1)
  fit.MDS1 <- mds(D1)
  d1 <- as.matrix(dist(fit.MDS1$conf))
  d2 <- diskDist(d1)

  CmatE <- get_cmat(d1, M, "Euclidean", max_nc, clustMethod)
  CmatH <- get_cmat(d2, M, "Hyperbolic", max_nc, clustMethod)
  list(CmatE=CmatE, CmatH=CmatH) # lists of adjacency matrices, for different numbers of clusters
}


#' Calculate a misplacement index.
#'
#' Takes a dissimilarity matrix from a pile-sorting study and a distance matrix
#' obtained from the 2-D map (from multidimensional scaling) and returns an
#' index which is close to 0 for statements which are positioned close (in 2-D)
#' to other statements with which they put in the same pile often (anchors),
#' and a value close to 1 for points which are positioned in the 2-D map
#' far from other  statements with which were put in the same pile (bridges).
#' @param D dissimilarity matrix obtained for the pile sorting data.
#' @param d distance matrix obtained for points in the 2-D projection.
#' @param cst The power parameter for the misplacement index function (default=1).
#' @param func the function applied to the pair of matrices, d and D. Default=median.
#' @param ... additional parameters to be passed to func.
#' @return A vector of size nrow(D) with the bridging/anchoring metric.
#' @export
MPindex <- function(D, d, cst=1, func=median, ...) {
  d <- d/max(d)
  D <- D/max(D)
  M0 <- (abs(D - d))^(1 + cst*pmin(d, D))
  M <- matrix(0, nrow = nrow(M0), ncol = ncol(M0)-1)
  for (i in 1:nrow(M)) { M[i, ] <- M0[i, -i] }
  return(apply(M, 1, func, ...))
}


#' Calculate reliability metrics.
#'
#' Takes a dissimilarity matrix from a pile-sorting study and a distance matrix
#' obtained from the 2-D map (from multidimensional scaling) and returns an
#' index which is close to 0 for statements which are positioned close (in 2-D)
#' to other statements with which they put in the same pile often (anchors),
#' and a value close to 1 for points which are positioned in the 2-D map
#' far from other  statements with which were put in the same pile (bridges).
#' @param adjMat pile-sorting data (n 0-1 matrices).
#' @param B the number of random splits (default=10).
#' @param disttype The distance metric to be used ("Hyperbolic", otherwise it's assumed to be Euclidean).
#' @param seed The random seed to be used.
#' @param plotit if TRUE, show a plot of the differences between the misplacement index between the two random halves.
#' @return A list with two elements: reliability, which contains the mean correlation between the splits, the distance metric being used, and the function used to calculate the misplacement index.
#' @export
splitHalf <- function(adjMat, B=10, disttype="Hyperbolic",
                      seed=968421, plotit=FALSE) {
  n <- length(adjMat)
  cors <- rep(0, B)
  for (i in 1:B) {
    set.seed(seed + i)
    grp1 <- sort(sample(n)[1:round(n/2)])
    grp2 <- setdiff(1:n, grp1)
    adjMattmp1 <- adjMat[grp1]
    D1 <- distanceMatrix(adjMattmp1)
    adjMattmp2 <- adjMat[grp2]
    D2 <- distanceMatrix(adjMattmp2)
    fit.MDS1 <- mds(D1)
    fit.MDS2 <- mds(D2)
    d1 <- as.matrix(dist(fit.MDS1$conf))
    d2 <- as.matrix(dist(fit.MDS2$conf))
    distmetric <- "Euclidean"
    if (disttype == "Hyperbolic") {
      d1 <- diskDist(d1)
      d2 <- diskDist(d2)
      distmetric <- "Hyperbolic"
    }
    yrng <- c(0, 1)
    cors[i] <- cor.test(as.vector(d1), as.vector(d2))$estimate
  }
  if(plotit) {
    plot(cors, cex=0.6, col="blue", main="Split-half correlations", pch=19,
         xlab="Split #", ylab="Correlation between halves", ylim=c(-1,1),
         axes=F)
    axis(1); axis(2); abline(h=0, col="grey", lwd=2); grid()
  }
  list(cors=cors, distmetric=distmetric)
}


#' Perform a leave-one-out analysis on the sorters and plot the MADD values (mean absolute distance distortion) in a dotchart.
#'
#' Takes a list of Jaccard-index matrices as input creates an interaction plot, showing the scores for all the statements, as well as the median, and the 25th and 75th percentiles.
#' @param eta_ij A list of adjacency matrix from all the sorters.
#' @export
plotjaccard <- function(eta_ij) {
  M <- matrix(0, nrow=length(eta_ij), ncol=ncol(eta_ij[[1]]))
  M2 <- matrix(0, nrow=0, ncol=3)
  for (i in 1:length(eta_ij)) {
    M[i, ] <- colMeans(eta_ij[[i]] < 0.3)
    M2 <- rbind(M2, cbind(rep(i+1, length(colMeans(eta_ij[[i]]))),
                          1:ncol(eta_ij[[1]]), colMeans(eta_ij[[i]])))
  }
  M2 <- as.data.frame(M2)
  colnames(M2) <- c("clusters", "statement", "eta")
  #interaction.plot(M2$statement, M2$clusters, M2$eta, legend=F, lwd=0.2)
  interaction.plot(M2$clusters, M2$statement, M2$eta, legend=F, lwd=0.2,
                   xlab="No. clusters", ylab="Jaccard", col="grey66",
                   axes=FALSE, main="")
  clusterstats <- by(M2$eta, M2$clusters, quantile, probs=c(0.25, 0.5, 0.75))
  clusterstats <- array2DF(clusterstats, simplify = TRUE)
  lines(unique(M2$clusters)-1, clusterstats$`25%`, col="blue", lwd=2)
  lines(unique(M2$clusters)-1, clusterstats$`50%`, col="navyblue", lwd=3)
  lines(unique(M2$clusters)-1, clusterstats$`75%`, col="blue", lwd=2)
  axis(1, at=seq(1, max(M2$clusters)-1), labels=seq(2, max(M2$clusters)))
  axis(2, at=seq(0.1, 1, by=0.1))
}


# Check that the file containing the statements is loaded properly:
read_statements <- function(filename, enc = "UTF-8", sep=",") {
  if(!file.exists(filename))
    return(list(issues=paste("read_statements: File", filename,
                             "doesn't exist.\n", troubleshooting["cardfilename"]),
                cardNames = NULL))
  cardNames <- read.csv(filename, header=TRUE,  encoding = enc, sep=sep)
  # The file must have exactly two columns: ID and statement text:
  if (ncol(cardNames) != 2)
    return(list(issues=paste("read_statements:", filename,
                             "must have exactly two columns.\n",
                             troubleshooting["cardfilesep"]),
                cardNames = NULL))
  # The file must have some data (we're only checking if it has at least 3 rows)
  if (nrow(cardNames) < 3)
    return(list(issues=paste("read_statements:", filename,
                             "must have some statements."),
                cardNames = NULL))
  # File structure appears to be correct (encoding is left for the user to set.)
  # Card numbers should be consecutive in the input file, but we're going to
  # assume that it may not be the case, and create a column called SeqID:
  colnames(cardNames) <- c("StatementID", "Statement")
  cardNames$SeqID <- 1:nrow(cardNames)
  cardNames$LC <- tolower(iconv(cardNames$Statement, from = "ISO-8859-1",
                                to = "UTF-8"))
  cardNames$LC <- gsub("^\\s+(.+)\\s+$", "\\1", cardNames$LC)
  return(list(issues=NULL, cardNames=cardNames))
}


# Check that the file containing the pile-sorting data is loaded properly:
read_piles <- function(filename, ncards, enc = "UTF-8", sep=",") {
  if(!file.exists(filename))
    return(list(issues=paste("read_piles: File", filename,
                             "doesn't exist.\n", troubleshooting["pilesfilename"]),
                cardDat = NULL))

  cardDat <- suppressWarnings(read.csv(filename, header = FALSE,
                                       encoding = enc, sep = sep,
                                       colClasses = rep("character", 2+ncards)))
  if (ncol(cardDat) <= 3)
    return(list(issues=paste("read_piles:", filename,
                             "must have at least three columns.\n",
                             troubleshooting["pilesfilesep"]),
                cardDat = NULL))

  cardDat <- cardDat[ ,which(colnames(cardDat) != "")]
  pileLabels <- data.frame(pileLabel=character(), cardNum=integer())
  for (i in 1:nrow(cardDat)) {
    cardDat[i,] <- iconv(cardDat[i,], from = "ISO-8859-1", to = "UTF-8")
    cardDat[i,] <- gsub("\\s+$", "", cardDat[i,])
    pile <- cardDat[i, ]
    pileLabel <- cardDat[i, 2]
    if(length(grep("\\D", cardDat[i, -c(1,2)])) > 0)
      return(list(issues=paste("read_piles:", filename, "Line ", i,
                               troubleshooting["non-numeric"]),
                  cardDat = NULL))
    cardsInPile <- na.omit(as.numeric(cardDat[i, -c(1,2)]))
    if (length(cardsInPile) == 0)
      next
    cardNum <- unique(sort(as.numeric(cardsInPile)))
    pileLabel <- rep(pileLabel, length(cardsInPile))
    pileLabels <- rbind(pileLabels, data.frame(pileLabel, cardNum=cardNum))
  }
  #  cardDat[ ,1] <- gsub("\\D", "", cardDat[ ,1])
  pileLabelsLC <- tolower(iconv(pileLabels[ ,1], from = "ISO-8859-1",
                                to = "UTF-8"))
  pileLabelsLC <- gsub("^\\s+(.+)\\s+$", "\\1", pileLabelsLC)
  return(list(issues=NULL, cardDat=cardDat, pileLabels=pileLabels,
              pileLabelsLC=pileLabelsLC))
}


# Check that the file containing the rating data is loaded properly:
read_ratings <- function(filename, statements, enc = "UTF-8", sep=",") {
  if(!file.exists(filename))
    return(list(issues=paste("read_ratings: File", filename,
                             "doesn't exist.\n", troubleshooting["ratingsfilename"]),
                cardDat = NULL))
  ratings <- read.csv(filename, header=TRUE)
  if (ncol(ratings) <= 2)
    return(list(issues=paste("read_ratings:", filename,
                             "must have at least three columns.\n",
                             troubleshooting["ratingsfilesep"]),
                ratings = NULL))
  notinstatements <- which(ratings[,2] %in% statements == FALSE)
  if (length(notinstatements) > 0)
    return(list(issues=paste("read_ratings:", filename,
                             "Some rated statements are not in the Statements.csv file:\n",
                             notinstatements),
                ratings = NULL))
  if (!(colnames(ratings)[1] %in% c("RaterID", "UserID", "SorterID")))
    return(list(issues=paste("read_ratings:", filename,
                             "The first column name must be RaterID\n"),
                ratings = NULL))
  colnames(ratings)[1] <- "RaterID"
  if (colnames(ratings)[2] != "StatementID")
    return(list(issues=paste("read_ratings:", filename,
                             "The second column names must be StatementID\n"),
                ratings = NULL))

  tmpmat <- as.matrix(ratings[ ,3:ncol(ratings)], nrow=nrow(ratings))
  nonlikert <- cbind(rowSums(tmpmat < 1), rowSums(tmpmat > 5))
  nonlikert <- colSums(nonlikert, na.rm = TRUE)
  if (max(nonlikert) > 0)
    return(list(issues=paste("read_ratings:", filename,
                             "The rating values must be between 1 and 5\n",
                             "Columns [", colnames(ratings)[2+which(nonlikert > 0)],
                             "] have invalid values.\n"),
                ratings = NULL))
  ratings[,2] <- as.factor(ratings[,2])
  return(list(issues=NULL, ratings=ratings))
}

# Check that the file containing the demographic data is loaded properly:
read_demographics <- function(filename, ratings, enc = "UTF-8", sep=",") {
  if(!file.exists(filename))
    return(list(issues=paste("read_demographics: File", filename,
                             "doesn't exist.\n", troubleshooting["demographicsfilename"]),
                demographics.dat = NULL))
  demographics.dat <- read.csv(filename, stringsAsFactors = TRUE,
                               header=TRUE)
  if (ncol(demographics.dat) < 2)
    return(list(issues=paste("read_demographics:", filename,
                             "must have at least two columns.\n",
                             troubleshooting["demographicsfilesep"]),
                demographics.dat = NULL))
  if (!(colnames(demographics.dat)[1] %in% c("RaterID", "UserID", "SorterID")))
    return(list(issues=paste("read_demographics:", filename,
                             "The first column name must be RaterID\n"),
                ratings = NULL))
  colnames(demographics.dat)[1] <- "RaterID"
  dframe <- data.frame(matrix(0, ncol=ncol(demographics.dat)-1, nrow=nrow(ratings)))
  cohorts <- data.frame("allRaters"=matrix(1, nrow=nrow(ratings), ncol=1))
  skipcols <- c()
  for(colnum in 2:ncol(demographics.dat)) {
    if(length(unique(demographics.dat[,colnum])) == 1) { # no variation in the column
      skipcols <- c(skipcols, colnum)
      next
    }
    if (class(demographics.dat[,colnum]) == "factor") {
      if(length(unique(demographics.dat[,colnum])) == nrow(demographics.dat)) { # all distinct
        skipcols <- c(skipcols, colnum)
        next
      }
      dframe[,colnum-1] <- factor(rep("", nrow(ratings)),
                                  levels=levels(demographics.dat[,colnum]))
    } else {
      dframe[,colnum-1] <- rep(0, nrow(ratings))
    }
    for (i in 1:nrow(ratings)) {
      matched <- which(demographics.dat$RaterID == ratings$RaterID[i])
      if (length(matched) > 1)
        return(list(issues=paste("read_demographics:", filename,
                                 "Rater ", ratings$RaterID[i],
                                 "appears more than once.\n"),
                    demographics.dat = NULL))
      if (length(matched) == 0)
        return(list(issues=paste("read_demographics:", filename,
                                 "Rater ", ratings$RaterID[i],
                                 "has no demographic data.\n"),
                    demographics.dat = NULL))

      dframe[i,colnum-1] <- demographics.dat[,colnum][which(demographics.dat$RaterID ==
                                                              ratings$RaterID[i])]
    }
  }
  if (length(skipcols) == ncol(dframe)) {
    return(list(issues=paste("read_demographics: File", filename,
                             "Invalid demographic data.\n", troubleshooting["demographicsnodata"]),
                demographics.dat = NULL))
  }
  if (length(skipcols) > 0) {
    dframe <- as.data.frame(dframe[, setdiff(1:ncol(dframe), skipcols-1)])
    colnames(dframe) <- colnames(demographics.dat)[-c(1,skipcols)]
  }
  for (colnum in 1:ncol(dframe)) {
    if(class(dframe[,colnum]) == "factor") {
      tmp <- model.matrix( ~ dframe[,colnum] - 1)
      colnames(tmp) <- sprintf("%s_%s",colnames(dframe)[colnum],
                               levels(dframe[,colnum]))
      cohorts <- data.frame(cohorts, tmp)
    } else { # type=numeric is assumed
      tmp <- model.matrix( ~ (dframe[,colnum] > median(dframe[,colnum])) - 1)
      colnames(tmp) <- paste0(colnames(dframe)[colnum],c("_H","_L"))
      cohorts <- data.frame(cohorts, tmp)
    }
  }
  return(list(issues=NULL, demographics.dat=demographics.dat,
              dframe=dframe, cohorts=cohorts))
}

#' Read the card-sorting data files and generate the MDS data.
#'
#' Takes a data directory as input. This directory must contain the four project files.
#' @param dataDir The project's directory.
#' @param enc Encoding (default is UTF-8).
#' @param sep Separator for read.csv (default is comma).
#' @return A list with 9 elements n.indiv=number of sorters, cardNames=the statements on the cards, adj.mat=the adjacency matrix, DS=the pairwise distances matrix in the original high-dimensional space, D2e=the Euclidean distance matrix in the 2-D space, D2h=the hyperbolic distance matrix in the 2-D space, x=the x coordinates in the MDS plot, y=the y coordinates, stress=the stress of the MDS plot.
#' @importFrom smacof mds
#' @export
initCMap <- function(dataDir, enc="UTF-8", sep=",") {
  # Check that the four input files are all there, and are formatted correctly:

  clustMethod <<- "ward.D2"
  savedfile <- paste0(dataDir,"/CMapSession.RData")
  if (file.exists(savedfile)) {
    cat(green("A previous project file was found.\n"))
    prmpt <- "Hit Enter to use saved project, or enter L to load the project from raw files: "
    usefile <- readline(prompt = prmpt)
    if (toupper(usefile) == "") {
      load(savedfile)
      return(cmapdat)
    }
  }
  cat(green("Loading project files. Please wait.\n\n"))
  cmapdat <<- list()
  # Statements:
  cardNames <- read_statements(paste0(dataDir, "/Statements.csv"), enc, sep)
  if(length(cardNames$issues) > 0)
    stop(cardNames$issues)
  # If no issues - keep only the data frame
  cardNames <- cardNames$cardNames
  card_name_dict <- cardNames

  # Piles:
  cardDat <- read_piles(paste0(dataDir, "/SortedCards.csv"), nrow(cardNames),
                        enc, sep)
  if(length(cardDat$issues) > 0)
    stop(cardDat$issues)
  # If no issues - keep only the data frame
  label_dict <- data.frame(pileLabels=cardDat$pileLabels[,1],
                           lcpilelabels=cardDat$pileLabelsLC)
  pileLabels <- cardDat$pileLabelsLC
  cardDat <- cardDat$cardDat

  sorters <- unique(cardDat[ ,1])
  tmpids <- rep(-1, length(sorters))
  for (i in 1:length(sorters)) {
    tmpids[which(cardDat[ ,1] == sorters[i])] <- i
  }
  sorters_dict <- data.frame(orig_ID=sorters, seq_ID=unique(tmpids))
  #cardDat[ ,1] <- tmpids

  # Create the statement adjacency matrix from the piles:
  tmp <- getAdjMatrices(cardDat, cardNames, sorters_dict)
  adj.mat <- tmp[[1]]
  issues <- tmp[[2]]
  n.indiv <- length(adj.mat)
  # multidimensional scaling, and clustering:
  DS <- distanceMatrix(adj.mat)
  rownames(DS) <- cardNames$StatementID
  fit.MDS <- mds(DS)
  x <- fit.MDS$conf[,1]
  y <- fit.MDS$conf[,2]
  stress <- fit.MDS$stress
  D2e <- as.matrix(dist(fit.MDS$conf))
  D2h <- diskDist(D2e)
  rownames(D2e) <- rownames(D2h) <- cardNames[ ,1]
  #cardNames <- cardNames$Statement

  # Ratings:
  ratings <- read_ratings(paste0(dataDir, "/Ratings.csv"),
                          cardNames$StatementID, enc, sep)
  if(length(ratings$issues) > 0)
    stop(ratings$issues)
  ratings <- ratings$ratings
  # Demographic data:
  demographics.dat <- read_demographics(paste0(dataDir, "/Demographics.csv"),
                                        ratings, enc, sep)
  if(length(demographics.dat$issues) > 0)
    stop(demographics.dat$issues)
  dframe <- demographics.dat$dframe
  cohorts <- demographics.dat$cohorts
  demographics.dat <- demographics.dat$demographics.dat

  if (!dir.exists(paste0(dataDir,"/output")))
    dir.create(paste0(dataDir,"/output"))

  M <- nrow(adj.mat[[1]]) # number of statements
  N <- length(adj.mat)    # number of sorters
  max_nc <- round(M/4)

  cat("\nNumber of sorters:", N, "\n")
  cat("Number of statements:", M, "\n")
  cat("Calculating the Jaccard index with number of clusters from 2 to", max_nc, "...\n\n")
  C0 <- clusterings(adj.mat, max_nc = max_nc, clustMethod)
  eta_ijk_e <- eta_ijk_h <- Cj <- list()
  for (m in 1:(max_nc-1)) {
    eta_ijk_e[[m]] <- eta_ijk_h[[m]] <- matrix(0, N, M)
  }
  for (j in 1:N) {
    Cj <- clusterings(adj.mat, j=j, max_nc = max_nc, clustMethod)
    c0tmp <- C0[[1]][[m]]
    cjtmp <- Cj[[1]][[m]]
    eta_tmp <- colSums(pmin(c0tmp, cjtmp))/colSums(pmax(c0tmp, cjtmp))
    eta_tmp[is.nan(eta_tmp)] <- 0
    for (m in 1:(max_nc-1)) {
      eta_ijk_e[[m]][j, ] <- eta_tmp
      c0tmp <- C0[[2]][[m]]
      cjtmp <- Cj[[2]][[m]]
      eta_tmp <- colSums(pmin(c0tmp, cjtmp))/colSums(pmax(c0tmp, cjtmp))
      eta_tmp[is.nan(eta_tmp)] <- 0
      eta_ijk_h[[m]][j, ] <- eta_tmp
    }
  }

  save(card_name_dict, sorters_dict, label_dict,
       file=paste0(dataDir,"/dictionaries.RData"))
  cmapdat <- list(n.indiv=n.indiv, cardNames=cardNames$Statement,
                  adj.mat=adj.mat, DS=DS,
                  D2e=D2e, D2h=D2h, x=x, y=y, stress=stress, ratings=ratings,
                  clustname=list(), pileLabels=pileLabels, issues=issues,
                  demographics=demographics.dat, ratingDemographics=dframe,
                  clustMethod=clustMethod,
                  distmetric="Euclidean", cohortCol=1,
                  nclust=3, splhalf=list(), cohorts=cohorts,
                  pcoordChoice=c(1, 2), clusMember=cmapdat$clusMember, #rep(1 ,length(cardNames$Statement)),
                  clrscm="rcmap", eta_ijk_e=eta_ijk_e, eta_ijk_h=eta_ijk_h,
                  dataDir=dataDir)
  save(cmapdat, file=paste0(dataDir,"/CMapSession.RData"))
  return(cmapdat)
}


#' Get ratings statistics (N, mean, sd, min, max) for each statement
#'
#' Raters provide Likert scale ratings for k criteria (e.g. feasibility). This function returns summary statistics for each statement.
#' @param ratingsDat The ratings data.
#' @return A matrix with summary statistics for all rating variables and all statements.
#' @export
ratingSummary <- function(ratingsDat) {
  nVar <- ncol(ratingsDat)-2
  ratingsDat$StatementID <- droplevels(ratingsDat$StatementID)
  summ <- matrix(0, nrow=length(unique(ratingsDat[,2]))+1, ncol=5*nVar)
  cnames <- rep("", nVar*5)
  for (colnum in 1:nVar) {
    summ[,(colnum-1)*5+ 1] <- c(rowSums(table(ratingsDat[,2], ratingsDat[,colnum+2])),
                                mean(rowSums(table(ratingsDat[,2], ratingsDat[,colnum+2]))))
    summ[,(colnum-1)*5+ 2] <- c(as.numeric(by(ratingsDat[,colnum+2], ratingsDat[,"StatementID"], mean, na.rm=TRUE)),
                                mean(ratingsDat[,colnum+2], na.rm = TRUE))
    summ[,(colnum-1)*5+ 3] <- c(as.numeric(by(ratingsDat[,colnum+2], ratingsDat[,"StatementID"], sd, na.rm=TRUE)),
                                sd(ratingsDat[,colnum+2], na.rm = TRUE))
    summ[,(colnum-1)*5+ 4] <- c(as.numeric(by(ratingsDat[,colnum+2], ratingsDat[,"StatementID"], min, na.rm=TRUE)),
                                min(ratingsDat[,colnum+2], na.rm = TRUE))
    summ[,(colnum-1)*5+ 5] <- c(as.numeric(by(ratingsDat[,colnum+2], ratingsDat[,"StatementID"], max, na.rm=TRUE)),
                                max(ratingsDat[,colnum+2], na.rm = TRUE))
    cnames[(colnum-1)*5 + 1:5] <- paste(colnames(ratingsDat)[colnum+2],
                                        c("N","Mean","SD","Min", "Max"), sep="_")
  }
  colnames(summ) <- cnames
  summ
}


#' Plot the 2-D representation of the pile-sorting data
#'
#' @param cols The colors of the data (default= all blue.)
#' @param metric The type of plot to show (NULL, the default is for the simple map.)
#' @param tau The threshold for the Jaccard index (0.3)
#' @export
showMDSPlot <- function(metric=NULL, tau=0.3) {
  ter.cols <- rev(c("#440154D0","#472D7BD0","#3B528BD0","#2C728ED0",
                    "#21908CD0","#27AD81B0","#5DC863D0","#AADC32F0","#FDE725F0"))
  if(!is.null(metric))
    if(metric == 0)
      return(NULL)
  if(is.null(metric)) {
    sz <- 0.25
    cols <- "blue"
    ttl <- ""
  } else {
    if (metric == "MPindex") {
      ttl <- "Misplacement index"
      if(cmapdat$distmetric == "Hyperbolic")
        sz <- colMeans(cmapdat$eta_ijk_h[[cmapdat$nclust-1]] < tau)
      if(cmapdat$distmetric == "Euclidean")
        sz <- colMeans(cmapdat$eta_ijk_e[[cmapdat$nclust-1]] < tau)
      lgdkey <- seq(0, 1, length=length(ter.cols))
    } else {
      ttl <- paste(colnames(cmapdat$ratings)[metric+2],
                   colnames(cmapdat$cohorts)[cmapdat$cohortCol])
      ratingstmp <- cmapdat$ratings[which(cmapdat$cohorts[,cmapdat$cohortCol] == 1),]
      rtgsmm <- ratingSummary(ratingstmp)
      rownames(rtgsmm) <- c(rownames(cmapdat$DS), "total")
      sz <- rtgsmm[-nrow(rtgsmm),
                   which(colnames(rtgsmm) ==
                           paste0(colnames(cmapdat$ratings)[metric+2],"_Mean"))]
      sz <- (sz-1)/4
      lgd <- seq(0, 1, length=length(ter.cols))
      lgdkey <- lgd*4 + 1
    }
    if (!is.null(metric)) {
      cols <- ter.cols[cut(sz, breaks = seq(0, 1, length=length(ter.cols)),
                           include.lowest = TRUE)]
      lgd <- seq(0, 1, length=length(ter.cols))
    }
  }
  plot(cmapdat$x, cmapdat$y, pch=19, col=cols,
       cex=sz+0.4, main=ttl, xlab="",
       ylab="", xlim=c(min(cmapdat$x)-0.25,max(cmapdat$x)*1.2),
       ylim=c(min(cmapdat$y),max(cmapdat$y)*1.2), axes=F, cex.main=1)
  text(cmapdat$x, cmapdat$y, pos=3, cex=0.8, labels=names(cmapdat$x))
  if(length(sz) > 1) {
    rect(min(cmapdat$x)-0.3, 0.2*max(cmapdat$y),
         min(cmapdat$x)-0.1, 1.3*max(cmapdat$y), col="whitesmoke", border="white")
    points(rep(min(cmapdat$x)-0.15, length(lgd)),
           seq(0.3*max(cmapdat$y), 1.2*max(cmapdat$y),length=length(lgd)),
           col=ter.cols, cex = lgd+0.4, pch=19)
    text(rep(min(cmapdat$x)-0.25, length(lgd)),
         seq(0.3*max(cmapdat$y), 1.2*max(cmapdat$y), length=length(lgd)),
         lgdkey, cex=0.6)
    if (metric == "MPindex") {
      x.cent <- rep(0, cmapdat$nclust)
      y.cent <- rep(0, cmapdat$nclust)
      for (i in 1:cmapdat$nclust)  {
        x.cent[i] <- mean(cmapdat$x[cmapdat$clusMember==i])
        y.cent[i] <- mean(cmapdat$y[cmapdat$clusMember==i])
        x.groups <- cmapdat$x[cmapdat$clusMember==i]
        y.groups <- cmapdat$y[cmapdat$clusMember==i]
        for (j in 1:length(x.groups)) {
          lines(c(x.cent[i], x.groups[j]), c(y.cent[i], y.groups[j]),
                lty=3, lwd=1, col="grey66")
        }
      }
    }
  }
}

#' Plot the clusters in the 2-D representation of the pile-sorting data
#'
#' Can show the clusters either as convex polygons, or as lines connected to the center of each cluster.
#' @param type Rays, or polygon.
#' @param metric cols The metric to be used (Euclidean or hyperbolic.)
#' @export
showClusterPlot <- function(type="rays", metric=NULL) {
  if(!is.null(metric))
    if(metric == 0)
      return(NULL)
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    sz <- rep(0, nclust)
    if(!is.null(metric)) {
      ttl <- paste(colnames(cmapdat$ratings)[metric+2],
                   colnames(cmapdat$cohorts)[cmapdat$cohortCol])
      ratingstmp <- ratings[which(cohorts[,cohortCol] == 1),]
      for (i in 1:length(unique(groups))) {
        cardsInCluster <- names(which(groups == i))
        cardsInCluster <- which(ratingstmp$StatementID %in% cardsInCluster)
        rtgsmm <- ratingSummary(ratingstmp[cardsInCluster,])
        rownames(rtgsmm) <- c(rownames(DS)[which(groups == i)], "total")
        sz[i] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_Mean"))]
      }
      sz <- 1+2*(sz-1)/4
      sz <- round(sz, digits = 1) # to be between 1 (for 1) and 2 (for 5)
      lgd <- seq(1, 3, by=0.5)
      lgdkey <- (lgd-1)*2+1
    }

    colpoints <- cols[groups]
    if (type=="polygon") {
      eps <- matrix(c(0,0.02,-0.02,0,0,-0.02,0.02,0),ncol=2,byrow = T)
      eps <- rbind(eps, eps%*%matrix(c(1,1,-1,1),ncol=2)*sqrt(2)/2)
      colpoints <- coltext <- "white"
    }
    plot(x,y,pch=15,col=colpoints,cex=0.6,main=ttl,xlab="",ylab="",
         xlim=c(min(x)-0.25,max(x)*1.2),
         ylim=c(min(y),max(y)*1.2),axes=F)
    x.cent <- rep(0, nclust)
    y.cent <- rep(0, nclust)
    for (i in 1:nclust)  {
      x.cent[i] <- mean(x[groups==i])
      y.cent[i] <- mean(y[groups==i])
      x.groups <- x[groups==i]
      y.groups <- y[groups==i]
      if (type == "rays") {
        for (j in 1:length(x.groups)) {
          lines(c(x.cent[i], x.groups[j]), c(y.cent[i], y.groups[j]), lty=3,
                lwd=2, col=cols[i])
        }
      } else {
        ptsx = x[which(groups==i)]; ptsy = y[which(groups==i)]
        hpts = chull(ptsx, ptsy)
        hpts <- c(hpts, hpts[1])
        extpts <- kronecker(cbind(ptsx[hpts],ptsy[hpts]),rep(1,8))+
          kronecker(rep(1,length(hpts)),eps)
        hpts <- chull(extpts)
        hpts <- c(hpts, hpts[1])
        polygon(extpts[hpts,], col=cols[i], border="white")
        points(x, y, pch=15, col=colpoints, cex=0.6)
      }
    }
    for (i in 1:nclust)  {
      points(x.cent[i], y.cent[i],pch=18,col="orange", cex=sz[i])
    }
    for (i in 1:length(cardNames)) {
      text(x[i], y[i], cex=0.7, groups[i]+2, labels=names(cmapdat$x)[i],
           pos=3)
    }
    if(!is.null(metric)) {
      rect(min(x)-0.3, 0.2*max(y),
           min(x)-0.1, 1.3*max(y), col="whitesmoke", border="white")
      points(rep(min(x)-0.15, length(lgd)),
             seq(0.3*max(y), 1.2*max(y),length=length(lgd)),
             col="orange", cex = lgd, pch=18)
      text(rep(min(x)-0.25, length(lgd)),
           seq(0.3*max(y), 1.2*max(y),length=length(lgd)),
           lgdkey,cex=0.6)
    }
    # Plot Cluster Centers
    lbl = clusterNames()
    text(x.cent,y.cent,cex=0.8,labels=lbl, font=3,col="blue")
  })
}


update_clusters <- function() {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    clus <- cutree(fit.Clust, nclust)
    cmapdat$clusMember <<- clus
    cat(cmapdat$clusMember,"\n")
  })
}

#' Plot the clusters as a dendrogram or a phylogenic tree plot.
#'
#' @param phylo Whether to show a dendrogram (a default) or a phylogenic tree.
#' @importFrom ape as.phylo
#' @export
showDendrogram <- function(phylo=FALSE) {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    clus <- cutree(fit.Clust, nclust)
    cmapdat$clusMember <<- clus
    if(phylo) {
      X <- as.phylo(fit.Clust)
      edge.clus <- sapply(1:cmapdat$nclust,function(i)max(which(X$edge[,2] %in% which(clus==i))))
      order     <- order(edge.clus)
      edge.clus <- c(min(edge.clus),diff(sort(edge.clus)))
      edge.clus <- rep(order,edge.clus)
      plot(X,type='fan',
           tip.color=cols[clus],edge.color=cols[edge.clus],
           label.offset=0,no.margin=FALSE, cex=0.70)
      text(rep(0.6,length(cardNames)), max(fit.Clust$height)*
             seq(-0.5, 0.5, length=length(clusterNames())), clusterNames(),
           pos=4, col=cols, cex=0.8, font=2)
    } else {
      clusDendro <- dendrapply(as.dendrogram(fit.Clust), colLab)
      plot(clusDendro, main = "")
      text(rep(2,length(cardNames)),  max(fit.Clust$height)*
             seq(0.3, 1, length=length(clusterNames())), clusterNames(),
           pos=4, col=cols, cex=0.8, font=2)
    }
  })
}

#' Show statement ratings as a dotchart.
#'
#' @param metric Euclidean or hyperbolic.
#' @export
showDotPlot <- function(metric=NULL) {
  if(metric == 0)
    return(NULL)
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 2) - 1) - 2
    y <- c()
    nm <- c()
    cls <- c()
    if(!is.null(metric)) {
      ttl <- paste(colnames(cmapdat$ratings)[metric+2],
                   colnames(cmapdat$cohorts)[cmapdat$cohortCol])
      ratingstmp <- ratings[which(cohorts[,cohortCol] == 1),]
      rtgsmm <- ratingSummary(ratingstmp)
      rownames(rtgsmm) <- c(rownames(DS), "total")
      y <- rtgsmm[-nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_Mean"))]
      #layout(matrix(c(1,2), 1, 2, byrow = TRUE))#,widths=c(2,1))
      dotchart(sort(y, decreasing = FALSE),
               labels = row.names(rtgsmm[order(y, decreasing = FALSE)]),
               groups = groups[order(y, decreasing = FALSE)], gcolor = cols,
               color = cols[groups[order(y, decreasing = FALSE)]], xlim = c(1,6),
               cex = 0.6,  pch = 19, xlab = ttl, frame.plot = FALSE, xaxt = "n")
      abline(v=1:4, lwd=2, lty=2, col="gray66")
      text((1:5), rep(max(gpos)+1, 5), 1:5, cex=0.8, font=2)
      rect(5.05,0,7,max(gpos)+2, col = "white", border="white")
      rect(0,0,0.95,max(gpos)+2, col = "white", border="white")
      text(rep(5.2, length(cardNames)), gpos, clusterNames(), cex=0.8,
           col=cols, font=2, pos = 4)
    }
  })
}

#' Show average ratings by clusters as a barplot.
#'
#' @param metric Euclidean or hyperbolic.
#' @export
showBarPlot <- function(metric=NULL) {
  if(metric == 0)
    return(NULL)
  barord <- menu(c("Alphabetical", "Height",
                   bold(magenta("Plot menu"))))
  if (barord == 3) {
    return(NULL)
  }
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    MN <- SD <- SS <- rep(0, nclust)
    if(!is.null(metric)) {
      ttl <- paste(colnames(cmapdat$ratings)[metric+2],
                   colnames(cmapdat$cohorts)[cmapdat$cohortCol])
      ratingstmp <- ratings[which(cohorts[,cohortCol] == 1),]
      for (i in 1:length(unique(groups))) {
        cardsInCluster <- names(which(groups == i))
        cardsInCluster <- which(ratingstmp$StatementID %in% cardsInCluster)
        rtgsmm <- ratingSummary(ratingstmp[cardsInCluster,])
        rownames(rtgsmm) <- c(rownames(DS)[which(groups == i)], "total")
        MN[i] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_Mean"))]
        SD[i] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_SD"))]
        SS[i] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_N"))]
      }
      cname <- clusterNames()
      if (barord == 1) {
        ord <- order(clusterNames())
      }
      if (barord == 2) {
        ord <- order(MN)
      }
      MN <- MN[ord]
      SD <- SD[ord]
      SS <- SS[ord]
      cname <- cname[ord]
      barcols <- cols[ord]
    }
    SEM <- SD/sqrt(SS)
    Yl = ceiling(max(MN+SEM+0.5, na.rm = TRUE))
    bp <- barplot(MN, main=ttl,col=barcols,ylim=c(0,Yl),
                  names.arg=cname,cex.names=0.8,las=3,
                  width=rep(0.8,nclust),space=0.2)
    segments(bp, MN - SEM, bp, MN + SEM, lwd=2, col="grey66")
    segments(bp - 0.1, MN - SEM, bp + 0.1, MN - SEM, lwd=2, col="grey66")
    segments(bp - 0.1, MN + SEM, bp + 0.1, MN + SEM, lwd=2, col="grey66")
  })
}

#' A parallel coordinate plot.
#'
#' @export
showParallelCoordinates <- function() {
  with(cmapdat,{
    choices <- expand.grid(colnames(ratings)[-(1:2)], colnames(cohorts))
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    MN <- matrix(0, nrow=nclust, ncol=length(pcoordChoice))
    for (jj in 1:length(pcoordChoice)) {
      vars <- unlist(strsplit(pcoordChoice[jj],"/"))
      metric <- which(colnames(ratings) == vars[1])
      cohortCol <- which(colnames(cohorts) == vars[2])
      ttl <- paste(colnames(ratings)[metric],
                   colnames(cmapdat$cohorts)[cmapdat$cohortCol])
      ratingstmp <- ratings[which(cohorts[,cohortCol] == 1),]
      for (i in 1:length(unique(groups))) {
        cardsInCluster <- names(which(groups == i))
        cardsInCluster <- which(ratingstmp$StatementID %in% cardsInCluster)
        rtgsmm <- ratingSummary(ratingstmp[cardsInCluster,])
        rownames(rtgsmm) <- c(rownames(DS)[which(groups == i)], "total")
        MN[i, jj] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric],"_Mean"))]
        #SD[i] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_SD"))]
        #SS[i] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric+2],"_N"))]
      }
      if (jj == 1) {
        plot(rep(1, nclust), MN[,1], ylim=c(0.5,5), xlim=c(0.2,length(pcoordChoice)+0.8),
             main="", col=cols, axes=F, xlab="", ylab="Rating", pch=19, cex=0.3)
        axis(2)
        grid()
      } else {
        points(rep(jj, nclust), MN[,jj], col=cols, pch=19, cex=0.3)
        for (k in 1:nclust)
          lines(c(jj-1,jj), c(MN[k,jj-1], MN[k,jj]), col=cols[k])
      }
      text(jj, 0.8, cmapdat$pcoordChoice[jj], cex=0.6)
    }
    cnms <- clusterNames()
    text(rep(length(cmapdat$pcoordChoice)+0.1, length(cnms)),
         seq(1,5,length=length(cnms)), cnms, col=cols, cex=0.8, pos=4,font=2)

  })
}


#' A Go-Zone plot.
#'
#' @param refline Show reference line at the mean (default) or median.
#' @export
showGoZone <- function() {
  #  with(cmapdat,{
  choices <- expand.grid(colnames(cmapdat$ratings)[-(1:2)],
                         colnames(cmapdat$cohorts))
  distM <- cmapdat$D2e
  if(cmapdat$distmetric == "Hyperbolic")
    distM <- cmapdat$D2h
  ttl <- ""
  cols <- clusterCols(cmapdat$nclust, cmapdat$clrscm)
  fit.Clust <- hclust(as.dist(distM), method=clustMethod)
  groups <- cutree(fit.Clust, k=cmapdat$nclust)
  MN <- matrix(0, nrow=1+nrow(distM), ncol=length(cmapdat$pcoordChoice))
  for (jj in 1:2) {
    vars <- unlist(strsplit(cmapdat$pcoordChoice[jj],"/"))
    metric <- which(colnames(cmapdat$ratings) == vars[1])
    cohortCol <- which(colnames(cmapdat$cohorts) == vars[2])
    ttl <- paste(colnames(cmapdat$ratings)[metric],
                 colnames(cmapdat$cohorts)[cohortCol])
    ratingstmp <- cmapdat$ratings[which(cmapdat$cohorts[,cohortCol] == 1),]
    #cardsInCluster <- names(which(groups == i))
    cardsInCluster <- ratingstmp$StatementID#) %in% cardsInCluster)
    rtgsmm <- ratingSummary(ratingstmp)
    rownames(rtgsmm) <- c(rownames(cmapdat$DS), "total")
    MN[, jj] <- rtgsmm[,which(colnames(rtgsmm) ==
                                paste0(colnames(cmapdat$ratings)[metric],"_Mean"))]
  }
  kendall <- cor.test(MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2], method = "kendall")
  Kendall <- paste0("Kendall's tau=", format(kendall$estimate, digits=2),
                    " (p=", format(kendall$p.value, digits=3),")")
  plot( MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2], ylim=c(1,5), xlim=c(1,5.9),
        main=Kendall, col=0, axes=F, xlab=cmapdat$pcoordChoice[1],
        ylab=cmapdat$pcoordChoice[2], pch=19, cex=0.6, cex.main=1)
  text(MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2],
       rownames(cmapdat$DS), cex=0.7,
       col=cols[groups], font=2)
  axis(1, at = 1:5); axis(2)
  grid()
  rect(xleft = 1, ybottom = 1, xright = MN[nrow(rtgsmm),1],
       ytop = MN[nrow(rtgsmm),2], col=rgb(0.5, 0, 0.5, .1),
       border = "white")
  rect(xleft = 1, ybottom = MN[nrow(rtgsmm),2], xright = MN[nrow(rtgsmm),1],
       ytop = 5, col=rgb(0.5, 0.5, 0, .1),
       border = "white")
  rect(xleft = MN[nrow(rtgsmm),1], ybottom = 1, xright = 5,
       ytop = MN[nrow(rtgsmm),2], col=rgb(0.5, 0, 0, .1),
       border = "white")
  rect(xleft = MN[nrow(rtgsmm),1], ybottom = MN[nrow(rtgsmm),2], xright = 5,
       ytop = 5, col=rgb(0, 0.5, 0.5, .1),
       border = "white")
  abline(v=MN[nrow(rtgsmm),1], col="grey", lwd=2)
  abline(h=MN[nrow(rtgsmm),2], col="grey", lwd=2)
  cnms <- clusterNames()
  rect(5.1, 0, 7, 7, col = "white", border= "white")
  text(rep(5.1, length(cnms)), seq(1,5,length=length(cnms)),
       cnms, col=cols, cex=0.8, pos=4,font=2)

  #  })
}


#' A report about the sorters.
#'
#' @export
sorterReport <- function() {
  with(cmapdat, {
    for (i in 1:length(adj.mat)) {
      cardsInPiles <- which(rowSums(adj.mat[[i]])>0)
      cardsSorted <- length(cardsInPiles)
      adjM <- adj.mat[[i]][cardsInPiles, cardsInPiles]
      diag(adjM) <- 1
      noPiles <- length(unique(apply(adjM, 1, paste0, collapse="")))
      cat(paste0("Sorter ", i, " sorted ", cardsSorted, " cards into ", noPiles," piles\n"))
    }
  })
  readline("Press any key to continue. ")
}

#' A report about the raters.
#'
#' @export
raterReport <- function() {
  with(cmapdat, {
    print(summary(cmapdat$demographics[,-1]))
  })
  readline("Press any key to continue. ")
}

#' A report about the statements.
#'
#' @export
statementReport <- function() {
  with(cmapdat, {
    distM <- D2e
    eta_ij <- eta_ijk_e[[nclust+1]]
    if(distmetric == "Hyperbolic") {
      distM <- D2h
      eta_ij <- eta_ijk_h[[nclust+1]]
    }
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    jidx <- colMeans(eta_ij)
    groups <- cutree(fit.Clust, k=nclust)
    retdf <- as.data.frame(matrix(0, ncol=4+5*(ncol(ratings)-2), nrow=nclust+nrow(DS)))
    k <- 0
    for (i in 1:nclust) {
      cardsInCluster <- names(which(groups == i))
      jidx_i <- c(jidx[which(groups == i)], mean(jidx[which(groups == i)]))
      cNames <- c(cardNames[which(groups == i)],"")
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      rtgsmm <- ratingSummary(ratings[cardsInCluster,])
      rtgsmm <- data.frame(cNames, c(names(cmapdat$x)[which(groups == i)],""),
                           rep(i, nrow(rtgsmm)), jidx_i, rtgsmm)
      colnames(rtgsmm)[1:4] <- c("Statement", "CardNo", "ClusterNo", "JI")
      rownames(rtgsmm) <- c(rownames(DS)[which(groups == i)], paste0("total_",i))
      retdf[(1+k):(k+nrow(rtgsmm)),] <- rtgsmm
      k <- k + nrow(rtgsmm)
    }
    colnames(retdf) <- colnames(rtgsmm)
    fn <- paste0(dataDir,"output/StatementSummary", nclust,".csv")
    write.csv(retdf, file = fn, row.names = FALSE)
    cat(paste0("Report saved to ", fn,".\n"))
    readline("Press any key to continue. ")
  })
}

#' Run analysis of variance by clusters.
#'
#' @export
clusterANOVA <- function() {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    fn <- paste0(dataDir,"output/ANOVA", nclust,".txt")
    Cluster <- rep(1, nrow(ratings))
    for (i in 1:length(groups)) {
      cardsInCluster <- names(which(groups == i))
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      Cluster[cardsInCluster] <- i
    }
    Cluster <- as.factor(Cluster)
    for (i in 3:ncol(ratings)) {
      ftest <- summary(aov(ratings[,i] ~ Cluster))
      cat("Analysis of Variance: Response=", colnames(ratings)[i],"\n")
      print(ftest)
      cat(paste0(rep("_",80),collapse=""),"\n")
    }
    sink(file=fn)
    for (i in 3:ncol(ratings)) {
      ftest <- summary(aov(ratings[,i] ~ Cluster))
      cat("Analysis of Variance: Response=", colnames(ratings)[i],"\n")
      print(ftest)
      cat(paste0(rep("_",80),collapse=""),"\n")
    }
    sink()
    cat(paste0("Output saved to ", fn,".\n"))
    readline("Press any key to continue. ")
  })
}

#' All pairwise comparisons with Tukey's method for multiple comparisons.
#'
#' @export
clusterTUKEY <- function() {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    Cluster <- rep(1, nrow(ratings))
    fn <- paste0(dataDir,"output/Tukey", nclust,".txt")
    for (i in 1:length(groups)) {
      cardsInCluster <- names(which(groups == i))
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      Cluster[cardsInCluster] <- i
    }
    Cluster <- as.factor(Cluster)

    for (i in 3:ncol(ratings)) {
      cat("Analysis of Variance: Response=", colnames(ratings)[i],"\n")
      print(TukeyHSD(aov(ratings[,i] ~ Cluster)))
      cat(paste0(rep("_",80),collapse=""),"\n")
    }
    sink(file=fn)
    for (i in 3:ncol(ratings)) {
      cat("Analysis of Variance: Response=", colnames(ratings)[i],"\n")
      print(TukeyHSD(aov(ratings[,i] ~ Cluster)))
      cat(paste0(rep("_",80),collapse=""),"\n")
    }
    sink()
    cat(paste0("Output saved to ", fn,".\n"))
    readline("Press any key to continue. ")
  })
}

#' Show the header of the menu
topLine <- function() {
  cat("\014" %+% blue(underline(bold("\nRCMap command-line interface.\n"))))
}

#' Convert menu choices to text
#' @param varType The type of variable (e.g., "clusteringMethod", "distance")
#' @param choice The number of the selected menu item
toText <- function(varType, choice) {
  if (varType == "clusteringMethod") {
    #    strs <- c("ward.D", "ward.D2", "single","complete", "average","mcquitty",
    #              "median", "centroid")
    strs <- "ward.D2"
  }
  if (varType == "distance") {
    strs <- c("Euclidean", "Hyperbolic")
  }
  if (varType == "MPparam") {
    strs <- c("median", "mean", "max")
  }
  if (varType == "colorScheme") {
    strs <- c("rcmap", "rainbow")
  }
  return(strs[choice])
}

#' Colors for the dendrograms.
#' @param dn The dendrogram object.
colLab <- function(dn) {
  clsmem <- cmapdat$clusMember
  labelColors <- clusterCols(length(unique(clsmem)), cmapdat$clrscm)
  if (is.leaf(dn)) {
    a <- attributes(dn)
    labCol <- labelColors[cmapdat$clusMember[which(names(cmapdat$clusMember) == a$label)]]
    attr(dn, "nodePar") <- c(a$nodePar, list(lab.col = labCol, pch=19,
                                             cex=0.5, col=labCol, lab.cex=0.7))
  }
  dn
}

#' Concept mapping text-based menu-driven interface.
#'
#' This is the main function in the package. Run RCMapMenu() to perform the concept mapping analysis via a hierarchical menu.
#' @importFrom crayon bold red blue green yellow cyan magenta italic silver %+% underline
#' @importFrom utils menu read.csv select.list stack
#' @importFrom graphics abline axis grid points rect text
#' @importFrom stats as.dendrogram as.dist cor.test cutree dendrapply dist hclust is.leaf median model.matrix na.exclude sd
#' @importFrom factoextra fviz_nbclust
#'
#' @param enc Encoding (default=UTF-8).
#' @param sep Column separator for csv files (default=,).
#' @export
#' @examples
#' \donttest{
#' #RCMapMenu()
#' }
RCMapMenu <- function(enc = "UTF-8", sep=",") {
  menuLevel <- 0
  while(menuLevel >= 0) {
    topLine()
    ## Main menu
    if(menuLevel == 0) {
      show_all <- ifelse(!exists("cmapdat"), FALSE,
                         ifelse(is.null(cmapdat$x), FALSE, TRUE))
      rcmapmenu <- menu(c(ifelse(show_all,"Choose the data folder",bold("Choose the data folder")),
                          ifelse(show_all,"Summary",silver("summary")),
                          ifelse(show_all,"Settings",silver("Settings")),
                          ifelse(show_all,"Plots",silver("Plots")),
                          ifelse(show_all,"Reports",silver("Reports")),
                          ifelse(show_all,"Analysis",silver("Analysis")),
                          bold(magenta("R prompt"))), title=bold("  Main menu"))
      if (rcmapmenu %in% c(0,7))
        menuLevel <- -1 ## get the R prompt back
      if (rcmapmenu == 1) {
        loadRCMapData(enc)
        menuLevel <- 0
      }
      if(!exists("cmapdat"))
        next
    }
    ## general properties of the data:
    if(rcmapmenu == 2) {
      topLine()
      # Rating variables
      cat(sprintf("Data folder: %s\nNumber of sorters: %d\nNumber of statements: %d\nIssues:\n%s\n",
                  cmapdat$dataDir, cmapdat$n.indiv, length(cmapdat$cardNames),cmapdat$issues))
      if(length(cmapdat$splhalf) == 2) {
        cat(blue("Mean correlation between split halves:",
                 bold(sprintf("%1.2f",mean(cmapdat$splhalf$cors))),
                 "(using 20 random splits, distance=", cmapdat$splhalf$distmetric,
                 #, func=",cmapdat$splhalf$func,
                 ")\n"))
      }
      summMenu <- menu(c("Perform split-half analysis",
                         "Perform leave-one-out analysis", # (Jaccard Index)",
                         bold(magenta("Main menu"))))
      if(summMenu == 1) {
        cmapdat$splhalf <<- splitHalf(cmapdat$adj.mat, B=20,
                                      disttype = cmapdat$distmetric,
                                      plotit=TRUE, seed=23456)
      }
      if(summMenu == 2) {
        if (cmapdat$distmetric == "Euclidean")
          plotjaccard(cmapdat$eta_ijk_e)
        else
          plotjaccard(cmapdat$eta_ijk_h)
      }
      menuLevel <- 1
      if(summMenu %in% c(3))
        menuLevel <- 0
      if(summMenu %in% c(0))
        menuLevel <- -1
    }
    ## settings
    if(rcmapmenu == 3) {
      topLine()
      showSettingMenu <- TRUE
      if(exists("settingmenu"))
        if(settingmenu == 6)
          showSettingMenu <- FALSE
      if(showSettingMenu)  {
        settingmenu <- menu(c("Choose the distance metric",
                              #"Choose the clustering method",
                              #"Choose the misplacement index options",
                              "Choose the number of clusters",
                              "Set cluster names",
                              "Choose color scheme",
                              bold(magenta("Main menu"))), title=bold("  Settings"))
        menuLevel <- 1
        if(settingmenu == 5) { menuLevel <- 0 }
        if(settingmenu == 0) { menuLevel <- -1 }
      }
      if(settingmenu == 1) {
        topLine()
        dtype <- rcmenu(toText("distance"),
                        title=bold("  Distance function"),
                        slctd=which(toText("distance") == cmapdat$distmetric))
        if(dtype > 0)
          cmapdat$distmetric <<- toText("distance", dtype)
      }
      # if(settingmenu == 2) {
      #   topLine()
      #   cmeth <- rcmenu(toText("clusteringMethod"),
      #                   title=bold("  Clustering method"),
      #                   slctd=which(toText("clusteringMethod") == cmapdat$clustMethod))
      #   if(cmeth > 0)
      #     cmapdat$clustMethod <<- toText("clusteringMethod", cmeth)
      # }
      # if(settingmenu == 3) { ## misplacement index
      #   topLine()
      #   mpfunc <- rcmenu(toText("MPparam"),
      #                    title=bold("  Misplacement index parameters"),
      #                    slctd=which(toText("MPparam") == cmapdat$MPfunc))
      #   if(mpfunc > 0)
      #     cmapdat$MPfunc <<- toText("MPparam", mpfunc)
      # }
      if(settingmenu == 2) { ## the number of clusters
        topLine()
        distM <- cmapdat$D2e
        if(cmapdat$distmetric == "Hyperbolic")
          distM <- cmapdat$D2h
        clustNo <- menu(c("Within cluster sums of squares",
                          "Average silhouette",
                          "Set manually",
                          bold(magenta("Main menu"))),
                        title=bold("  Method used to select the number of clusters"))
        menuLevel <- ifelse(clustNo == 4, 0, 1)
        if(clustNo == 1) {
          nbcwss <- fviz_nbclust(cbind(cmapdat$x, cmapdat$y),
                                 factoextra::hcut,
                                 method = c("wss"), diss=distM,
                                 k.max = length(cmapdat$cardNames)/3)
          difflogs <- -diff(log(nbcwss$data$y))
          plot(nbcwss$data$y, pch=19, col=4, cex=0.6,
               ylim=c(0,1.05*max(nbcwss$data$y)),
               xlab="No. of clusters", ylab="Within cluster sums of squares",
               main="", cex.main=1, axes=F)
          grid(); axis(1); axis(2)
          abline(h=0, lwd=2)
          if(all(difflogs >= 0.1)) {
            cmapdat$nclust <<- length(nbcwss$data$y)
          } else {
            cmapdat$nclust <<- min(which(difflogs < 0.1))
          }
          update_clusters()
          points(cmapdat$nclust, nbcwss$data$y[cmapdat$nclust], cex=1.5, pch=18,
                 col="navyblue")
        }
        if(clustNo == 2) {
          nbcsil <- fviz_nbclust(cbind(cmapdat$x, cmapdat$y),
                                 factoextra::hcut,
                                 method = c("sil"), diss=distM,
                                 k.max = length(cmapdat$cardNames)/3)
          plot(nbcsil)
          cmapdat$nclust <<- which.max(nbcsil$data$y)
          update_clusters()
        }
        if(clustNo == 3) {
          fit.Clust <- hclust(as.dist(distM), method=clustMethod)
          #cmapdat$clusMember <- cutree(fit.Clust, k=cmapdat$nclust)
          showDendrogram()
          curNC <- cmapdat$nclust
          tmpnclust <- rcmenu(1:floor(length(cmapdat$cardNames)/3),
                              title = bold("Set the number of clusters"), slctd = curNC)
          if(tmpnclust > 0)
            cmapdat$nclust <<- tmpnclust
          update_clusters()
          #cmapdat$clusMember <- cutree(fit.Clust, k=cmapdat$nclust)
          showDendrogram()
        }
      }
      if(settingmenu == 3) { ## the cluster names
        topLine()
        clusterNum <- menu(clusterNames(),
                           title = bold("  Select cluster number (or 0 to return to the main menu)"))
        if (clusterNum == 0) {
          settingmenu = 1
          next
        }
        cat(clusterNum, "[",clusterNames()[clusterNum],"]\n")
        groups <- clusterConfig()
        cat(bold("Statements in the cluster"),"\n")
        cardsInCluster <- which(groups == clusterNum)
        cat(blue(cmapdat$cardNames[cardsInCluster]),sep="\n")
        tmplbl <- cmapdat$pileLabels[which(cmapdat$pileLabels$cardNum %in% cardsInCluster),]
        if(nrow(tmplbl) > 0) {
          tmplbl <- sort(table(tmplbl[,1]), decreasing = TRUE)
          toShow <- min(5, nrow(tmplbl))
          cat(bold("Suggested names"),"\n")
          cat(green(names(tmplbl[1:toShow])),sep="\n")
        }
        clNames <- clusterNames()
        cat("\nEnter a name [", clNames[clusterNum],"]\n")
        clName <- scan("",quiet = TRUE, what = "chr", n=1, sep = "\n")
        if (length(clName) == 0)
          clName <- clNames[clusterNum]
        key <- sprintf("clust_%02d_%02d", length(clNames), clusterNum)
        cmapdat$clustname[[key]] <<- clName
      }
      if(settingmenu == 4) {
        topLine()
        clrs <- rcmenu(toText("colorScheme"),
                       title=bold("  Color scheme"),
                       slctd=which(toText("colorScheme") == cmapdat$clrscm))
        if(clrs > 0)
          cmapdat$clrscm <<- toText("colorScheme", clrs)
      }
    }
    ## plots
    if(rcmapmenu == 4) {
      topLine()
      menuLevel <- 1
      plotmenu <- menu(c(blue("Point map (MDS)"),
                         blue("Clusters (rays)"),
                         blue("Clusters (polygons)"),
                         blue("Dendrogram"),
                         blue("Phylogenic tree"),
                         blue("Misplacement"),
                         crayon::green("Statement Rating (Map)"),
                         crayon::green("Statement Rating (Dot chart)"),
                         crayon::green("Cluster Rating (Map)"),
                         crayon::green("Cluster Rating (Bar chart)"),
                         cyan("Parallel Coordinates"),
                         cyan("GoZone"),
                         bold(magenta("Main menu"))), title=bold("  Plots"))
      if(plotmenu  == 13) { menuLevel <- 0 }
      if(plotmenu  == 0) { menuLevel <- -1 }
      if(plotmenu %in% 7:10) {
        cmapdat$cohortCol <<- menu(colnames(cmapdat$cohorts))
        if (cmapdat$cohortCol == 0)
          next
      }
      if(plotmenu %in% c("11", "12")) {
        ormore <- ""
        if(plotmenu == "11")
          ormore <- " or more "
        pcoordChoice <-
          select.list(apply(expand.grid(colnames(cmapdat$ratings)[-(1:2)],
                                        colnames(cmapdat$cohorts)), 1, paste,
                            collapse="/"), multiple=TRUE,
                      title=bold("Choose 2" %+% ormore %+% " cohort/variable combinations (space separated)"))
        if(length(unique(setdiff(pcoordChoice, 0))) < 2)
          next
        cmapdat$pcoordChoice <<- pcoordChoice
      }
      switch(plotmenu,
             "1" = showMDSPlot(),
             "2" = showClusterPlot(),
             "3" = showClusterPlot(type="polygon"),
             "4" = showDendrogram(),
             "5" = showDendrogram(phylo = TRUE),
             "6" = showMDSPlot(metric = "MPindex"),
             "7" = showMDSPlot(metric = menu(c(colnames(cmapdat$ratings)[-(1:2)]),
                                             title = "Rating variable")),
             "8" = showDotPlot(metric = menu(c(colnames(cmapdat$ratings)[-(1:2)]),
                                             title = "Rating variable")),
             "9" = showClusterPlot(metric = menu(c(colnames(cmapdat$ratings)[-(1:2)]),
                                                 title = "Rating variable")),
             "10" = showBarPlot(metric = menu(c(colnames(cmapdat$ratings)[-(1:2)]),
                                              title = "Rating variable")),
             "11" = showParallelCoordinates(),
             "12" = showGoZone()
      )
    }
    if(rcmapmenu == 5) { # reports
      topLine()
      menuLevel <- 1
      reportsmenu <- menu(c("Sorters", "Raters" ,"Statements",
                            bold(magenta("Main menu"))), title=bold("  Report"))
      if(reportsmenu == 0) { menuLevel <- -1 }
      if(reportsmenu == 4) { menuLevel <- 0 }
      switch (reportsmenu,
              "1" = sorterReport(),
              "2" = raterReport(),
              "3" = statementReport()
      )
    }
    if(rcmapmenu == 6) { # statistical inference methods
      topLine()
      menuLevel <- 1
      analysismenu <- menu(c("Between-cluster ANOVA", "Tukey - all cluster pairs",
                             bold(magenta("Main menu"))), title=bold("  Analysis"))
      menuLevel <- 1
      if(analysismenu == 3) { menuLevel <- 0 }
      if(analysismenu == 0) { menuLevel <- -1 }
      if(analysismenu %in% c(1, 2))
        switch (analysismenu,
                "1" = clusterANOVA(),
                "2" = clusterTUKEY(),
        )
    }
  }
  save(cmapdat, file=paste0(cmapdat$dataDir,"/CMapSession.RData"))
}


#' Load the input files for a concept mapping project.
loadRCMapData <- function(enc = "UTF-8", sep=",") {
  topLine()
  dataDir <- paste0(choose_directory(),"/")
  cat(bold(crayon::green("Opening project folder:" %+% dataDir %+% "...\n\n")))
  cmapdat <<- initCMap(dataDir, enc, sep)
  update_clusters()
}

#' Choose the project's directory.
#'
#' From the easycsv package.
#' @return The user's selected directory.
#' @export
choose_directory = function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

choose_dir <- function(){
  os <- getOS()
  if(os == "windows"){
    directory <- utils::choose.dir()
    directory <- paste0(directory,"\\")
  }
  if(os == "Linux"){
    directory <- system("zenity --file-selection --directory", intern = TRUE)
  }
  if(os == "MacOSX"){
    system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
           intern = FALSE, ignore.stderr = TRUE)
    directory <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  }
  return(directory)
}

#' Return the operating system info
getOS <- function() {
  pl <- tolower(.Platform$OS.type)
  if(pl == "windows")
    return(pl)
  si <- as.list(Sys.info())
  if(tolower(si$sysname) == "linux")
    return("Linux")
  if(tolower(si$sysname) == "darwin")
    return("MacOSX")
}


#' Return the clusters for the concept mapping project
clusterConfig <- function() {
  distM <- cmapdat$D2e
  if(cmapdat$distmetric == "Hyperbolic")
    distM <- cmapdat$D2h
  fit.Clust <- hclust(as.dist(distM), method=clustMethod)
  cutree(fit.Clust, k=cmapdat$nclust)
}

#' Return the cluster names.
clusterNames <- function() {
  lbl <- rep("",cmapdat$nclust)
  for(cl in 1:cmapdat$nclust) {
    key <- sprintf("clust_%02d_%02d",cmapdat$nclust,cl)
    if (length(cmapdat$clustname[[key]]) > 0) {
      lbl[cl] <- cmapdat$clustname[[key]]
    } else {
      lbl[cl] <- cl
    }
  }
  lbl
}

#' Set cluster colors.
#'
#' Either use the RCMap scheme which has 21 colors, or the rainbow function.
#' @param n Number of colors.
#' @param type Either rcmap (default) or rainbow.
#' @importFrom grDevices rainbow
#' @export
clusterCols <- function(n, type="rcmap") {
  if(type == "rainbow")
    return(rainbow(n))
  rcmappal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#A65628", "#F781BF", "#999999", "#8000FF", "#FF00BF",
                "#00FF40", "#0040FF", "#FFBF00",
                "#E41A1CAA", "#377EB8AA", "#4DAF4AAA", "#984EA3AA",
                "#FF7F00AA", "#A65628AA", "#F781BFAA", "#999999AA")
  if(n <= length(rcmappal))
    return(rcmappal[1:n])
  # ideally should not get here - too many clusters can't be distinguished
  # by colors:
  return(rainbow(n))
}


rcmenu <- function (choices, graphics = FALSE, title = NULL, slctd=NULL) {
  if (!interactive())
    stop("rcmenu() cannot be used non-interactively")
  nc <- length(choices)
  if (length(title) && nzchar(title[1L]))
    cat(title[1L], "\n")
  op <- paste0(format(seq_len(nc)), ": ", choices)
  if (nc > 10L) {
    fop <- format(op)
    nw <- nchar(fop[1L], "w") + 2L
    ncol <- getOption("width")%/%nw
    if (ncol > 1L)
      op <- paste0(fop, c(rep.int("  ", min(nc, ncol) -
                                    1L), "\n"), collapse = "")
  }
  cat("", op, "", sep = "\n")
  repeat {
    prmpt <- "Selection: "
    if (!is.null(slctd))
      prmpt <- "Selection [" %+% as.character(slctd) %+% "] :"
    ind <- readline(prmpt)
    if(nchar(ind) == 0)
      return(0)
    if(ind %in% as.character(0:nc))
      return(as.numeric(ind))
    cat(gettext("Choose an item number from the menu, or hit Enter to exit\n"))
  }
}


troubleshooting <- c(
  "cardfilename" = "Make sure the file is spelled correctly (Statements.csv, case-sesitive), and that the data folder contains all the input files.",
  "cardfilesep" = "Make sure you use the correct column separator. If your statements contain commas, use tab as a column separator, and invoke the program with RCMapMenu(sep='\t')",
  "pilesfilename" = "Make sure the file is spelled correctly (SortedCards.csv, case-sesitive), and that the data folder contains all the input files.",
  "pilesfilesep" = "Make sure you use the correct column separator. The first column in SortedCards.csv must contain the sorter IDs, the second column has to contain pile labels (text), and columns 3,4,... contain the cards in each pile.",
  "ratingsfilename" = "Make sure the file is spelled correctly (Ratings.csv, case-sesitive), and that the data folder contains all the input files.",
  "ratingsfilesep" = "The first two columns must be named RaterID	and StatementID (case-sensitive) and StatementID must match the values in the Statements.csv file",
  "demographicsfilename" = "Make sure the file is spelled correctly (Demographics.csv, case-sesitive), and that the data folder contains all the input files.",
  "demographicsfilesep" = "The first column must be named RaterID (case-sensitive) and must match the values in the Ratings.csv file",
  "demographicsnodata" = "Either no columns in the file, or all columns contain categorical data with distinct values, or all values in each column are identical.",
  "non-numeric" = "has a non-numeric value where a card number is expected"
)
