
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
#' @param piledat The pile sorting data frame as read from \code{SortedCards.csv}:
#'   first column sorter ID, second column pile label, columns 3+ card numbers.
#' @param cardNames Data frame of statements from \code{read_statements()}, with
#'   at least a \code{StatementID} column used to map card numbers to matrix rows/cols.
#' @param sorters_dict Data frame with columns \code{orig_ID} (original sorter
#'   identifiers as they appear in \code{piledat}) and \code{seq_ID} (sequential
#'   integer index 1..n used internally).
#' @param showWarnings If \code{TRUE} (default), print warnings about potential
#'   data problems such as unsorted cards, all cards in one pile, or cards
#'   appearing in multiple piles.
#' @return A list with two elements: \code{mat.list} (a list of n SxS 0/1
#'   adjacency matrices, one per sorter) and \code{issues} (a character string
#'   accumulating any warnings found, empty string if none).
#' @export
getAdjMatrices <- function(piledat, cardNames, sorters_dict, showWarnings=TRUE) {
  issues <- ""
  if (ncol(piledat) < 3)
    issues <- issues %+% red("Pile sorting file must have at least 3 columns!\n")
  sorters <- unique(piledat[ ,1])
  # Size the adjacency matrix by the total number of statements in cardNames,
  # not just the number of cards that appear in piles. Using length(cardsused)
  # would make nCards < nrow(cardNames) when some cards are never sorted, and
  # since which(cardNames[,1] %in% cards) returns positions within cardNames
  # (not ID values), those positions can exceed length(vct) and go out of bounds.
  nCards <- nrow(cardNames)
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


#' Add adjacency matrices from pile-sorting data and return a distance matrix.
#'
#' Sums n SxS adjacency matrices to create a similarity matrix M, then computes
#' the dissimilarity as \code{n*(J - I) - M}, where J is the all-ones matrix and
#' I is the identity. A tiny jitter is applied to break ties when two statements
#' have identical pile distributions across all sorters.
#' @param adjMats A list of n SxS 0/1 adjacency matrices (one per sorter), as
#'   returned by \code{getAdjMatrices()}.
#' @return An SxS dissimilarity matrix. The diagonal is zero; off-diagonal entry
#'   (i,j) is the number of sorters who did NOT place cards i and j in the same
#'   pile (plus a negligible jitter for numerical stability).
#' @param seed Integer random seed set before applying jitter, for
#'   reproducibility (default 154204). Affects the global RNG state.
#' @details Sets the random seed before applying jitter. This ensures
#'   reproducibility but will affect the global RNG state.
#' @export
distanceMatrix <- function(adjMats, seed=154204) {
  M <- Reduce("+", adjMats)
  n <- length(adjMats)
  S <- nrow(M)
  Imat <- diag(1, S, S)
  Jmat <- matrix(1, S, S)
  set.seed(seed)
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


#' Convert 2-D coordinates in the unit disk to Poincare (hyperbolic) distances.
#'
#' Takes a matrix of 2-D coordinates of points inside the unit disk and returns
#' the pairwise Poincare distances in hyperbolic space.
#' @param x A matrix of 2-D coordinates (one point per row) for points lying
#'   inside the unit disk. This is NOT a distance matrix — it is the raw
#'   coordinate output from an MDS projection.
#' @return A symmetric matrix of pairwise hyperbolic (Poincare) distances.
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

#' Build co-membership matrices for a range of cluster solutions.
#'
#' Runs hierarchical clustering on a distance matrix and, for each number of
#' clusters k from 2 to \code{max_nc}, builds an MxM binary co-membership
#' matrix where entry (i,j) is 1 if statements i and j are in the same cluster.
#' @param dists An MxM distance matrix (e.g. Euclidean or hyperbolic 2-D distances).
#' @param M Number of statements (rows/columns of \code{dists}).
#' @param disttype Character label for the distance type (e.g. \code{"Euclidean"}
#'   or \code{"Hyperbolic"}). Currently used for labelling only; does not alter
#'   the computation.
#' @param max_nc Maximum number of clusters to consider.
#' @param clustMethod Hierarchical clustering method passed to \code{hclust}
#'   (default \code{"ward.D2"}).
#' @return A list of length \code{max_nc - 1}. Element \code{[[k-1]]} is an
#'   MxM 0/1 co-membership matrix for the k-cluster solution (diagonal is 0).
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

#' Run MDS and build co-membership matrices for Euclidean and hyperbolic distances.
#'
#' Optionally drops one sorter (leave-one-out), computes the dissimilarity
#' matrix, runs 2-D MDS, then calls \code{get_cmat} for both Euclidean and
#' hyperbolic distance metrics across all cluster solutions up to \code{max_nc}.
#' @param adjmat A list of n SxS adjacency matrices (one per sorter).
#' @param j Index of the sorter to leave out (0 = use all sorters, default).
#'   If \code{j} exceeds the number of sorters it is silently reset to 0.
#' @param max_nc Maximum number of clusters to consider (upper bound for k).
#' @param clustMethod Hierarchical clustering method passed to \code{hclust}
#'   (default \code{"ward.D2"}).
#' @return A named list with two elements: \code{CmatE} and \code{CmatH}, each
#'   a list of \code{max_nc - 1} co-membership matrices (see \code{get_cmat})
#'   based on Euclidean and hyperbolic 2-D distances respectively.
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


#' Assess reliability using split-half analysis.
#'
#' Randomly splits the sorters into two halves B times and computes the
#' Pearson correlation between the pairwise distance matrices derived from each
#' half. Higher mean correlation indicates greater reliability of the map.
#' @param adjMat A list of n SxS adjacency matrices (one per sorter).
#' @param B The number of random splits (default=10).
#' @param disttype The distance metric to be used: "Hyperbolic" for Poincare
#'   distances, anything else for Euclidean (default="Hyperbolic").
#' @param seed The base random seed; each split uses seed+i.
#' @param plotit If TRUE, show a plot of the per-split correlations (default=FALSE).
#' @return A list with two elements: \code{cors} (numeric vector of length B
#'   with the correlation for each split) and \code{distmetric} (the distance
#'   metric used as a character string).
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


#' Plot Jaccard index stability across cluster solutions.
#'
#' Takes a list of per-sorter Jaccard index matrices (as returned by the
#' leave-one-out analysis in \code{initCMap}) and produces an interaction plot
#' showing the Jaccard index for each statement across different numbers of
#' clusters. The median, 25th and 75th percentiles are highlighted.
#' @param eta_ij A list of matrices (one per number-of-clusters solution),
#'   where each matrix has rows = sorters and columns = statements, and cell
#'   values are Jaccard indices measuring clustering stability.
#' @param threshold Jaccard values below this are counted as "unstable"
#'   when computing the proportion of sorters with low stability per statement
#'   (default 0.3).
#' @export
plotjaccard <- function(eta_ij, threshold=0.3) {
  M <- matrix(0, nrow=length(eta_ij), ncol=ncol(eta_ij[[1]]))
  M2 <- matrix(0, nrow=0, ncol=3)
  for (i in 1:length(eta_ij)) {
    M[i, ] <- colMeans(eta_ij[[i]] < threshold)
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


## Check that the file containing the statements is loaded properly
#'
#' Read statements file with basic validation and normalization.
#'
#' The statements file should be a CSV with two columns (ID and statement text).
#' The function is forgiving about header names (it will rename the first two
#' columns to `StatementID` and `Statement`) and will normalize encodings and
#' whitespace. It validates that StatementID values are numeric and unique.
#'
#' @param filename Path to the statements CSV file.
#' @param enc File encoding (default UTF-8).
#' @param sep Field separator (default comma).
#' @return A list with elements `issues` (NULL if none) and `cardNames` (data.frame)
#' @examples
#' read_statements('Statements.csv')
#' @export
read_statements <- function(filename, enc = "UTF-8", sep=",") {
  if (!file.exists(filename))
    return(list(issues = paste("read_statements: File", filename,
                               "doesn't exist.\n", troubleshooting["cardfilename"]),
                cardNames = NULL))

  # Read with a safe try/catch to produce clearer error messages
  cardNames <- tryCatch(
    read.csv(filename, header = TRUE, encoding = enc, sep = sep, stringsAsFactors = FALSE),
    error = function(e) {
      return(structure(list(error = TRUE, message = e$message), class = "read_error"))
    }
  )
  if (inherits(cardNames, "read_error")) {
    return(list(issues = paste("read_statements: Error reading file:", cardNames$message), cardNames = NULL))
  }

  # If there are more than 2 columns, warn but continue using the first two
  if (ncol(cardNames) < 2) {
    return(list(issues = paste("read_statements:", filename,
                              "must have at least two columns (ID and statement).\n",
                              troubleshooting["cardfilesep"]), cardNames = NULL))
  }
  if (ncol(cardNames) > 2) {
    warning(sprintf("read_statements: %s has more than two columns; only the first two will be used.", filename))
    cardNames <- cardNames[, 1:2]
  }

  # Normalize column names and trim whitespace
  colnames(cardNames)[1:2] <- c("StatementID", "Statement")
  # Trim whitespace in Statement text
  cardNames$Statement <- trimws(iconv(as.character(cardNames$Statement), from = "ISO-8859-1", to = "UTF-8"))

  # Convert StatementID to character, then to numeric with checks
  raw_ids <- as.character(cardNames$StatementID)
  stmt_ids <- suppressWarnings(as.numeric(raw_ids))
  if (any(is.na(stmt_ids))) {
    bad_idx <- which(is.na(stmt_ids))
    bad_vals <- raw_ids[bad_idx]
    issues_msg <- paste0(sprintf("read_statements: 'StatementID' must be numeric. Offending values (up to 5): %s (rows: %s).",
                                 paste(head(bad_vals, 5), collapse = ","), paste(head(bad_idx, 5), collapse = ",")),
                         "\n", troubleshooting["cardfilesep"])
    return(list(issues = issues_msg, cardNames = NULL))
  }
  if (any(duplicated(stmt_ids))) {
    dupvals <- unique(stmt_ids[duplicated(stmt_ids)])
    issues_msg <- sprintf("read_statements: Duplicate StatementID values found. Examples (up to 5): %s.", paste(head(dupvals,5), collapse = ","))
    return(list(issues = issues_msg, cardNames = NULL))
  }

  cardNames$StatementID <- stmt_ids
  # Add SeqID and lowercase normalized statement for matching
  cardNames$SeqID <- seq_len(nrow(cardNames))
  cardNames$LC <- tolower(cardNames$Statement)
  # Basic content check: at least 3 statements
  if (nrow(cardNames) < 3) {
    return(list(issues = paste("read_statements:", filename, "must have some statements."), cardNames = NULL))
  }

  return(list(issues = NULL, cardNames = cardNames))
}


#' Read and validate the pile-sorting data file.
#'
#' Reads \code{SortedCards.csv} and validates that card numbers are numeric and
#' that the file has the correct structure. Also extracts pile labels for each
#' card, which are used later to suggest cluster names.
#' @param filename Path to the pile-sorting CSV file (\code{SortedCards.csv}).
#'   Expected format: column 1 = sorter ID, column 2 = pile label (text),
#'   columns 3+ = card numbers (integers matching \code{StatementID}).
#' @param ncards Total number of cards (statements), used to pre-size the
#'   column class vector when reading the file.
#' @param enc File encoding passed to \code{read.csv} (default \code{"UTF-8"}).
#' @param sep Column separator passed to \code{read.csv} (default \code{","}).
#' @return A list with elements: \code{issues} (NULL if none, otherwise an error
#'   string), \code{cardDat} (the raw data frame), \code{pileLabels} (data
#'   frame of pile label / card number pairs), and \code{pileLabelsLC}
#'   (lowercase, trimmed version of pile labels for matching).
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


#' Read and validate the ratings data file.
#'
#' Reads \code{Ratings.csv} and checks that the first two column names are
#' \code{RaterID} (or \code{UserID}/\code{SorterID}) and \code{StatementID},
#' that all StatementIDs match those in \code{Statements.csv}, and that all
#' rating values fall within \code{[1, ratingscale]}.
#' @param filename Path to the ratings CSV file (\code{Ratings.csv}).
#' @param statements Numeric vector of valid statement IDs (from
#'   \code{cardNames$StatementID}) used to cross-validate the ratings file.
#' @param ratingscale Maximum allowed rating value; ratings must be integers in
#'   \code{[1, ratingscale]} (default 5).
#' @param enc File encoding passed to \code{read.csv} (default \code{"UTF-8"}).
#' @param sep Column separator passed to \code{read.csv} (default \code{","}).
#' @return A list with elements: \code{issues} (NULL if none, otherwise an error
#'   string) and \code{ratings} (the validated data frame with \code{RaterID}
#'   as character and \code{StatementID} as factor).
read_ratings <- function(filename, statements, ratingscale=5, enc = "UTF-8", sep=",") {
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
  tmpmat <- apply(ratings[ ,3:ncol(ratings), drop=FALSE], 2,
                  function(col) suppressWarnings(as.numeric(as.character(col))))
  nonlikert <- cbind(rowSums(tmpmat < 1, na.rm=TRUE), rowSums(tmpmat > ratingscale, na.rm=TRUE))
  nonlikert <- colSums(nonlikert, na.rm = TRUE)
  if (max(nonlikert) > 0)
    return(list(issues=paste("read_ratings:", filename,
                             "The rating values must be between 1 and ", ratingscale,"\n",
                             "Columns [", colnames(ratings)[2+which(nonlikert > 0)],
                             "] have invalid values.\n"),
                ratings = NULL))
  ratings[,2] <- as.factor(ratings[,2])
  return(list(issues=NULL, ratings=ratings))
}

#' Read and validate the demographics data file.
#'
#' Reads \code{Demographics.csv} and aligns it to the raters in \code{ratings}.
#' Columns with no variation (all identical values) or where every value is
#' distinct (e.g. free-text IDs) are silently dropped. Categorical columns are
#' one-hot encoded and numeric columns are split into above/below-median binary
#' indicators to form a \code{cohorts} matrix used for subgroup comparisons.
#' @param filename Path to the demographics CSV file (\code{Demographics.csv}).
#'   The first column must be named \code{RaterID} (or \code{UserID}/
#'   \code{SorterID}) and must match the RaterIDs in \code{ratings}.
#' @param ratings The validated ratings data frame returned by
#'   \code{read_ratings()}, used to align demographic rows to raters.
#' @param enc File encoding passed to \code{read.csv} (default \code{"UTF-8"}).
#' @param sep Column separator passed to \code{read.csv} (default \code{","}).
#' @return A list with elements: \code{issues} (NULL on success, otherwise an
#'   error string that halts loading), \code{warnings} (NULL if all raters have
#'   demographic data, otherwise a message listing RaterIDs found in
#'   \code{Ratings.csv} but absent from \code{Demographics.csv} — these raters
#'   are included in the \code{allRaters} cohort with NA for all demographic
#'   variables), \code{demographics.dat} (the raw demographics data frame),
#'   \code{dframe} (cleaned and aligned demographic variables, with NA rows for
#'   missing raters), and \code{cohorts} (binary indicator matrix with one
#'   column per subgroup; missing raters have 0 for all subgroup columns).
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

  # Pre-check for duplicate RaterIDs in the demographics file (done once, not per row)
  demo_ids <- as.character(demographics.dat$RaterID)
  dup_ids <- unique(demo_ids[duplicated(demo_ids)])
  if (length(dup_ids) > 0)
    return(list(issues=paste("read_demographics:", filename,
                             "RaterID(s)", paste(dup_ids, collapse=", "),
                             "appear more than once.\n"),
                demographics.dat = NULL))

  # Identify raters in Ratings.csv that have no row in Demographics.csv.
  # These are a warning (not an error): they will be assigned NA for all
  # demographic variables and will only appear in the 'allRaters' cohort.
  ratings_ids <- as.character(unique(ratings$RaterID))
  missing_raters <- setdiff(ratings_ids, demo_ids)

  # Initialise dframe with NA so missing raters stay NA rather than 0/empty
  dframe <- data.frame(matrix(NA, ncol=ncol(demographics.dat)-1, nrow=nrow(ratings)))
  cohorts <- data.frame("allRaters"=matrix(1, nrow=nrow(ratings), ncol=1))
  skipcols <- c()
  for(colnum in 2:ncol(demographics.dat)) {
    if(length(unique(demographics.dat[,colnum])) == 1) { # no variation in the column
      skipcols <- c(skipcols, colnum)
      next
    }
    if (inherits(demographics.dat[,colnum], "factor")) {
      if(length(unique(demographics.dat[,colnum])) == nrow(demographics.dat)) { # all distinct
        skipcols <- c(skipcols, colnum)
        next
      }
      dframe[,colnum-1] <- factor(rep(NA_character_, nrow(ratings)),
                                  levels=levels(demographics.dat[,colnum]))
    } else {
      dframe[,colnum-1] <- rep(NA_real_, nrow(ratings))
    }
    for (i in 1:nrow(ratings)) {
      matched <- which(demo_ids == as.character(ratings$RaterID[i]))
      if (length(matched) == 0) next  # missing rater — leave as NA
      dframe[i, colnum-1] <- demographics.dat[matched, colnum]
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
    if (inherits(dframe[,colnum], "factor")) {
      # na.action=na.pass keeps rows with NA; fill resulting NAs with 0
      # so missing raters don't belong to any subgroup
      tmp <- model.matrix(~ dframe[,colnum] - 1, na.action = na.pass)
      tmp[is.na(tmp)] <- 0
      colnames(tmp) <- sprintf("%s_%s", colnames(dframe)[colnum],
                               levels(dframe[,colnum]))
      cohorts <- data.frame(cohorts, tmp)
    } else { # type=numeric is assumed
      tmp <- model.matrix(~ (dframe[,colnum] > median(dframe[,colnum], na.rm=TRUE)) - 1,
                          na.action = na.pass)
      tmp[is.na(tmp)] <- 0
      colnames(tmp) <- paste0(colnames(dframe)[colnum], c("_H","_L"))
      cohorts <- data.frame(cohorts, tmp)
    }
  }

  warn_msg <- NULL
  if (length(missing_raters) > 0)
    warn_msg <- sprintf(
      "The following RaterID(s) appear in Ratings.csv but have no demographic data: %s\nThese raters will be included in the 'allRaters' cohort only; subgroup analyses will exclude them.",
      paste(missing_raters, collapse=", ")
    )

  return(list(issues=NULL, warnings=warn_msg, demographics.dat=demographics.dat,
              dframe=dframe, cohorts=cohorts))
}

#' Print the expected input file formats for an RCMap project.
#'
#' Prints a checklist of required and optional project files with their expected
#' column layouts. Useful for diagnosing input errors.
#' @export
print_input_checklist <- function() {
  cat("Required project files and expected columns:\n")
  cat(" - Statements.csv : two columns [StatementID, Statement]\n")
  cat("     * StatementID should be numeric and unique (e.g. 1,2,3...)\n")
  cat(" - SortedCards.csv : first column sorter ID, second column pile label, columns 3.. contain card numbers (integers referencing StatementID)\n")
  cat(" - Ratings.csv : first two columns must be RaterID and StatementID, remaining columns are numeric rating variables (1..ratingscale)\n")
  cat(" - Demographics.csv : first column RaterID matching Ratings.csv, other columns demographic variables (factors or numeric)\n")
  cat("Optional files:\n - config.txt (optional)\n - Weights.csv (optional)\n - CMapSession.RData / dictionaries.RData (previous session)\n")
  invisible(NULL)
}

#' Read the card-sorting data files and generate the MDS data.
#'
#' Loads all project input files from \code{dataDir}, validates them, builds
#' per-sorter adjacency matrices, runs multidimensional scaling (MDS), computes
#' Euclidean and hyperbolic 2-D distances, and performs the leave-one-out
#' Jaccard index analysis. If a previously saved session (\code{CMapSession.RData})
#' is found in \code{dataDir}, the user is prompted to reuse it.
#'
#' @param dataDir Path to the project folder containing \code{Statements.csv},
#'   \code{SortedCards.csv}, \code{Ratings.csv}, and \code{Demographics.csv}.
#'   An optional \code{config.txt} (for \code{ratingscale}) and
#'   \code{Weights.csv} (per-sorter weights) are also read if present.
#' @param enc File encoding passed to \code{read.csv} (default \code{"UTF-8"}).
#' @param sep Column separator passed to \code{read.csv} (default \code{","}).
#' @return A named list (\code{cmapdat}) with elements: \code{n.indiv} (number
#'   of sorters), \code{cardNames} (statements data frame), \code{adj.mat}
#'   (list of adjacency matrices), \code{DS} (SxS dissimilarity matrix),
#'   \code{D2e}/\code{D2h} (Euclidean and hyperbolic 2-D distance matrices),
#'   \code{x}/\code{y} (MDS coordinates), \code{stress} (MDS stress value),
#'   \code{ratings}, \code{demographics}, \code{cohorts}, \code{eta_ijk_e}/
#'   \code{eta_ijk_h} (Jaccard index arrays), and settings for the interactive
#'   menu. The list is also saved to \code{CMapSession.RData} in \code{dataDir}.
#' @details Uses \code{<<-} to set \code{clustMethod} and \code{cmapdat} in the
#'   calling environment, as required by the menu-driven interface. After
#'   reading pile labels, similar labels (across sorters) are merged using
#'   Jaro-Winkler fuzzy matching (threshold 0.15) and a \code{canonical} label
#'   is assigned to each group; the mapping is stored in \code{label_dict} and
#'   saved to \code{dictionaries.RData}. Run
#'   \code{print_input_checklist()} to review the expected file formats.
#' @importFrom smacof mds
initCMap <- function(dataDir, enc="UTF-8", sep=",") {
  # Check that the four input files are all there, and are formatted correctly:

  if (is.null(dataDir) || length(dataDir) != 1 || is.na(dataDir) || !nzchar(dataDir))
    stop("initCMap: 'dataDir' must be a non-empty string pointing to the project directory.")
  if (!dir.exists(dataDir))
    stop(sprintf("initCMap: data directory '%s' does not exist.", dataDir))

  required_files <- c("Statements.csv", "SortedCards.csv", "Ratings.csv", "Demographics.csv")
  missing_files <- required_files[!file.exists(file.path(dataDir, required_files))]
  if (length(missing_files) > 0) {
    stop(sprintf("initCMap: Missing required file(s) in '%s': %s.\nPlease ensure the project folder contains these CSV files (case-sensitive).",
                 dataDir, paste(missing_files, collapse=", ")))
  }

  savedfile <- paste0(dataDir,"/CMapSession.RData")
  dictfile <- paste0(dataDir,"/dictionaries.RData")
  if (file.exists(savedfile)) {
    cat(green("A previous project file was found.\n"))
    prmpt <- "Hit Enter to use saved project, or enter L to load the project from raw files: "
    usefile <- readline(prompt = prmpt)
    if (toupper(usefile) == "") {
      load(savedfile)
      if (file.exists(dictfile)) load(dictfile)
      return(cmapdat)
    }
  }

  cat(green("Loading project files. Please wait.\n\n"))
  cmapdat <<- list()
  conffile <- paste0(dataDir, "/config.txt")

  # All recognised config keys with their defaults and inline help text.
  config_defaults <- list(
    ratingscale            = list(val = "5",        type = "numeric",
      comment = "number of Likert scale points (e.g. 5, 7)"),
    clust_method           = list(val = "ward.D2",  type = "character",
      comment = "hierarchical clustering method (ward.D2, ward.D, single, complete, average, mcquitty, median, centroid)"),
    dist_metric            = list(val = "Euclidean", type = "character",
      comment = "initial distance metric for plots (Euclidean, Hyperbolic)"),
    color_scheme           = list(val = "rcmap",    type = "character",
      comment = "initial color scheme for plots (rcmap, rainbow)"),
    n_clusters             = list(val = "3",        type = "integer",
      comment = "initial number of clusters"),
    mds_seed               = list(val = "154204",   type = "integer",
      comment = "random seed for MDS distance-matrix jitter (reproducibility)"),
    splithalf_seed         = list(val = "23456",    type = "integer",
      comment = "random seed for split-half reliability analysis"),
    splithalf_B            = list(val = "20",       type = "integer",
      comment = "number of split-half replications"),
    fuzzy_label_threshold  = list(val = "0.15",     type = "numeric",
      comment = "Jaro-Winkler distance threshold for pile label fuzzy matching (0-1, lower = stricter)"),
    jaccard_threshold      = list(val = "0.3",      type = "numeric",
      comment = "Jaccard index values below this are flagged as unstable in the stability plot")
  )

  # Parse existing config file (skip comment lines and blank lines).
  cfg <- list()
  if (file.exists(conffile)) {
    for (line in readLines(conffile)) {
      line <- trimws(line)
      if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
      key <- trimws(sub("=.*", "", line))
      val <- trimws(sub("^[^=]+=", "", line))
      if (nzchar(key) && nzchar(val)) cfg[[key]] <- val
    }
  }

  # Write defaults for any keys absent from the file.
  missing_keys <- setdiff(names(config_defaults), names(cfg))
  if (length(missing_keys) > 0) {
    header <- if (!file.exists(conffile))
      c("# RCMap project configuration file",
        "# Edit the values below to override defaults.", "")
    else
      c("", "# The following settings were added automatically with default values.")
    lines_to_add <- header
    for (key in missing_keys) {
      lines_to_add <- c(lines_to_add,
                        paste0("# ", key, ": ", config_defaults[[key]]$comment),
                        paste0(key, "=", config_defaults[[key]]$val),
                        "")
    }
    cat(paste(lines_to_add, collapse = "\n"),
        file = conffile, append = file.exists(conffile))
  }

  # Merge: user-supplied values override defaults.
  cfg_merged <- config_defaults
  for (key in names(cfg)) {
    if (key %in% names(cfg_merged)) cfg_merged[[key]]$val <- cfg[[key]]
  }

  # Helper: extract a typed scalar from the merged config.
  cfg_get <- function(key) {
    entry <- cfg_merged[[key]]
    v <- entry$val
    switch(entry$type,
           numeric   = { n <- suppressWarnings(as.numeric(v));   if (is.na(n)) config_defaults[[key]]$val else n },
           integer   = { n <- suppressWarnings(as.integer(v));   if (is.na(n)) config_defaults[[key]]$val else n },
           character = v)
  }

  ratingscale           <- cfg_get("ratingscale")
  clustMethod           <- cfg_get("clust_method")
  clustMethod           <<- clustMethod          # global used by update_clusters()
  distmetric_default    <- cfg_get("dist_metric")
  clrscm_default        <- cfg_get("color_scheme")
  nclust_default        <- cfg_get("n_clusters")
  mds_seed              <- cfg_get("mds_seed")
  splithalf_seed        <- cfg_get("splithalf_seed")
  splithalf_B           <- cfg_get("splithalf_B")
  fuzzy_label_threshold <- cfg_get("fuzzy_label_threshold")
  jaccard_threshold     <- cfg_get("jaccard_threshold")

  # Statements:
  cardNames <- tryCatch(
    read_statements(paste0(dataDir, "/Statements.csv"), enc, sep),
    error = function(e) {
      message("initCMap: Error reading 'Statements.csv': ", e$message)
      message("Run print_input_checklist() for the expected file/column formats.")
      stop(e)
    }
  )
  if (!is.null(cardNames$issues) && length(cardNames$issues) > 0)
    stop(sprintf("initCMap: Problems in 'Statements.csv': %s\nRun print_input_checklist() for expected formats.", cardNames$issues))
  # If no issues - keep only the data frame
  cardNames <- cardNames$cardNames
  if (is.null(cardNames) || nrow(cardNames) == 0)
    stop("initCMap: 'Statements.csv' appears empty or malformed. Ensure it contains statement IDs and text. Run print_input_checklist() for details.")
  card_name_dict <- cardNames
  # Validate StatementID is numeric and unique
  stmt_ids <- suppressWarnings(as.numeric(as.character(cardNames$StatementID)))
  if (any(is.na(stmt_ids)))
    {
      bad_idx <- which(is.na(stmt_ids))
      bad_vals <- as.character(cardNames$StatementID[bad_idx])
      stop(sprintf("initCMap: 'Statements.csv' StatementID column must be numeric. Offending values (up to 5): %s (rows: %s). Run print_input_checklist() for expected format.",
                   paste(head(bad_vals, 5), collapse=","), paste(head(bad_idx,5), collapse=",")))
    }
  if (any(duplicated(stmt_ids)))
    {
      dupvals <- unique(stmt_ids[duplicated(stmt_ids)])
      stop(sprintf("initCMap: Duplicate StatementID values found in 'Statements.csv'. Examples (up to 5): %s. StatementID must be unique.", paste(head(dupvals,5), collapse=",")))
    }
  cardNames$StatementID <- stmt_ids

  # Piles:
  cardDat <- tryCatch(
    read_piles(paste0(dataDir, "/SortedCards.csv"), nrow(cardNames), enc, sep),
    error = function(e) {
      message("initCMap: Error reading 'SortedCards.csv': ", e$message)
      message("Run print_input_checklist() for the expected file/column formats.")
      stop(e)
    }
  )
  if (!is.null(cardDat$issues) && length(cardDat$issues) > 0)
    stop(sprintf("initCMap: Problems in 'SortedCards.csv': %s\nRun print_input_checklist() for expected formats.", cardDat$issues))
  # If no issues - keep only the data frame
  if (is.null(cardDat$cardDat)) stop("initCMap: 'SortedCards.csv' did not return expected structure from read_piles().")
  label_dict <- data.frame(pileLabels=cardDat$pileLabels[,1],
                           lcpilelabels=cardDat$pileLabelsLC,
                           stringsAsFactors=FALSE)
  pileLabels <- cardDat$pileLabels
  cardDat <- cardDat$cardDat

  # Fuzzy pile label matching: group similar labels across sorters and assign
  # a canonical form to each cluster.  Labels are compared using Jaro-Winkler
  # distance (p=0.1), which handles transpositions, abbreviations, and minor
  # typos.  Labels whose maximum pairwise distance within a group is <= 0.15
  # are merged; the canonical label is the most frequently used lowercase label
  # in that group.
  unique_lc <- unique(label_dict$lcpilelabels)
  if (length(unique_lc) > 1) {
    dist_mat  <- stringdist::stringdistmatrix(unique_lc, unique_lc,
                                              method = "jw", p = 0.1)
    hc        <- hclust(as.dist(dist_mat), method = "complete")
    clusters  <- cutree(hc, h = fuzzy_label_threshold)
    lc_to_canonical <- setNames(unique_lc, unique_lc)  # default: identity
    merged_info <- character(0)
    for (cid in unique(clusters)) {
      members <- unique_lc[clusters == cid]
      if (length(members) > 1) {
        freq_tab  <- sort(table(label_dict$lcpilelabels[
                                  label_dict$lcpilelabels %in% members]),
                          decreasing = TRUE)
        canonical <- names(freq_tab)[1]
        lc_to_canonical[members] <- canonical
        merged_info <- c(merged_info,
                         sprintf("  {%s} -> '%s'",
                                 paste(members, collapse = ", "), canonical))
      }
    }
    label_dict$canonical <- lc_to_canonical[label_dict$lcpilelabels]
    if (length(merged_info) > 0) {
      cat(yellow("Note: The following pile labels were merged based on similarity:\n"))
      cat(yellow(paste(merged_info, collapse = "\n")), "\n")
    }
  } else {
    label_dict$canonical <- label_dict$lcpilelabels
  }
  # Replace pile labels in pileLabels with canonical forms so that
  # label-frequency tables across sorters aggregate correctly.
  pileLabels$pileLabel <- label_dict$canonical
  # Validate that card numbers referenced in SortedCards.csv exist in Statements.csv
  cards_used <- tryCatch({
    vals <- suppressWarnings(as.numeric(na.exclude(unique(as.vector(unlist(cardDat[ ,3:ncol(cardDat)]))))))
    sort(vals)
  }, error = function(e) NA)
  if (!all(is.na(cards_used))) {
    invalid_cards <- setdiff(cards_used, cardNames$StatementID)
    if (length(invalid_cards) > 0) {
      # find sample rows in SortedCards that reference invalid cards
      bad_row_idx <- which(apply(cardDat[, 3:ncol(cardDat), drop=FALSE], 1, function(r) any(na.omit(suppressWarnings(as.numeric(r))) %in% invalid_cards)))
      sample_rows <- head(bad_row_idx, 5)
      samples <- sapply(sample_rows, function(i) paste(na.omit(as.character(cardDat[i, 3:ncol(cardDat)])), collapse=","))
      stop(sprintf("initCMap: SortedCards.csv references card ID(s) not present in Statements.csv: %s. Example offending rows (up to 5): %s (row numbers: %s). Run print_input_checklist() to review expected formats.",
                   paste(head(invalid_cards,5), collapse=","), paste(samples, collapse=" | "), paste(sample_rows, collapse=",")))
    }
  }

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
  wgts <- rep(1, n.indiv)
  wgtfile <- paste0(dataDir,"/Weights.csv")
  if (file.exists(wgtfile)) {
    wgtdat <- tryCatch(
      read.csv(wgtfile),
      error = function(e) stop(sprintf("initCMap: Error reading 'Weights.csv': %s", e$message))
    )
    if (!"Weight" %in% colnames(wgtdat))
      stop("initCMap: 'Weights.csv' must contain a column named 'Weight'.")
    if (nrow(wgtdat) != n.indiv)
      stop(sprintf("initCMap: 'Weights.csv' has %d row(s) but there are %d sorter(s). Each sorter must have exactly one weight.", nrow(wgtdat), n.indiv))
    wgt_vals <- suppressWarnings(as.numeric(wgtdat$Weight))
    if (any(is.na(wgt_vals)))
      stop("initCMap: 'Weights.csv' Weight column contains missing or non-numeric values.")
    if (any(wgt_vals <= 0))
      stop("initCMap: 'Weights.csv' Weight column must contain strictly positive values.")
    wgts <- n.indiv * wgt_vals / sum(wgt_vals)
  }
  for (i in 1:n.indiv) {
    adj.mat[[i]] <- adj.mat[[i]]*wgts[i]
  }
  # multidimensional scaling, and clustering:
  DS <- distanceMatrix(adj.mat, seed = mds_seed)
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
  ratings <- tryCatch(
    read_ratings(paste0(dataDir, "/Ratings.csv"), cardNames$StatementID, ratingscale, enc, sep),
    error = function(e) {
      message("initCMap: Error reading 'Ratings.csv': ", e$message)
      message("Run print_input_checklist() for the expected file/column formats.")
      stop(e)
    }
  )
  if (!is.null(ratings$issues) && length(ratings$issues) > 0)
    stop(sprintf("initCMap: Problems in 'Ratings.csv': %s\nRun print_input_checklist() for expected formats.", ratings$issues))
  ratings <- ratings$ratings
  # Validate that rating columns (3..) are numeric and within range
  if (ncol(ratings) >= 3) {
    for (ci in 3:ncol(ratings)) {
      colvals <- ratings[,ci]
      conv <- suppressWarnings(as.numeric(as.character(colvals)))
      # If there are non-missing entries that become NA after coercion -> invalid
      if (any(!is.na(colvals) & is.na(conv))) {
        bad_rows <- which(!is.na(colvals) & is.na(conv))
        sample_bad <- head(bad_rows, 5)
        sample_strs <- paste(sapply(sample_bad, function(r) sprintf("row %d: RaterID=%s, StatementID=%s, value='%s'", r, as.character(ratings$RaterID[r]), as.character(ratings$StatementID[r]), as.character(colvals[r]))), collapse=" | ")
        stop(sprintf("initCMap: Ratings column '%s' contains non-numeric values. Examples: %s. Run print_input_checklist() for expected formats.", colnames(ratings)[ci], sample_strs))
      }
      # replace in-place with numeric conversion (preserve NAs)
      ratings[,ci] <- conv
      # check range
      if (any(!is.na(conv) & (conv < 1 | conv > ratingscale))) {
        bad_rows <- which(!is.na(conv) & (conv < 1 | conv > ratingscale))
        sample_bad <- head(bad_rows, 5)
        sample_strs <- paste(sapply(sample_bad, function(r) sprintf("row %d: RaterID=%s, StatementID=%s, value=%s", r, as.character(ratings$RaterID[r]), as.character(ratings$StatementID[r]), as.character(ratings[r,ci]))), collapse=" | ")
        stop(sprintf("initCMap: Ratings column '%s' contains values outside 1..%d. Examples: %s", colnames(ratings)[ci], ratingscale, sample_strs))
      }
    }
  }
  for (i in 1:n.indiv) {
    # Match by original sorter ID, not sequential index i, so weights are
    # applied correctly when RaterIDs are non-consecutive or non-numeric.
    rownums <- which(as.character(ratings$RaterID) == as.character(sorters_dict$orig_ID[i]))
    ratings[rownums, -c(1,2)] <- ratings[rownums, -c(1,2)]*wgts[i]
  }
  # Demographic data:
  demographics.dat <- tryCatch(
    read_demographics(paste0(dataDir, "/Demographics.csv"), ratings, enc, sep),
    error = function(e) stop(sprintf("initCMap: Error reading 'Demographics.csv': %s", e$message))
  )
  if (!is.null(demographics.dat$issues) && length(demographics.dat$issues) > 0)
    stop(sprintf("initCMap: Problems in 'Demographics.csv': %s", demographics.dat$issues))

  # Case 1 (warning): raters in Ratings.csv with no demographic row.
  # read_demographics already filled their rows with NA and built cohorts
  # treating them as 'allRaters' only. Inform the user and allow proceeding.
  if (!is.null(demographics.dat$warnings)) {
    cat(yellow(bold("\nWarning: ")) %+% yellow(demographics.dat$warnings) %+% "\n")
    readline("Press Enter to continue. ")
  }

  dframe <- demographics.dat$dframe
  cohorts <- demographics.dat$cohorts
  demographics.dat <- demographics.dat$demographics.dat

  # Case 2 (exclusion): raters in Demographics.csv with no data in Ratings.csv.
  # These raters cannot contribute to any analysis; exclude them from
  # demographics.dat so they don't appear in reports or summaries.
  ratings_raters <- as.character(unique(ratings$RaterID))
  demo_raters    <- as.character(unique(demographics.dat$RaterID))
  extra_in_demo  <- setdiff(demo_raters, ratings_raters)
  if (length(extra_in_demo) > 0) {
    cat(yellow(bold("\nWarning: ")) %+%
        yellow(sprintf("%d rater(s) in Demographics.csv have no rating data and will be excluded: %s\n",
                       length(extra_in_demo), paste(extra_in_demo, collapse=", "))))
    demographics.dat <- demographics.dat[
      !as.character(demographics.dat$RaterID) %in% extra_in_demo, , drop=FALSE]
  }

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
  cmapdat <- list(n.indiv=n.indiv, cardNames=cardNames,
                  adj.mat=adj.mat, DS=DS,
                  D2e=D2e, D2h=D2h, x=x, y=y, stress=stress, ratings=ratings,
                  clustname=list(), pileLabels=pileLabels, issues=issues,
                  demographics=demographics.dat, ratingDemographics=dframe,
                  clustMethod=clustMethod,
                  distmetric=distmetric_default, cohortCol=1,
                  nclust=nclust_default, splhalf=list(), cohorts=cohorts,
                  pcoordChoice=c(1, 2), clusMember=cmapdat$clusMember,
                  clrscm=clrscm_default, eta_ijk_e=eta_ijk_e, eta_ijk_h=eta_ijk_h,
                  ratingscale=ratingscale, dataDir=dataDir,
                  sorters_dict=sorters_dict, label_dict=label_dict,
                  mds_seed=mds_seed, splithalf_seed=splithalf_seed,
                  splithalf_B=splithalf_B, jaccard_threshold=jaccard_threshold)
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


#' Plot the 2-D MDS representation of the pile-sorting data.
#'
#' @param metric Controls what is displayed. \code{NULL} (default) shows a plain
#'   point map with all statements in blue. \code{"MPindex"} colors and sizes
#'   points by the misplacement index (Jaccard-based bridging/anchoring score).
#'   An integer \code{k} (offset by 2 into the ratings data frame) colors and
#'   sizes points by the mean rating for that variable.
#' @param tau The Jaccard index threshold used to compute the misplacement
#'   index proportion (default=0.3). Only used when \code{metric="MPindex"}.
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
      sz <- (sz-1)/(cmapdat$ratingscale-1)
      lgd <- seq(0, 1, length=length(ter.cols))
      lgdkey <- lgd*4*(cmapdat$ratingscale/5) + 1
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

#' Plot the clusters in the 2-D representation of the pile-sorting data.
#'
#' Can show the clusters either as convex hull polygons or as rays connecting
#' each statement to its cluster centroid.
#' @param type Display style: \code{"rays"} (default) draws lines from each
#'   point to the cluster center; \code{"polygon"} draws a filled convex hull
#'   around each cluster.
#' @param metric Controls cluster-center sizing. \code{NULL} (default) shows
#'   cluster centers as equal-sized markers. An integer \code{k} (offset by 2
#'   into the ratings data frame) sizes the cluster-center marker by the mean
#'   rating for that variable in the cluster.
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
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
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
        sz[i] <- rtgsmm[nrow(rtgsmm),
                        which(colnames(rtgsmm) ==
                                paste0(colnames(ratings)[metric+2],"_Mean"))]
      }
      #sz <- 1+2*(sz-1)/(ratingscale-1) # 1+2*(sz-1)/4
      #sz <- round(sz, digits = 1) # to be between 1 (for 1) and 2 (for 5)
      sz <- 3*round(sz/cmapdat$ratingscale, digits = 1)
      #lgd <- seq(1, 3, by=0.5)
      #lgdkey <- (lgd-1)*(ratingscale/5)+1
      lgd <- 3*seq(0, 4/5, length=5)
      lgdkey <- lgd*cmapdat$ratingscale/2
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
    for (i in 1:length(cardNames[,1])) {
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


#' Recompute cluster membership and update the global session object.
#'
#' Runs hierarchical clustering on the currently selected distance matrix
#' (\code{D2e} or \code{D2h} depending on \code{cmapdat$distmetric}) with
#' \code{cmapdat$nclust} clusters, and stores the result in
#' \code{cmapdat$clusMember} via \code{<<-}.
update_clusters <- function() {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
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
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
    clus <- cutree(fit.Clust, nclust)
    cmapdat$clusMember <<- clus
    cnames <- clusterNames()
    # Compute legend layout: how many columns fit across the figure width
    item_width_in <- (max(nchar(cnames)) + 3) * par("cin")[1] * 0.8
    ncol_leg <- max(1, floor(par("fin")[1] / item_width_in))
    nrow_leg <- ceiling(length(cnames) / ncol_leg)
    # Expand bottom margin to fit the legend rows plus a gap line
    bot_mar <- nrow_leg * 2 + 4
    old_mar <- par(mar = c(bot_mar, par("mar")[2:4]))
    on.exit(par(old_mar), add = TRUE)
    if(phylo) {
      X <- as.phylo(fit.Clust)
      edge.clus <- sapply(seq_len(cmapdat$nclust),
                          function(i) max(which(X$edge[,2] %in% which(clus==i))))
      ord       <- order(edge.clus)
      edge.clus <- c(min(edge.clus), diff(sort(edge.clus)))
      edge.clus <- rep(ord, edge.clus)
      plot(X, type = "fan",
           tip.color = cols[clus], edge.color = cols[edge.clus],
           label.offset = 0, no.margin = FALSE, cex = 0.70)
    } else {
      clusDendro <- dendrapply(as.dendrogram(fit.Clust), colLab)
      plot(clusDendro, main = "")
    }
    # Place legend in the bottom margin (below the plot region)
    inset_down <- (nrow_leg * par("cin")[2] * 2 + 0.1) / par("pin")[2]
    legend("bottom", legend = cnames, fill = cols, ncol = ncol_leg,
           bty = "n", cex = 0.8, xpd = NA, inset = c(0, -inset_down))
  })
}

#' Show statement ratings as a dotchart.
#'
#' Displays mean ratings per statement grouped by cluster, sorted in ascending
#' order, with statements color-coded by their cluster membership.
#' @param metric An integer selecting which rating variable to display (column
#'   index into the ratings data frame, offset by 2). Returns \code{NULL}
#'   invisibly if \code{metric} is 0.
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
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
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
               color = cols[groups[order(y, decreasing = FALSE)]],
               xlim = c(1,ratingscale+1),
               cex = 0.6,  pch = 19, xlab = ttl, frame.plot = FALSE,
               xaxt = "n")
      abline(v=1:(ratingscale-1), lwd=2, lty=2, col="gray66")
      text((1:ratingscale), rep(max(gpos)+1, ratingscale), 1:ratingscale, cex=0.8, font=2)
      rect(ratingscale+.05,0,ratingscale+2,max(gpos)+2, col = "white", border="white")
      rect(0,0,0.95,max(gpos)+2, col = "white", border="white")
      text(rep(ratingscale+.2, length(cardNames[,1])), gpos, clusterNames(), cex=0.8,
           col=cols, font=2, pos = 4)
    }
  })
}

#' Show average ratings by cluster as a bar chart.
#'
#' Displays mean ratings (with standard error bars) per cluster for a selected
#' rating variable. Prompts the user to choose bar ordering (alphabetical or by
#' height) before plotting.
#' @param metric An integer selecting which rating variable to display (column
#'   index into the ratings data frame, offset by 2). Returns \code{NULL}
#'   invisibly if \code{metric} is 0.
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
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
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
    bot_mar <- max(5, ceiling(max(nchar(cname)) * 0.7) + 2)
    old_mar <- par(mar = c(bot_mar, par("mar")[2:4]))
    on.exit(par(old_mar), add = TRUE)
    bp <- barplot(MN, main=ttl,col=barcols,ylim=c(0,Yl),
                  names.arg=cname,cex.names=0.8,las=3,
                  width=rep(0.8,nclust),space=0.2)
    segments(bp, MN - SEM, bp, MN + SEM, lwd=2, col="grey66")
    segments(bp - 0.1, MN - SEM, bp + 0.1, MN - SEM, lwd=2, col="grey66")
    segments(bp - 0.1, MN + SEM, bp + 0.1, MN + SEM, lwd=2, col="grey66")
  })
}

#' Show a parallel coordinates plot of cluster mean ratings.
#'
#' Plots mean ratings per cluster across multiple rating-variable / cohort
#' combinations, with one axis position per combination and one colored line per
#' cluster. The variable / cohort combinations are taken from
#' \code{cmapdat$pcoordChoice}, which the user selects before this function is
#' called. Cluster names are shown at the right edge of the plot.
#' @export
showParallelCoordinates <- function() {
  with(cmapdat,{
    choices <- expand.grid(colnames(ratings)[-(1:2)], colnames(cohorts))
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    cnms <- clusterNames()
    # Bottom margin to fit horizontal legend
    item_width_in <- (max(nchar(cnms)) + 3) * par("cin")[1] * 0.8
    ncol_leg <- max(1, floor(par("fin")[1] / item_width_in))
    nrow_leg <- ceiling(length(cnms) / ncol_leg)
    bot_mar <- nrow_leg * 2 + 4
    old_mar <- par(mar = c(bot_mar, par("mar")[2:4]))
    on.exit(par(old_mar), add = TRUE)
    MN <- matrix(0, nrow=nclust, ncol=length(pcoordChoice))
    for (jj in seq_along(pcoordChoice)) {
      vars <- unlist(strsplit(pcoordChoice[jj],"/"))
      metric <- which(colnames(ratings) == vars[1])
      cohortCol <- which(colnames(cohorts) == vars[2])
      ratingstmp <- ratings[which(cohorts[,cohortCol] == 1),]
      for (i in seq_along(unique(groups))) {
        cardsInCluster <- names(which(groups == i))
        cardsInCluster <- which(ratingstmp$StatementID %in% cardsInCluster)
        rtgsmm <- ratingSummary(ratingstmp[cardsInCluster,])
        rownames(rtgsmm) <- c(rownames(DS)[which(groups == i)], "total")
        MN[i, jj] <- rtgsmm[nrow(rtgsmm),which(colnames(rtgsmm) == paste0(colnames(ratings)[metric],"_Mean"))]
      }
      if (jj == 1) {
        plot(rep(1, nclust), MN[,1], ylim=c(0.5, ratingscale), xlim=c(0.2,length(pcoordChoice)+0.8),
             main="", col=cols, axes=FALSE, xlab="", ylab="Rating", pch=19, cex=0.3)
        axis(2)
        grid()
      } else {
        points(rep(jj, nclust), MN[,jj], col=cols, pch=19, cex=0.3)
        for (k in seq_len(nclust))
          lines(c(jj-1,jj), c(MN[k,jj-1], MN[k,jj]), col=cols[k])
      }
      text(jj, 0.8, cmapdat$pcoordChoice[jj], cex=0.6)
    }
    inset_down <- (nrow_leg * par("cin")[2] * 2 + 0.1) / par("pin")[2]
    legend("bottom", legend = cnms, fill = cols, ncol = ncol_leg,
           bty = "n", cex = 0.8, xpd = NA, inset = c(0, -inset_down))
  })
}


#' A Go-Zone plot.
#'
#' Plots two rating variables against each other (one per axis), with each
#' statement labeled by its ID and colored by cluster. Reference lines are drawn
#' at the overall mean of each variable, dividing the plot into four quadrants.
#' Kendall's tau and its p-value are shown in the title. The two rating
#' variables to compare are taken from \code{cmapdat$pcoordChoice}.
#' @export
showGoZone <- function() {
  #  with(cmapdat,{
  choices <- expand.grid(colnames(cmapdat$ratings)[-(1:2)],
                         colnames(cmapdat$cohorts))
  distM <- cmapdat$D2e
  if(cmapdat$distmetric == "Hyperbolic")
    distM <- cmapdat$D2h
  cols <- clusterCols(cmapdat$nclust, cmapdat$clrscm)
  cnms <- clusterNames()
  item_width_in <- (max(nchar(cnms)) + 3) * par("cin")[1] * 0.8
  ncol_leg <- max(1, floor(par("fin")[1] / item_width_in))
  nrow_leg <- ceiling(length(cnms) / ncol_leg)
  bot_mar <- nrow_leg * 2 + 4
  old_mar <- par(mar = c(bot_mar, par("mar")[2:4]))
  on.exit(par(old_mar), add = TRUE)
  fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
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
  plot( MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2],
        ylim=c(1,cmapdat$ratingscale), xlim=c(1,cmapdat$ratingscale),
        main=Kendall, col=0, axes=FALSE, xlab=cmapdat$pcoordChoice[1],
        ylab=cmapdat$pcoordChoice[2], pch=19, cex=0.6, cex.main=1)
  text(MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2],
       rownames(cmapdat$DS), cex=0.7,
       col=cols[groups], font=2)
  axis(1, at = 1:cmapdat$ratingscale); axis(2, at = 1:cmapdat$ratingscale)
  grid()
  rect(xleft = 1, ybottom = 1, xright = MN[nrow(rtgsmm),1],
       ytop = MN[nrow(rtgsmm),2], col=rgb(0.5, 0, 0.5, .1),
       border = "white")
  rect(xleft = 1, ybottom = MN[nrow(rtgsmm),2], xright = MN[nrow(rtgsmm),1],
       ytop = cmapdat$ratingscale, col=rgb(0.5, 0.5, 0, .1),
       border = "white")
  rect(xleft = MN[nrow(rtgsmm),1], ybottom = 1, xright = cmapdat$ratingscale,
       ytop = MN[nrow(rtgsmm),2], col=rgb(0.5, 0, 0, .1),
       border = "white")
  rect(xleft = MN[nrow(rtgsmm),1], ybottom = MN[nrow(rtgsmm),2], xright = cmapdat$ratingscale,
       ytop = cmapdat$ratingscale, col=rgb(0, 0.5, 0.5, .1),
       border = "white")
  abline(v=MN[nrow(rtgsmm),1], col="grey", lwd=2)
  abline(h=MN[nrow(rtgsmm),2], col="grey", lwd=2)
  inset_down <- (nrow_leg * par("cin")[2] * 2 + 0.1) / par("pin")[2]
  legend("bottom", legend = cnms, fill = cols, ncol = ncol_leg,
         bty = "n", cex = 0.8, xpd = NA, inset = c(0, -inset_down))

  #  })
}


#' Print a summary report about each sorter.
#'
#' For each sorter, prints the number of cards sorted and the number of piles
#' created. Pauses for user input before returning to the menu.
#' @export
sorterReport <- function() {
  with(cmapdat, {
    for (i in 1:length(adj.mat)) {
      cardsInPiles <- which(rowSums(adj.mat[[i]])>0)
      cardsSorted <- length(cardsInPiles)
      adjM <- (adj.mat[[i]][cardsInPiles, cardsInPiles] > 0)
      diag(adjM) <- 1
      noPiles <- length(unique(apply(adjM, 1, paste0, collapse="")))
      cat(paste0("Sorter ", sorters_dict$orig_ID[i], " sorted ", cardsSorted, " cards into ", noPiles," piles\n"))
    }
  })
  readline("Press any key to continue. ")
}

#' Print a summary report about the raters.
#'
#' Prints a statistical summary of all demographic variables for the raters.
#' Pauses for user input before returning to the menu.
#' @export
raterReport <- function() {
  with(cmapdat, {
    print(summary(cmapdat$demographics[,-1]))
  })
  readline("Press any key to continue. ")
}

#' Generate and save a per-statement summary report.
#'
#' For each statement, reports its cluster assignment, Jaccard index (clustering
#' stability), and rating summary statistics (N, mean, SD, min, max) for every
#' rating variable. Results are written to
#' \code{<dataDir>/output/StatementSummary<nclust>.csv} and the path is printed
#' to the console. Pauses for user input before returning to the menu.
#' @export
statementReport <- function() {
  with(cmapdat, {
    distM <- D2e
    eta_ij <- eta_ijk_e[[nclust+1]]
    if(distmetric == "Hyperbolic") {
      distM <- D2h
      eta_ij <- eta_ijk_h[[nclust+1]]
    }
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
    jidx <- colMeans(eta_ij)
    groups <- cutree(fit.Clust, k=nclust)
    retdf <- as.data.frame(matrix(0, ncol=5+5*(ncol(ratings)-2), nrow=nclust+nrow(DS)))
    clusLabels <- clusterNames()
    k <- 0
    for (i in 1:nclust) {
      cardsInCluster <- names(which(groups == i))
      jidx_i <- c(jidx[which(groups == i)], mean(jidx[which(groups == i)]))
      cNames <- c(cardNames[which(groups == i), "Statement"],"")
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      rtgsmm <- ratingSummary(ratings[cardsInCluster,])
      rtgsmm <- data.frame(cNames, c(names(cmapdat$x)[which(groups == i)],""),
                           rep(i, nrow(rtgsmm)), jidx_i, rtgsmm,
                           rep(clusLabels[i], nrow(rtgsmm)))
      colnames(rtgsmm)[1:4] <- c("Statement", "CardNo", "ClusterNo", "JI")
      colnames(rtgsmm)[ncol(rtgsmm)] <- "ClusterName"
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

#' Run one-way ANOVA comparing rating variables across clusters.
#'
#' For each rating variable, fits a one-way ANOVA with cluster membership as the
#' factor and prints the F-test table to the console. Results are also written
#' to \code{<dataDir>/output/ANOVA<nclust>.txt}. Pauses for user input before
#' returning to the menu.
#' @export
clusterANOVA <- function() {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    fn <- paste0(dataDir,"output/ANOVA", nclust,".txt")
    clusLabels <- clusterNames()
    Cluster <- rep(clusLabels[1], nrow(ratings))
    for (i in 1:length(groups)) {
      cardsInCluster <- names(which(groups == i))
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      Cluster[cardsInCluster] <- clusLabels[i]
    }
    Cluster <- factor(Cluster, levels = clusLabels)
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

#' Run Tukey HSD post-hoc tests for all pairwise cluster comparisons.
#'
#' For each rating variable, fits a one-way ANOVA and then runs Tukey's Honest
#' Significant Difference test for all pairs of clusters. Results are printed to
#' the console and written to \code{<dataDir>/output/Tukey<nclust>.txt}. Pauses
#' for user input before returning to the menu.
#' @export
clusterTUKEY <- function() {
  with(cmapdat,{
    distM <- D2e
    if(distmetric == "Hyperbolic")
      distM <- D2h
    ttl <- ""
    cols <- clusterCols(nclust, clrscm)
    fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    clusLabels <- clusterNames()
    Cluster <- rep(clusLabels[1], nrow(ratings))
    fn <- paste0(dataDir,"output/Tukey", nclust,".txt")
    for (i in 1:length(groups)) {
      cardsInCluster <- names(which(groups == i))
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      Cluster[cardsInCluster] <- clusLabels[i]
    }
    Cluster <- factor(Cluster, levels = clusLabels)

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

#' Clear the console and print the RCMap menu header.
#'
#' Clears the screen using a form-feed character (works in RStudio / RGui) and
#' a system \code{clear}/\code{cls} call (works in terminal sessions). Screen
#' clearing is skipped when the package-level flag
#' \code{.rcmap_clear_screen} is \code{FALSE}.
topLine <- function() {
  if (isTRUE(getOption("rcmap_clear_screen", default = TRUE))) {
    if (getOS() == "windows") {
      cat("\014")           # form-feed: works in RGui / RStudio on Windows
    } else {
      cat("\033[2J\033[H")  # ANSI: clear screen + cursor home (terminal/Mac)
      cat("\014")           # form-feed: also clears RStudio console pane
    }
  }
  cat(blue(underline(bold("\nRCMap command-line interface.\n"))))
}

#' Convert a menu selection index to its text label.
#'
#' Maps an integer menu selection to the corresponding string value for a given
#' settings category. When called with only \code{varType}, returns the full
#' vector of options for that category (useful for building menus).
#' @param varType Category of setting. One of \code{"clusteringMethod"},
#'   \code{"distance"}, \code{"MPparam"}, or \code{"colorScheme"}.
#' @param choice Integer index into the options vector for \code{varType}.
#' @return A character string (or vector when \code{choice} is missing) giving
#'   the label(s) for the selected setting.
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

#' Color dendrogram leaves by cluster membership.
#'
#' A helper passed to \code{dendrapply} that sets the label color and point
#' style of each leaf node to the color of its assigned cluster.
#' @param dn A dendrogram node (leaf or internal), as passed by \code{dendrapply}.
#' @return The same dendrogram node with updated \code{nodePar} attributes for
#'   leaf nodes; internal nodes are returned unchanged.
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
#' @importFrom utils menu read.csv select.list stack head
#' @importFrom graphics abline axis grid points rect text lines
#' @importFrom grDevices rgb
#' @importFrom stats as.dendrogram as.dist cor.test cutree dendrapply dist hclust is.leaf median model.matrix na.exclude na.omit na.pass quantile sd setNames interaction.plot
#' @importFrom cluster silhouette
#' @importFrom stringdist stringdistmatrix
#' @importFrom tcltk tk_choose.dir
#'
#' @param enc Encoding (default=UTF-8).
#' @param sep Column separator for csv files (default=,).
#' @param clear_screen Whether to clear the console before each menu
#'   (default=TRUE). Set to FALSE if your terminal does not support screen
#'   clearing, or to keep a scrollable history of menu interactions.
#' @export
#' @examples
#' \donttest{
#' #RCMapMenu()
#' }
RCMapMenu <- function(enc = "UTF-8", sep=",", clear_screen = TRUE) {
  old_opt <- options(rcmap_clear_screen = isTRUE(clear_screen))
  on.exit(options(old_opt), add = TRUE)
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
        tryCatch(
          loadRCMapData(enc, sep),
          error = function(e) {
            cat(red(bold("\nError loading project: ")) %+% red(conditionMessage(e)) %+% "\n")
            readline("Press Enter to return to the main menu. ")
          }
        )
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
                  cmapdat$dataDir, cmapdat$n.indiv, length(cmapdat$cardNames[,1]),cmapdat$issues))
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
        cmapdat$splhalf <<- splitHalf(cmapdat$adj.mat,
                                      B       = cmapdat$splithalf_B,
                                      disttype = cmapdat$distmetric,
                                      plotit  = TRUE,
                                      seed    = cmapdat$splithalf_seed)
      }
      if(summMenu == 2) {
        thr <- cmapdat$jaccard_threshold
        if (is.null(thr)) thr <- 0.3
        if (cmapdat$distmetric == "Euclidean")
          plotjaccard(cmapdat$eta_ijk_e, threshold = thr)
        else
          plotjaccard(cmapdat$eta_ijk_h, threshold = thr)
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
      settingmenu <- menu(c("Choose the distance metric",
                            #"Choose the clustering method",
                            #"Choose the misplacement index options",
                            "Choose the number of clusters",
                            "Set cluster names",
                            "Choose color scheme",
                            "Edit pile label canonical names",
                            bold(magenta("Main menu"))), title=bold("  Settings"))
      menuLevel <- 1
      if(settingmenu == 6) { menuLevel <- 0 }
      if(settingmenu == 0) { menuLevel <- -1 }
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
          xy      <- cbind(cmapdat$x, cmapdat$y)
          kmax    <- max(2, floor(length(cmapdat$cardNames[, 1]) / 3))
          kvals   <- 2:kmax
          fit_hc  <- hclust(as.dist(distM), method = cmapdat$clustMethod)
          wss_vals <- sapply(kvals, function(k) {
            grps <- cutree(fit_hc, k = k)
            sum(sapply(unique(grps), function(g) {
              pts <- xy[grps == g, , drop = FALSE]
              if (nrow(pts) < 2) return(0)
              sum(sweep(pts, 2, colMeans(pts), "-")^2)
            }))
          })
          difflogs <- -diff(log(wss_vals))
          plot(kvals, wss_vals, pch=19, col=4, cex=0.6,
               ylim=c(0, 1.05*max(wss_vals)),
               xlab="No. of clusters",
               ylab="Within cluster sums of squares",
               main="", cex.main=1, axes=FALSE)
          grid(); axis(1); axis(2)
          abline(h=0, lwd=2)
          if (all(difflogs >= 0.1)) {
            cmapdat$nclust <<- kmax
          } else {
            cmapdat$nclust <<- kvals[min(which(difflogs < 0.1))]
          }
          update_clusters()
          points(cmapdat$nclust, wss_vals[cmapdat$nclust - 1],
                 cex=1.5, pch=18, col="navyblue")
        }
        if(clustNo == 2) {
          kmax     <- max(2, floor(length(cmapdat$cardNames[, 1]) / 4))
          kvals    <- 2:kmax
          fit_hc   <- hclust(as.dist(distM), method = cmapdat$clustMethod)
          sil_vals <- sapply(kvals, function(k) {
            grps <- cutree(fit_hc, k = k)
            mean(cluster::silhouette(grps, as.dist(distM))[, 3])
          })
          plot(kvals, sil_vals, pch=19, col=4, cex=0.6, type="b",
               ylim=c(0, 1.05*max(sil_vals)),
               xlab="No. of clusters", ylab="Average silhouette width",
               main="", axes=FALSE)
          grid(); axis(1); axis(2)
          cmapdat$nclust <<- kvals[which.max(sil_vals)]
          update_clusters()
          points(cmapdat$nclust, sil_vals[which.max(sil_vals)],
                 cex=1.5, pch=18, col="navyblue")
        }
        if(clustNo == 3) {
          fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
          #cmapdat$clusMember <- cutree(fit.Clust, k=cmapdat$nclust)
          showDendrogram()
          curNC <- cmapdat$nclust
          tmpnclust <- rcmenu(1:floor(length(cmapdat$cardNames[,1])/4),
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
        cat(clusterNum, "[", clusterNames()[clusterNum], "]\n")
        groups <- clusterConfig()
        cat(bold("Statements in the cluster"), "\n")
        cardsInCluster <- which(groups == clusterNum)
        cat(blue(cmapdat$cardNames[cardsInCluster, "Statement"]), sep="\n")
        tmplbl <- cmapdat$pileLabels[which(cmapdat$pileLabels$cardNum %in% cardsInCluster),]
        if(nrow(tmplbl) > 0) {
          tmplbl <- sort(table(tmplbl[,1]), decreasing = TRUE)
          toShow <- min(5, nrow(tmplbl))
          cat(bold("Suggested names"), "\n")
          cat(green(names(tmplbl[1:toShow])), sep="\n")
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
      if(settingmenu == 5) {  ## edit pile label canonical names
        topLine()
        ld <- cmapdat$label_dict
        if (is.null(ld) || nrow(ld) == 0) {
          cat(yellow("No pile label data available.\n"))
          readline("Press Enter to continue. ")
        } else {
          repeat {
            topLine()
            ld <- cmapdat$label_dict
            # Build one display row per unique canonical label, listing all
            # original lowercase labels that map to it.
            unique_can <- unique(ld$canonical)
            display_items <- vapply(unique_can, function(can) {
              orig <- unique(ld$lcpilelabels[ld$canonical == can])
              if (length(orig) == 1 && orig == can)
                can
              else
                yellow(sprintf("%s  [originally: %s]",
                               can, paste(orig, collapse=", ")))
            }, character(1))
            cat(bold("  Pile label canonical names\n"))
            cat(yellow("(highlighted entries were standardized by fuzzy matching)\n\n"))
            choice <- menu(c(display_items, bold(magenta("Settings menu"))),
                           title = "Select a label to rename (or choose Settings menu to return)")
            if (choice == 0 || choice == length(unique_can) + 1) break
            old_can <- unique_can[choice]
            cat(sprintf("\nCurrent canonical label: '%s'\n", old_can))
            new_can <- trimws(readline("New canonical label (Enter to keep current): "))
            if (nzchar(new_can) && new_can != old_can) {
              cmapdat$label_dict$canonical[
                cmapdat$label_dict$canonical == old_can] <<- new_can
              cmapdat$pileLabels$pileLabel[
                cmapdat$pileLabels$pileLabel == old_can] <<- new_can
              # Persist the updated dictionaries and session
              card_name_dict <- cmapdat$cardNames
              sorters_dict   <- cmapdat$sorters_dict
              label_dict     <- cmapdat$label_dict
              save(card_name_dict, sorters_dict, label_dict,
                   file = paste0(cmapdat$dataDir, "/dictionaries.RData"))
              save(cmapdat, file = paste0(cmapdat$dataDir, "/CMapSession.RData"))
              cat(green(sprintf("  Updated: '%s' -> '%s'\n", old_can, new_can)))
              readline("Press Enter to continue. ")
            }
          }
        }
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


#' Interactively select a project folder and load its data.
#'
#' Prompts the user to choose a directory via a GUI dialog, then calls
#' \code{initCMap()} to load all project files and \code{update_clusters()} to
#' compute the initial cluster assignment. The resulting \code{cmapdat} object
#' is stored in the calling environment via \code{<<-}.
#' @param enc File encoding passed through to \code{initCMap} (default \code{"UTF-8"}).
#' @param sep Column separator passed through to \code{initCMap} (default \code{","}).
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
#' @param caption Title shown in the folder-selection dialog (default
#'   \code{'Select data directory'}).
#' @return The user's selected directory.
#' @export
choose_directory = function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

#' Open a GUI folder-selection dialog (OS-specific implementation).
#'
#' A lower-level alternative to \code{choose_directory}. Uses
#' \code{utils::choose.dir()} on Windows, \code{zenity} on Linux, and an
#' AppleScript call on macOS.
#' @return A character string with the selected directory path, or an empty
#'   string if the user cancels.
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

#' Detect the current operating system.
#'
#' @return A character string: \code{"windows"}, \code{"Linux"}, or
#'   \code{"MacOSX"}.
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


#' Return the current cluster assignment vector.
#'
#' Runs hierarchical clustering on the active distance matrix using the current
#' settings in \code{cmapdat} and cuts the tree at \code{cmapdat$nclust} clusters.
#' @return A named integer vector of length S (number of statements) giving the
#'   cluster index (1 to \code{nclust}) for each statement.
clusterConfig <- function() {
  distM <- cmapdat$D2e
  if(cmapdat$distmetric == "Hyperbolic")
    distM <- cmapdat$D2h
  fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
  cutree(fit.Clust, k=cmapdat$nclust)
}

#' Return the display names for the current clusters.
#'
#' Looks up user-assigned names from \code{cmapdat$clustname} (keyed by cluster
#' count and index). Falls back to the cluster index number for any unnamed cluster.
#' @return A character vector of length \code{cmapdat$nclust} with the display
#'   name for each cluster.
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


#' Display an interactive selection menu with an optional pre-selected default.
#'
#' A replacement for \code{utils::menu()} that supports showing a pre-selected
#' (current) value in the prompt. Pressing Enter without a selection returns 0.
#' Only usable in interactive sessions.
#' @param choices Character vector of menu item labels.
#' @param graphics Ignored (kept for interface compatibility with \code{menu}).
#' @param title Optional title string printed above the menu.
#' @param slctd Integer index of the currently active selection, shown in the
#'   prompt as \code{[slctd]}. \code{NULL} (default) shows no default.
#' @return An integer: the index of the chosen item (1 to \code{length(choices)}),
#'   or 0 if the user presses Enter without selecting.
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
  "cardfilename" = "Make sure the file is spelled correctly (Statements.csv, case-sensitive), and that the data folder contains all the input files.",
  "cardfilesep" = "Make sure you use the correct column separator. If your statements contain commas, use tab as a column separator, and invoke the program with RCMapMenu(sep='\\t')",
  "pilesfilename" = "Make sure the file is spelled correctly (SortedCards.csv, case-sensitive), and that the data folder contains all the input files.",
  "pilesfilesep" = "Make sure you use the correct column separator. The first column in SortedCards.csv must contain the sorter IDs, the second column has to contain pile labels (text), and columns 3,4,... contain the cards in each pile.",
  "ratingsfilename" = "Make sure the file is spelled correctly (Ratings.csv, case-sensitive), and that the data folder contains all the input files.",
  "ratingsfilesep" = "The first two columns must be named RaterID	and StatementID (case-sensitive) and StatementID must match the values in the Statements.csv file",
  "demographicsfilename" = "Make sure the file is spelled correctly (Demographics.csv, case-sensitive), and that the data folder contains all the input files.",
  "demographicsfilesep" = "The first column must be named RaterID (case-sensitive) and must match the values in the Ratings.csv file",
  "demographicsnodata" = "Either no columns in the file, or all columns contain categorical data with distinct values, or all values in each column are identical.",
  "non-numeric" = "has a non-numeric value where a card number is expected"
)
