
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
getAdjMatrices <- function(piledat, cardNames, showWarnings=TRUE) {
  sorters <- sort(as.numeric(unique(piledat[,1])))
  cardsused <- na.exclude(unique(stack(piledat[,3:ncol(piledat)])[,1]))
  cardsused <- sort(as.numeric(setdiff(cardsused, "")))
  nCards <- length(cardsused)
  issues <- ""
  mat.list <- list()
  # Making a matrix for each individual and storing them in a list
  for (i in 1:length(sorters)) {
    sorter.dat <- piledat[which(piledat[,1] == sorters[i]),]
    adj.mat <- matrix(0, nrow=nCards, ncol=nCards)
    for (j in 1:nrow(sorter.dat)) {
      cards <- sorter.dat[j, which(!is.na(sorter.dat[j,]))]
      cards <- as.numeric(sort(unlist(cards[-(1:2)])))
      vct   <- rep(0, nCards)
      vct[which(cardNames[,1] %in% cards)] <- 1
      adj.mat <- adj.mat + vct%*%t(vct)
    }
    rsum <- rowSums(adj.mat)
    if(any(rsum == 0) & showWarnings) {
      didnotsort <- paste(cardNames[which(rsum == 0),1], collapse=",")
      issues <- issues %+% yellow(sprintf("Sorter  %d did not sort card(s) %s\n",
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
    if(length(gt1) > 0) { cat(sorters[i], "\n", gt1,"\n\n", issues,"\n") }
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
  dd <- matrix(0, P,P)
  for (i in 1:P) {
    for (j in 1:P) {
      dd[i,j] <- dd[j,i] <- acosh(1+2*(sumsq(x[i,]-x[j,]))/
                                    ((1-sqnorms[i])*(1-sqnorms[j])))
    }
  }
  return(dd)
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
  M0 <- (abs(D-d))^(1+cst*pmin(d,D))
  M <- matrix(0,nrow=nrow(M0), ncol=ncol(M0)-1)
  for (i in 1:nrow(M)) { M[i,] <- M0[i,-i] }
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
  cors <- rep(0,B)
  for (i in 1:B) {
    set.seed(seed+i)
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
    axis(1);axis(2); abline(h=0, col="grey", lwd=2); grid()
  }
  list(cors=cors, distmetric=distmetric)
}

#' Read the card-sorting data files and generate the MDS data.
#'
#' Takes two input file names.
#' @param dataDir The project's directory.
#' @return A list with 9 elements n.indiv=number of sorters, cardNames=the statements on the cards, adj.mat=the adjacency matrix, DS=the pairwise distances matrix in the original high-dimensional space, D2e=the Euclidean distance matrix in the 2-D space, D2h=the hyperbolic distance matrix in the 2-D space, x=the x coordinates in the MDS plot, y=the y coordinates, stress=the stress of the MDS plot.
#' @importFrom smacof mds
#' @export
initCMap <- function(dataDir) {
  cardNamesFile <- paste0(dataDir, "/Statements.csv")
  cardDatFile <- paste0(dataDir, "/SortedCards.csv")
  raitingsFile <- paste0(dataDir, "/Ratings.csv")
  demographicsFile <- paste0(dataDir, "/Demographics.csv")
  # if any is missing, show an error message, hit any key to go back, return

  if(!file.exists(cardNamesFile))
    stop("Card name File:", cardNamesFile, "Doesn't exist" )
  if(!file.exists(cardDatFile))
    stop("Card sortsing File:", cardDatFile, "Doesn't exist" )
  cardNames <- read.csv(cardNamesFile, header=TRUE)
  cardDat <- suppressWarnings(read.csv(cardDatFile, header=FALSE,
                                       encoding = "UTF-8",
                                       colClasses = rep("character", 2+nrow(cardNames))))
  pileLabels <- data.frame(pileLabel=character(), cardNum=integer())
  for (i in 1:nrow(cardDat)) {
    cardDat[i,] <- gsub("\\s+$", "", cardDat[i,])
    pile <- cardDat[i,]
    pileLabel <- cardDat[i, 2]
    cardsInPile <- na.omit(as.numeric(cardDat[i, -c(1,2)]))
    if (length(cardsInPile) == 0)
      next
    cardNum <- as.numeric(cardsInPile)
    pileLabel <- rep(pileLabel, length(cardsInPile))
    pileLabels <- rbind(pileLabels, data.frame(pileLabel, cardNum=cardNum))
  }
  cardDat[,1] <- gsub("\\D","",cardDat[,1])
  tmp <- getAdjMatrices(cardDat, cardNames)
  adj.mat <- tmp[[1]]
  issues <- tmp[[2]]
  n.indiv <- length(adj.mat)
  # multidimensional scaling, and clustering:
  DS <- distanceMatrix(adj.mat)
  rownames(DS) <- cardNames[,1]
  fit.MDS <- mds(DS)
  x <- fit.MDS$conf[,1]
  y <- fit.MDS$conf[,2]
  stress <- fit.MDS$stress
  D2e <- as.matrix(dist(fit.MDS$conf))
  D2h <- diskDist(D2e)
  rownames(D2e) <- rownames(D2h) <- cardNames[,1]
  cardNames <- cardNames[,2]

  if(!file.exists(raitingsFile))
    stop("Ratings File:", raitingsFile, "Doesn't exist" )
  ratings <- read.csv(raitingsFile, header=TRUE)
  ratings[,2] <- as.factor(ratings[,2])

  if(!file.exists(demographicsFile))
    stop("Demographics File:", demographicsFile, "Doesn't exist" )
  demographics.dat <- read.csv(demographicsFile, stringsAsFactors = TRUE)
  dframe <- data.frame(matrix(0,ncol=ncol(demographics.dat)-1, nrow=nrow(ratings)))
  cohorts <- data.frame("allRaters"=matrix(1, nrow=nrow(ratings), ncol=1))
  for(colnum in 2:ncol(demographics.dat)) {
    if (class(demographics.dat[,colnum]) == "factor") {
      dframe[,colnum-1] <- factor(rep("", nrow(ratings)),
                                  levels=levels(demographics.dat[,colnum]))
    } else {
      dframe[,colnum-1] <- rep(0, nrow(ratings))
    }
    for (i in 1:nrow(ratings)) {
      dframe[i,colnum-1] <- demographics.dat[,colnum][which(demographics.dat$RaterID ==
                                                              ratings$RaterID[i])]
    }
  }
  colnames(dframe) <- colnames(demographics.dat)[-1]
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
  if (!dir.exists(paste0(dataDir,"/output")))
    dir.create(paste0(dataDir,"/output"))
  list(n.indiv=n.indiv, cardNames=cardNames, adj.mat=adj.mat, DS=DS,
       D2e=D2e, D2h=D2h, x=x, y=y, stress=stress, ratings=ratings,
       clustname=list(), pileLabels=pileLabels, issues=issues,
       demographics=demographics.dat, ratingDemographics=dframe,
       clustMethod="ward.D2", distmetric="Hyperbolic",cohortCol=1,
       MPfunc="median", nclust=3, splhalf=list(), cohorts=cohorts,
       pcoordChoice=c(1,2), clusMember=rep(1,length(cardNames)),
       clrscm="rcmap", dataDir=dataDir)
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
    cnames[(colnum-1)*5 + 1:5] <- paste(colnames(ratingsDat)[colnum+2], c("N","Mean","SD","Min", "Max"), sep="_")
  }
  colnames(summ) <- cnames
  summ
}


#' Plot the 2-D representation of the pile-sorting data
#'
#' @param cols The colors of the data (default= all blue.)
#' @param metric The metric to be used (Euclidean or hyperbolic.)
#' @export
showMDSPlot <- function(cols="blue", metric=NULL) {
  if(!is.null(metric))
    if(metric == 0)
      return(NULL)
  if(is.null(metric)) {
    sz <- 0.7
    ttl <- ""
  } else{
    if (metric == "MPindex") {
      ttl <- "Misplacement index"
      distM <- cmapdat$D2e
      if(cmapdat$distmetric == "Hyperbolic")
        distM <- cmapdat$D2h
      sz <- MPindex(cmapdat$DS, distM, func=cmapdat$MPfunc)
      sz <- round(sz, digits = 1) + 0.3 # to be between 0.3 (for 0) and 1.3 (for 1)
      lgd <- seq(0.3, 1.3, by=0.2)
      lgdkey <- lgd-0.3
    } else {
      ttl <- paste(colnames(cmapdat$ratings)[metric+2],
                   colnames(cmapdat$cohorts)[cmapdat$cohortCol])
      ratingstmp <- cmapdat$ratings[which(cmapdat$cohorts[,cmapdat$cohortCol] == 1),]
      rtgsmm <- ratingSummary(ratingstmp)
      rownames(rtgsmm) <- c(rownames(cmapdat$DS), "total")
      sz <- rtgsmm[-nrow(rtgsmm),
                   which(colnames(rtgsmm) ==
                           paste0(colnames(cmapdat$ratings)[metric+2],"_Mean"))]
      sz <- round((sz-1)/4, digits = 1) + 0.3 # to be between 0.3 (for 1) and 1.3 (for 5)
      lgd <- seq(0.3, 1.1, by=0.2)
      lgdkey <- 1+(lgd-0.3)*5
    }
  }
  plot(cmapdat$x, cmapdat$y,pch=19,col=cols, cex=sz, main=ttl,xlab="",
       ylab="", xlim=c(min(cmapdat$x)-0.25,max(cmapdat$x)*1.2),
       ylim=c(min(cmapdat$y),max(cmapdat$y)*1.2),axes=F, cex.main=1)
  text( cmapdat$x, cmapdat$y, pos=3, cex=0.8,
        labels=names(cmapdat$x))
  if(length(sz) > 1) {
    rect(min(cmapdat$x)-0.3, 0.2*max(cmapdat$y),
         min(cmapdat$x)-0.1, 1.3*max(cmapdat$y), col="whitesmoke", border="white")
    points(rep(min(cmapdat$x)-0.15, length(lgd)),
           seq(0.3*max(cmapdat$y), 1.2*max(cmapdat$y),length=length(lgd)),
           col=cols, cex = lgd, pch=19)
    text(rep(min(cmapdat$x)-0.25, length(lgd)),
         seq(0.3*max(cmapdat$y), 1.2*max(cmapdat$y),length=length(lgd)),
         lgdkey,cex=0.6)
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
    }
    SEM <- SD/sqrt(SS)
    Yl = ceiling(max(MN+SEM+0.5, na.rm = TRUE))
    bp <- barplot(MN, main=ttl,col=cols,ylim=c(0,Yl),
                  names.arg=clusterNames(),cex.names=0.8,las=3,
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
  plot( MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2], ylim=c(1,5), xlim=c(1,5.9),
        main="", col=0, axes=F, xlab=cmapdat$pcoordChoice[1],
        ylab=cmapdat$pcoordChoice[2], pch=19, cex=0.6)
  text(MN[-nrow(rtgsmm),1], MN[-nrow(rtgsmm),2],
       rownames(cmapdat$DS), cex=0.7,
       col=cols[groups], font=2)
  axis(1, at = 1:5); axis(2)
  grid()
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
      noPiles <- length(unique(apply(adjM,1,paste0, collapse="")))
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
    if(distmetric == "Hyperbolic")
      distM <- D2h
    fit.Clust <- hclust(as.dist(distM), method=clustMethod)
    groups <- cutree(fit.Clust, k=nclust)
    retdf <- as.data.frame(matrix(0, ncol=3+5*(ncol(ratings)-2), nrow=nclust+nrow(DS)))
    k <- 0
    for (i in 1:nclust) {
      cardsInCluster <- names(which(groups == i))
      cNames <- c(cardNames[which(groups == i)],"")
      cardsInCluster <- which(ratings$StatementID %in% cardsInCluster)
      rtgsmm <- ratingSummary(ratings[cardsInCluster,])
      rtgsmm <- data.frame(cNames, c(names(cmapdat$x)[which(groups == i)],""),
                           rep(i,nrow(rtgsmm)), rtgsmm)
      colnames(rtgsmm)[1:3] <- c("Statement", "CardNo", "ClusterNo")
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
    strs <- c("ward.D", "ward.D2", "single","complete", "average","mcquitty",
              "median", "centroid")
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
#' @export
#' @examples
#' \donttest{
#' #RCMapMenu()
#' }
RCMapMenu <- function() {
  menuLevel <- 0
  while(menuLevel >= 0) {
    topLine()
    ## Main menu
    if(menuLevel == 0) {
      rcmapmenu <- menu(c(ifelse(exists("cmapdat"),"Choose the data folder",bold("Choose the data folder")),
                          ifelse(exists("cmapdat"),"Summary",silver("summary")),
                          ifelse(exists("cmapdat"),"Settings",silver("Settings")),
                          ifelse(exists("cmapdat"),"Plots",silver("Plots")),
                          ifelse(exists("cmapdat"),"Reports",silver("Reports")),
                          ifelse(exists("cmapdat"),"Analysis",silver("Analysis")),
                          bold(magenta("R prompt"))), title=bold("  Main menu"))
      if (rcmapmenu %in% c(0,7))
        menuLevel <- -1 ## get the R prompt back
      if (rcmapmenu == 1) {
        loadRCMapData()
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
                         bold(magenta("Main menu"))))
      if(summMenu == 1) {
        cmapdat$splhalf <<- splitHalf(cmapdat$adj.mat, B=20,
                                      disttype = cmapdat$distmetric,
                                      plotit=TRUE, seed=23456)
      }
      menuLevel <- 1
      if(summMenu %in% c(0,2))
        menuLevel <- 0
    }
    ## settings
    if(rcmapmenu == 3) {
      topLine()
      showSettingMenu <- TRUE
      if(exists("settingmenu"))
        if(settingmenu == 5)
          showSettingMenu <- FALSE
      if(showSettingMenu)  {
        settingmenu <- menu(c("Choose the distance metric",
                              "Choose the clustering method",
                              "Choose the misplacement index options",
                              "Choose the number of clusters",
                              "Set cluster names",
                              "Choose color scheme",
                              bold(magenta("Main menu"))), title=bold("  Settings"))
        menuLevel <- ifelse(settingmenu %in% c(0,7), 0, 1)
      }
      if(settingmenu == 1) {
        topLine()
        dtype <- rcmenu(toText("distance"),
                        title=bold("  Distance function"),
                        slctd=which(toText("distance") == cmapdat$distmetric))
        if(dtype > 0)
          cmapdat$distmetric <<- toText("distance", dtype)
      }
      if(settingmenu == 2) {
        topLine()
        cmeth <- rcmenu(toText("clusteringMethod"),
                        title=bold("  Clustering method"),
                        slctd=which(toText("clusteringMethod") == cmapdat$clustMethod))
        if(cmeth > 0)
          cmapdat$clustMethod <<- toText("clusteringMethod", cmeth)
      }
      if(settingmenu == 3) { ## misplacement index
        topLine()
        mpfunc <- rcmenu(toText("MPparam"),
                         title=bold("  Misplacement index parameters"),
                         slctd=which(toText("MPparam") == cmapdat$MPfunc))
        if(mpfunc > 0)
          cmapdat$MPfunc <<- toText("MPparam", mpfunc)
      }
      if(settingmenu == 4) { ## the number of clusters
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
        }
        if(clustNo == 3) {
          fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
          #cmapdat$clusMember <- cutree(fit.Clust, k=cmapdat$nclust)
          showDendrogram()
          curNC <- cmapdat$nclust
          tmpnclust <- rcmenu(1:floor(length(cmapdat$cardNames)/3),
                              title = bold("Set the number of clusters"), slctd = curNC)
          if(tmpnclust > 0)
            cmapdat$nclust <<- tmpnclust
          #cmapdat$clusMember <- cutree(fit.Clust, k=cmapdat$nclust)
          showDendrogram()
        }
      }
      if(settingmenu == 5) { ## the cluster names
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
        cat("\nEnter a name [",clNames[clusterNum],"]\n")
        clName <- scan("",quiet = TRUE, what = "chr", n=1, sep = "\n")
        if (length(clName) == 0)
          clName <- clNames[clusterNum]
        key <- sprintf("clust_%02d_%02d", length(clNames), clusterNum)
        cmapdat$clustname[[key]] <<- clName
      }
      if(settingmenu == 6) {
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
      menuLevel <- ifelse(plotmenu %in% c(0,13), 0, 1)
      if(plotmenu %in% 7:10) {
        cmapdat$cohortCol <<- menu(colnames(cmapdat$cohorts))
        if (cmapdat$cohortCol == 0)
          next
      }
      if(plotmenu %in% c("11","12")) {
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
      menuLevel <- ifelse(reportsmenu %in% c(0,4), 0, 1)
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
      menuLevel <- ifelse(analysismenu%in% c(0,3), 0, 1)
      if(analysismenu < 3)
        switch (analysismenu,
                "1" = clusterANOVA(),
                "2" = clusterTUKEY(),
        )
    }
  }
}


#' Load the input files for a concept mapping project.
loadRCMapData <- function() {
  topLine()
  dataDir <- paste0(choose_directory(),"/")
  cat(bold(crayon::green("Opening project folder:" %+% dataDir %+% "...\n\n")))
  cmapdat <<- initCMap(dataDir)
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
  fit.Clust <- hclust(as.dist(distM), method=cmapdat$clustMethod)
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

