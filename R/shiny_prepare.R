#' Find DE from HC with recommended k
#' @import genefilter
#' @keywords internal


FindHcDe <- function(cell.hc.clust, U.score.non.pca, cell.name){
  # prepare the process bar
  cat("\nCalculate p value of non-overlap clusters...\n")
  k.hc.clu <- dim(cell.hc.clust)[2]
  dot.num.non <- 2
  process.bar <- c(".")
  bar.count <- 1
  gene.name <- colnames(U.score.non.pca)

  while(bar.count < dot.num.non){
    process.bar <- paste(process.bar, ".", sep = ".")
    bar.count <- bar.count + 1
  }

  delete.ind.clu <- vector()

  de.p.non.overlap <- lapply(1:k.hc.clu, function(x){
    labels <- as.factor(cell.hc.clust[,x])
    ps <- colFtests(as.matrix(U.score.non.pca), labels)$p.value
    ps.ind <- which(ps < 0.05)
    ps.ind.order <- ps.ind[order(ps[ps.ind])]
    gene.name.choose <- gene.name[ps.ind.order]

    cat(process.bar)

    if(length(gene.name.choose) == 0) {
      delete.ind.clu <- c(delete.ind.clu, x)
      next
    }

    if(length(gene.name.choose) > 100) gene.name.choose <- gene.name.choose[1:100]

    hc.clu.x <- cell.hc.clust[,x]
    reorder.cell.ind <- order(hc.clu.x)

    data.de <- U.score.non.pca[reorder.cell.ind, gene.name.choose]
    return(data.de)
  })

  cat("done\n")

  return(list(hc.de.plot = de.p.non.overlap))
}

#' Find DE from Overlap clusters
#' @import genefilter
#' @keywords internal

FindOverlapDe <- function(clu.merge, U.score.non.pca, cell.name, n){
  k.group <- length(clu.merge)
  overlap.de.plot <- list()
  delete.clu <- vector()
  n = n
  gene.name <- colnames(U.score.non.pca)
  for(i in 1:k.group){
    cat(paste("\nCalculate p value of overlap clusters for the no.",
              i, " group clustering result...\n", sep = ""))

    pure.clu <- clu.merge[[i]][,-1]
    k.sub <- dim(pure.clu)[2]

    # Set the process bar
    dot.num <- 2
    process.bar <- c(".")
    bar.count <- 1
    while(bar.count < dot.num){
      process.bar <- paste(process.bar, ".", sep = ".")
      bar.count <- bar.count + 1
    }

    overlap.de.data <- list()
    for(j in 1:k.sub){
      y <- pure.clu[,j]
      cell.order <- order(y, decreasing = TRUE)

      labels <- as.factor(y)
      ps <- colttests(as.matrix(U.score.non.pca), labels)$p.value
      ps <- p.adjust(ps)

      ps.ind <- which(ps < 0.05)
      ps.ind.order <- ps.ind[order(ps[ps.ind])]
      gene.name.choose <- gene.name[ps.ind.order]

      if(length(gene.name.choose) < 10){
        delete.clu <- rbind(c(i,j), delete.clu)
        next
      }

      if(length(gene.name.choose) > 100) gene.name.choose <- gene.name.choose[1:100]

      data.de <- U.score.non.pca[cell.order, gene.name.choose]
      data.de <- as.data.frame(data.de)

      cat(process.bar)

      overlap.de.data[[j]] <- data.de
    }

    overlap.de.plot[[i]] <- overlap.de.data
  }

  if(length(delete.clu) > 0){
    whole.del <- vector()
    del.big <- unique(delete.clu[,1])
    for(i in del.big){
      del.samll <- delete.clu[delete.clu[,1] == i,2] + 1
      k.sub <- dim(clu.merge[[i]])[2] - 1
      if(length(del.samll) == k.sub) {
        clu.merge.whole.del <- c(whole.del, i) # delete the whole cluster
      } else {
        for(j in del.samll)
        clu.merge[[i]] <- clu.merge[[i]][, -j]
        }
    }
    clu.merge <- clu.merge[-whole.del]
  }

  cat("done\n")
  return(list(overlap.de.plot = overlap.de.plot, clu.merge.1 = clu.merge))
}

#' Find DE from original overlap clusters
#' @import genefilter
#' @keywords internal

FindBefDe <- function(mer.clu, overlap.clu, U.score.non.pca, cell.name){
  mer.clu.sum <- length(mer.clu)
  mer.clu.new <- vector()
  gene.name <- colnames(U.score.non.pca)
  bef.de.plot <- list()
  cat(paste("\nCalculate p value of potential overlap part ...\n", sep = ""))

  dot.num <- 2
  process.bar <- c(".")
  bar.count <- 1
  while(bar.count < dot.num){
    process.bar <- paste(process.bar, ".", sep = ".")
    bar.count <- bar.count + 1
  }

  for(i in 1:mer.clu.sum){
    x <- mer.clu[i]
    x <- unlist(strsplit(x, ", ", fixed = TRUE))
    labels.left.ind <- which(overlap.clu[,x][, 1] - overlap.clu[,x][, 2] == -1)
    labels.right.ind <- which(overlap.clu[,x][, 1] == 1 & overlap.clu[,x][, 2] == 1)
    data.bef.merge <- U.score.non.pca[c(labels.left.ind, labels.right.ind),]

    left.len <- length(labels.left.ind)
    right.len <- length(labels.right.ind)
    if(left.len < 5 |right.len < 5) next

    cell.order <- c(labels.left.ind, labels.right.ind)
    my.labels <- as.factor(c(rep(1, left.len), rep(0, right.len)))

    ps <- colttests(as.matrix(data.bef.merge), my.labels)$p.value
    ps <- p.adjust(ps)
    ps.ind <- which(ps < 0.05)
    if(length(ps.ind) < 20){
      next
    }

    ps.ind.order <- ps.ind[order(ps[ps.ind])]
    gene.name.choose <- gene.name[ps.ind.order]


    if(length(gene.name.choose) > 100) gene.name.choose <- gene.name.choose[1:100]

    data.de <- U.score.non.pca[cell.order, gene.name.choose]
    data.de <- as.data.frame(data.de)
    mer.clu.new <- c(mer.clu.new, mer.clu[i])

    now.ind <- length(bef.de.plot) + 1
    bef.de.plot[[now.ind]] <- data.de

    cat(process.bar)

  }
  cat("done\n")
  return(list(bef.de.plot = bef.de.plot, mer.clu = mer.clu.new))
}

#' Prepare the data for visualization
#' @import reshape2
#' @keywords internal

OverlapMelt <- function(K, overlap.clu, data.tsne){
  overlap.res <- overlap.clu[[K]]
  cell.overlap.size <- apply(overlap.res[,-1], 1, sum)
  cell.overlap.size <- as.data.frame(cell.overlap.size)
  X1 = data.tsne[,1]
  X2 = data.tsne[,2]
  overlap.res <- cbind(overlap.res, cell.overlap.size, X1, X2)

  overlap.res.melt <- melt(overlap.res,
                           id.vars = c("cell", "cell.overlap.size", "X1", "X2"))
  overlap.res.melt$value = factor(overlap.res.melt$value)
  return(overlap.res.melt)
}


## Find DE gene by using DESeq
