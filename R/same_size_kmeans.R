# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' get K centers with kmeans++ algorithm
#' 
#' @param x a matrix like input data
#' @param k total number of centers
#' 
#' @return the output is a 0-1 vector, each 1 indicate a center 

KmeansppInital <-function (x, k) {
  data.num = nrow(x)
  centers = sample(1:data.num,1)
  while (length(centers) < k) {
    dist = apply(x, 1, function(i) {
      #distance to the nearest center
      min(sapply(centers, function(cid) {
       GetSingleDist(i, x[cid,]) ^ 2
      }))
    })
    #cat(dist/sum(dist))
    #cat("\n")
    if (sum(dist) == 0){dist = rep(1, data.num)}
    centers = c(centers, sample(1:data.num, 1, prob = dist/sum(dist)))
  }
  sort(centers)
}

CalMeans <- function(x, clusters){
  clsts = unique(clusters)
  res = matrix(0, length(clsts),ncol(x))
  for(i in clsts){
    res[i,] = colMeans(x[clusters == i,])
  }
  res
}

GetSingleDist <-function(x, y) {
  #sqrt(sum((x - y) ^ 2))
  1 - (x %*% y / sqrt(x%*%x * y%*%y))
}

GetMultiDist  <- function (x, centers) {
  apply(centers, 1,function(i) {GetSingleDist(x, i)}) 
}

InitalAssign <- function(x, centers.idx, cluster.size) {
  clst.num = length(centers.idx)
  clst.lst = rep(1, clst.num)
  data.num = nrow(x)
  clst.asg = rep(0, data.num)
  
  for (i in 1:clst.num){
    clst.asg[centers.idx[i]] = i
  }
  
  centers = x[centers.idx,]
  
  max.size = ceiling(data.num / clst.num)
  
  dist.mat = t(apply(x, 1, GetMultiDist, centers = centers))
  for (idx in 1:data.num) {
    if (idx %in% centers.idx) {
      next
    }
    ord.lst = order(dist.mat[idx,])
    for (ord in 1:clst.num){
      clst.idx = ord.lst[ord]
      if (clst.lst[clst.idx] < max.size) {
        clst.asg[idx] = clst.idx
        clst.lst[clst.idx] = clst.lst[clst.idx] + 1
        break
      }
    }
  }
  clst.asg
}


GetPrioLst <- function(x, clusters, means){
  ids = 1:nrow(x)
  cur = clusters
  dist = t(apply(x, 1, GetMultiDist, centers = means))
  cur.dist = sapply(ids, function(i){ dist[i,clusters[i]]})
  delta = sapply(ids, function(i){ dist[i,clusters[i]] - min(dist[i,])})
  #delta current assign distance minus best alternative, which is improvement
  alt = sapply(ids, function(i){ which.min(dist[i,])})
  res = data.frame(ids, cur, delta, alt, cur.dist)
  res[order(res$delta, decreasing = TRUE),]
}
#' get K clusters with same size kmeans algorithm
#' 
#' @param x a matrix like input data
#' @param k total number of clusters
#' 
#' @return the output is a vector, each element indicate a row's cluster number 
SameSizeKmeans <- function(x, k, max.iter = 100) {
  data.num = nrow(x)
  cluster.size = ceiling(data.num / k)
  centers.idx = KmeansppInital(x, k)
  clusters = InitalAssign(x, centers.idx, k)
  size.lst = sapply(1:k, function (i) {length(clusters[clusters == i])})
  
  iter = 0
  while (iter < max.iter) {
    active = 0
    tranlst = rep(-1, k)
    means = CalMeans(x, clusters)
    prio.lst = GetPrioLst(x, clusters, means)
    for (pid in 1:nrow(prio.lst)) {
      transferred = FALSE
      if (prio.lst[pid,]$delta <= 0) {
	      next
      }
      for (other.clst in 1:k){
        if(other.clst == prio.lst[pid,]$cur) {
		next
	}
        gain = prio.lst[pid,]$cur.dist - GetSingleDist(x[prio.lst[pid,]$ids,], means[other.clst,])
        if(size.lst[other.clst] < cluster.size & gain > 0){
          clusters[prio.lst[pid,]$ids] = other.clst
          size.lst[other.clst] = size.lst[other.clst] + 1
          active = active + 1
          transferred = TRUE
          break
        }
        if(tranlst[other.clst] != -1){#target cluster has point  waiting for transfer
          trans.id = tranlst[other.clst]
          tgain = GetSingleDist(x[trans.id,], means[other.clst,]) - GetSingleDist(x[trans.id, ], means[prio.lst[pid,]$cur,])
          if(tgain + gain > 0){
            #exchange two point
            active = active + 2
            transferred = TRUE
            clusters[prio.lst[pid,]$ids] = other.clst
            clusters[trans.id] = prio.lst[pid,]$cur
            tranlst[other.clst] = -1
            break
          }
        }
      }
      if (transferred == FALSE) {
        tranlst[prio.lst[pid,]$cur] = prio.lst[pid,]$ids
        }
    }
    if(active == 0){
      break
    }
    iter = iter + 1
    cat(paste("number of transfers in iter", iter, active, "\n"))
    #cat(paste(clusters))
  }
  
  cat(paste("K-means clustering converged in", iter, "iters", "\n"))
  clusters
}

