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

kmeanspp <-function (x, k) {
  data_num = nrow(x)
  centers = sample(1:data_num,1)
  while (length(centers) < k) {
    dist = apply(x, 1, function(i) {
      #distance to the nearest center
      min(sapply(centers, function(cid) {
       sum((i - x[cid,]) ^ 2) 
      }))
    })
    #cat(dist/sum(dist))
    #cat("\n")
    centers = c(centers, sample(1:data_num, 1, prob = dist/sum(dist)))
  }
  sort(centers)
}

cal_means <- function(x, clusters){
  clsts = unique(clusters)
  res = matrix(0, length(clsts),2)
  for(i in clsts){
    res[i,] = colMeans(x[clusters == i,])
  }
  res
}

dist_s <-function(x, y) {
  #sqrt(sum((x - y) ^ 2))
  1 - (x %*% y / sqrt(x%*%x * y%*%y))
}
dist_f  <- function (x, centers) {
  apply(centers, 1,function(i) {dist_s(x, i)}) 
}

inital_assign <- function(x, centers_idx, cluster_size) {
  clst_num = length(centers_idx)
  clst_lst = rep(1, clst_num)
  data_num = nrow(x)
  clst_asg = rep(0, data_num)
  
  for (i in 1:clst_num){
    clst_asg[centers_idx[i]] = i
  }
  
  centers = x[centers_idx,]
  
  max_size = ceiling(data_num / clst_num)
  
  dist_mat = t(apply(x, 1, dist_f, centers = centers))
  for (idx in 1:data_num) {
    if (idx %in% centers_idx) {
      next
    }
    ord_lst = order(dist_mat[idx,])
    for (ord in 1:clst_num){
      clst_idx = ord_lst[ord]
      if (clst_lst[clst_idx] < max_size) {
        clst_asg[idx] = clst_idx
        clst_lst[clst_idx] = clst_lst[clst_idx] + 1
        break
      }
    }
  }
  clst_asg
}


get_prio_lst <- function(x, clusters, means){
  ids = 1:nrow(x)
  cur = clusters
  dist = t(apply(x, 1, dist_f, centers = means))
  cur_dist = sapply(ids, function(i){ dist[i,clusters[i]]})
  delta = sapply(ids, function(i){ dist[i,clusters[i]] - min(dist[i,])})
  #delta current assign distance minus best alternative, which is improvement
  alt = sapply(ids, function(i){ which.min(dist[i,])})
  res = data.frame(ids, cur, delta, alt, cur_dist)
  res[order(res$delta, decreasing = TRUE),]
}
#' get K clusters with same size kmeans algorithm
#' 
#' @param x a matrix like input data
#' @param k total number of clusters
#' 
#' @return the output is a vector, each element indicate a row's cluster number 
SameSizeKmeans <- function(x, k, max_iter = 100) {
  data_num = nrow(x)
  cluster_size = ceiling(data_num / k)
  centers_idx = kmeanspp(x, k)
  clusters = inital_assign(x, centers_idx, k)
  size_lst = sapply(1:k, function (i) {length(clusters[clusters == i])})
  
  iter = 0
  while (iter < max_iter) {
    active = 0
    tranlst = rep(-1, k)
    means = cal_means(x, clusters)
    prio_lst = get_prio_lst(x, clusters, means)
    for (pid in 1:nrow(prio_lst)) {
      transferred = FALSE
      if (prio_lst[pid,]$delta <= 0) {next}
      for (other_clst in 1:k){
        if(other_clst == prio_lst[pid,]$cur) {next}
        gain = prio_lst[pid,]$cur_dist - dist_s(x[prio_lst[pid,]$ids,], means[other_clst,])
        if(size_lst[other_clst] < cluster_size & gain > 0){
          clusters[prio_lst[pid,]$ids] = other_clst
          size_lst[other_clst] = size_lst[other_clst] + 1
          active = active + 1
          transferred = TRUE
          break
        }
        if(tranlst[other_clst] != -1){#target cluster has point  waiting for transfer
          trans_id = tranlst[other_clst]
          tgain = dist_s(x[trans_id,], means[other_clst,]) - dist_s(x[trans_id, ], means[prio_lst[pid,]$cur,])
          if(tgain + gain > 0){
            #exchange two point
            active = active + 2
            transferred = TRUE
            clusters[prio_lst[pid,]$ids] = other_clst
            clusters[trans_id] = prio_lst[pid,]$cur
            tranlst[other_clst] = -1
            break
          }
        }
      }
      if (transferred == FALSE) {
        tranlst[prio_lst[pid,]$cur] = prio_lst[pid,]$ids
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

