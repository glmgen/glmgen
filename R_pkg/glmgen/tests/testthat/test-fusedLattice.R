library(genlasso)
library(glmgen)
library(testthat)

set.seed(0)
EPS = 1e-4

# (1) Test fusedLattice and fusedGraph over a grid; compare to
#   the genlasso::fusedlasso2d
set.seed(1)
y = matrix(rnorm(4*12), ncol=4, nrow=12)
y[3:4,1:4] = y[3:4,1:4] + 2

D = getD2d(nrow(y), ncol(y))
fp = fusedlasso2d(y)

mat = matrix(1:length(y), ncol=ncol(y), nrow=nrow(y))
el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fl1 = fusedLattice(y, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="prox")
  fl2 = fusedLattice(y, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="dp")
  fg1 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="prox")
  fg2 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="dp")

  spath = matrix(coef(fp, lambda=thisLam)$beta, nrow=nrow(y), ncol=ncol(y))
  fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y))
  fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y))
  fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y))
  fg2_beta = matrix(fg2$beta, nrow=nrow(y), ncol=ncol(y))

  expect_true(max(abs(spath - fl1_beta)) < EPS)
  expect_true(max(abs(spath - fl2_beta)) < EPS)
  expect_true(max(abs(spath - fg1_beta)) < EPS)
  expect_true(max(abs(spath - fg2_beta)) < EPS)
}

# (2) Test weights of fusedLattice and fusedGraph over a grid; compare to
#   the genlasso::fusedlasso2d
set.seed(1)
y = matrix(rnorm(4*12), ncol=4, nrow=12)
y[3:4,1:4] = y[3:4,1:4] + 2

D = getD2d(nrow(y), ncol(y))
fp = fusedlasso2d(y)

mat = matrix(1:length(y), ncol=ncol(y), nrow=nrow(y))
el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))
ew = lapply(el, function(v) rep(1,length(v)))
ew = lapply(ew, function(v) v*5)
ewL = rep(5, nrow(y)*(ncol(y)-1) + ncol(y)*(nrow(y)-1))

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fl1 = fusedLattice(y*5, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="prox")
  fl2 = fusedLattice(y*5, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="dp")
  fl3 = fusedLattice(y, lambda=thisLam/5, eps=1e-06,
                      maxIter=250, method="prox", edgeWeights=ewL)
  fg1 = fusedGraph(y*5, edges=el, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="prox")
  fg2 = fusedGraph(y*5, edges=el, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="dp")
  fg3 = fusedGraph(y, edges=el, lambda=thisLam/5, eps=1e-06,
                      maxIter=250, method="prox", edgeWeights=ew)

  spath = matrix(coef(fp, lambda=thisLam)$beta, nrow=nrow(y), ncol=ncol(y))
  fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fl3_beta = matrix(fl3$beta, nrow=nrow(y), ncol=ncol(y))
  fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fg2_beta = matrix(fg2$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fg3_beta = matrix(fg3$beta, nrow=nrow(y), ncol=ncol(y))

  expect_true(max(abs(spath - fl1_beta)) < EPS)
  expect_true(max(abs(spath - fl2_beta)) < EPS)
  expect_true(max(abs(spath - fl3_beta)) < EPS)
  expect_true(max(abs(spath - fg1_beta)) < EPS)
  expect_true(max(abs(spath - fg2_beta)) < EPS)
  expect_true(max(abs(spath - fg3_beta)) < EPS)
}

# (3) Test  fusedLattice and fusedGraph over a grid with missing values;
#  compare to the genlasso::fusedlasso (with properly defined graph)
set.seed(2)
y = matrix(rnorm(4*5),ncol=5, nrow=4)
y[1:3,3:4] = NA
indexNA = !is.na(y)
y2 = y[indexNA]

el = cbind(c(1,2,3,5,6,7,11,12,13,1,2,3,4,8,9,10),
           c(2,3,4,6,7,8,12,13,14,5,6,7,8,9,10,14))
gr = graph.edgelist(el, FALSE)
fp = fusedlasso(y2, gr=gr)

mat = matrix(1:length(y), ncol=ncol(y), nrow=nrow(y))
el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fl1 = fusedLattice(y, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="prox")
  fl2 = fusedLattice(y, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="dp")
  fg1 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="prox")
  fg2 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="dp")

  spath = matrix(NA, nrow=nrow(y), ncol=ncol(y))
  spath[indexNA] = coef(fp, lambda=thisLam)$beta
  fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y))
  fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y))
  fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y))
  fg2_beta = matrix(fg2$beta, nrow=nrow(y), ncol=ncol(y))

  expect_true(max(abs(spath - fl1_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fl2_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg1_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg2_beta),na.rm=TRUE) < EPS)
}

# (4) Test  fusedLattice and fusedGraph over a grid with missing values
#  AND weights; compare to the genlasso::fusedlasso (with a properly
#  defined graph)
set.seed(2)
y = matrix(rnorm(4*5),ncol=5, nrow=4)
y[1:3,3:4] = NA
indexNA = !is.na(y)
y2 = y[indexNA]

el = cbind(c(1,2,3,5,6,7,11,12,13,1,2,3,4,8,9,10),
           c(2,3,4,6,7,8,12,13,14,5,6,7,8,9,10,14))
gr = graph.edgelist(el, FALSE)
fp = fusedlasso(y2, gr=gr)

mat = matrix(1:length(y), ncol=ncol(y), nrow=nrow(y))
el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))
ew = lapply(el, function(v) rep(1,length(v)))
ew = lapply(ew, function(v) v*5)
ewL = rep(5, nrow(y)*(ncol(y)-1) + ncol(y)*(nrow(y)-1))

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fl1 = fusedLattice(y*5, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="prox")
  fl2 = fusedLattice(y*5, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="dp")
  fl3 = fusedLattice(y, lambda=thisLam/5, eps=1e-06,
                      maxIter=250, method="prox", edgeWeights=ewL)
  fg1 = fusedGraph(y*5, edges=el, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="prox")
  fg2 = fusedGraph(y*5, edges=el, lambda=thisLam*5, eps=1e-06,
                      maxIter=250, method="dp")
  fg3 = fusedGraph(y, edges=el, lambda=thisLam/5, eps=1e-06,
                      maxIter=250, method="prox", edgeWeights=ew)

  spath = matrix(NA, nrow=nrow(y), ncol=ncol(y))
  spath[indexNA] = coef(fp, lambda=thisLam)$beta
  fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fl3_beta = matrix(fl3$beta, nrow=nrow(y), ncol=ncol(y))
  fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fg2_beta = matrix(fg2$beta, nrow=nrow(y), ncol=ncol(y)) / 5
  fg3_beta = matrix(fg3$beta, nrow=nrow(y), ncol=ncol(y))

  expect_true(max(abs(spath - fl1_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fl2_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fl3_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg1_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg2_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg3_beta),na.rm=TRUE) < EPS)
}

# # (5) Test fusedLattice and fusedGraph at scale, comparing to
# #   each other (should give the same results)
# url = "http://photogrammar.research.yale.edu/photos/service/pnp/fsa"
# url = paste0(url, "/8a21000/8a21500/8a21552v.jpg")
# download.file(url, tf <- tempfile())
# y = jpeg::readJPEG(tf)

# mat = matrix(1:length(y), ncol=ncol(y), nrow=nrow(y))
# el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))

# for (thisLam in 10^(2:(-4))) {
#   fl1 = fusedLattice(y, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="prox")
#   fl2 = fusedLattice(y, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="dp")
#   fg1 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="prox")
#   fg2 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="dp")

#   fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y))
#   fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y))
#   fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y))
#   fg2_beta = matrix(fg2$beta, nrow=nrow(y), ncol=ncol(y))

#   expect_true(max(abs(fl1_beta - fl2_beta)) < EPS)
#   expect_true(max(abs(fl1_beta - fg1_beta)) < EPS)
#   expect_true(max(abs(fl1_beta - fg2_beta)) < EPS)
# }

# # Now with NA's
# y[y < 0.05] = NA
# for (thisLam in 10^(2:(-4))) {
#   fl1 = fusedLattice(y, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="prox")
#   fl2 = fusedLattice(y, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="dp")
#   fg1 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="prox")
#   fg2 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
#                       maxIter=250, method="dp")

#   fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y))
#   fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y))
#   fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y))
#   fg2_beta = matrix(fg2$beta, nrow=nrow(y), ncol=ncol(y))

#   expect_true(max(abs(fl1_beta - fl2_beta),na.rm=TRUE) < EPS)
#   expect_true(max(abs(fl1_beta - fg1_beta),na.rm=TRUE) < EPS)
#   expect_true(max(abs(fl1_beta - fg2_beta),na.rm=TRUE) < EPS)
# }

# (6) Test just one chain, comparing fusedGraph to genlasso::fusedlasso
set.seed(1)
n = 100
y = sin(seq(0,2*pi,length.out=n)) + rnorm(n,sd=0.25)

fp = fusedlasso1d(y)
el = list(1:length(y))

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fl1 = fusedLattice(matrix(y,nrow=1), lambda=thisLam, eps=1e-06,
                      maxIter=250, method="prox")
  fl2 = fusedLattice(matrix(y,nrow=1), lambda=thisLam, eps=1e-06,
                      maxIter=250, method="dp")
  fg1 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="prox")
  fg2 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-06,
                      maxIter=250, method="dp")

  spath = matrix(coef(fp, lambda=thisLam)$beta, nrow=1, ncol=length(y))
  fl1_beta = matrix(fl1$beta, nrow=1, ncol=length(y))
  fl2_beta = matrix(fl2$beta, nrow=1, ncol=length(y))
  fg1_beta = matrix(fg1$beta, nrow=1, ncol=length(y))
  fg2_beta = matrix(fg2$beta, nrow=1, ncol=length(y))

  expect_true(max(abs(spath - fl1_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fl2_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg1_beta),na.rm=TRUE) < EPS)
  expect_true(max(abs(spath - fg2_beta),na.rm=TRUE) < EPS)
}

# (7) Test a grid using an alternative edge set (two unbroken
#   chains that snake through the graph).
set.seed(1)
y = matrix(rnorm(4*12), ncol=4, nrow=12)
y[3:4,1:4] = y[3:4,1:4] + 2

D = getD2d(nrow(y), ncol(y))
fp = fusedlasso2d(y)

mat = matrix(1:length(y), ncol=ncol(y), nrow=nrow(y))
tmat = t(mat)
el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))

el = replicate(2, list(0))
el[[1]] = rep(NA, length(mat))
index = (0:(length(mat)-1) %/% nrow(mat)) %% 2
el[[1]][index == 0] = as.numeric(mat[,(0:(ncol(mat)-1) %% 2) == 0])
el[[1]][index == 1] = as.numeric(mat[nrow(mat):1,(0:(ncol(mat)-1) %% 2) == 1])
el[[2]] = rep(NA, length(tmat))
index = (0:(length(mat)-1) %/% nrow(tmat)) %% 2
el[[2]][index == 0] = as.numeric(tmat[,(0:(ncol(tmat)-1) %% 2) == 0])
el[[2]][index == 1] = as.numeric(tmat[nrow(tmat):1,(0:(ncol(tmat)-1) %% 2) == 1])

ew = lapply(el, function(v) rep(1, length(v)))
ew[[1]][(1:nrow(mat))[(0:(nrow(mat)-1)) %% 2  == 1]] = 0.5
ew[[1]][(1:nrow(mat))[(0:(nrow(mat)-1)) %% 2  == 0] + length(mat)-nrow(mat)] = 0.5
ew[[1]][which(abs(diff(el[[1]])) != 1)] = 0.5
ew[[2]][(1:ncol(mat))[(0:(ncol(mat)-1)) %% 2  == 1]] = 0.5
ew[[2]][(1:ncol(mat))[(0:(ncol(mat)-1)) %% 2  == 0] + length(mat)-ncol(mat)] = 0.5
ew[[2]][which(abs(diff(el[[2]])) == 1)] = 0.5

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fl1 = fusedLattice(y, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="prox")
  fg1 = fusedGraph(y, edges=el, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="prox", edgeWeights=ew)

  spath = matrix(coef(fp, lambda=thisLam)$beta, nrow=nrow(y), ncol=ncol(y))
  fl1_beta = matrix(fl1$beta, nrow=nrow(y), ncol=ncol(y))
  fg1_beta = matrix(fg1$beta, nrow=nrow(y), ncol=ncol(y))

  expect_true(max(abs(spath - fl1_beta)) < EPS*10)
  expect_true(max(abs(spath - fg1_beta)) < EPS*10)
  expect_true(max(abs(fl1_beta - fg1_beta)) < EPS*10)
}

# (8) Test just a single chain using an E matrix to
#     append an additional data point on at the end
#     which is constrainted to equal to the original
#     last data point; sample weights are adjusted
#     as well
set.seed(4)
n = 100
y = sin(seq(0,2*pi,length.out=n)) + rnorm(n,sd=0.25)
y2 = c(y,y[length(y)])

fp = fusedlasso1d(y)
el = list(1:length(y2))

w = rep(1, length(y2))
w[length(y)] = 0.5
w[length(y2)] = 0.5

E = Matrix::spMatrix(1, length(y2), c(1,1), c(n,n+1), c(1,-1))

lams = fp$lambda
for (thisLam in c(lams,0)) {
  fg1 = fusedGraph(y2, w, edges=el, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="prox", E=E)
  fg2 = fusedGraph(y2, w, edges=el, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="dp", E=E)

  spath = coef(fp, lambda=thisLam)$beta
  fg1_beta = fg1$beta[1:length(y)]
  fg2_beta = fg2$beta[1:length(y)]

  expect_true(max(abs(spath - fg1_beta)) < EPS)
  expect_true(max(abs(spath - fg2_beta)) < EPS)
}

# (10) Now, append an entire extra column onto a grid of
#      data equal to the original last column. Use E to
#      constrain the solution set, sample weights to preserve
#      the original problem, and set edge weights between the
#      rows of the final column to zero (these are already
#      fully penalized).
set.seed(4)
y = matrix(rnorm(4*12), ncol=4, nrow=12)
y[3:4,1:4] = y[3:4,1:4] + 2
y2 = cbind(y, y[,ncol(y)])

D = getD2d(nrow(y), ncol(y))
fp = fusedlasso2d(y)

mat = matrix(1:length(y2), ncol=ncol(y2), nrow=nrow(y2))
el = list(as.numeric(rbind(mat, NA)), as.numeric(rbind(t(mat), NA)))
w = rep(1, length(y2))
w[(nrow(y2)*(ncol(y2)-2) + 1):length(y2)] = 0.5
ew = lapply(el, function(v) rep(1,length(v)))
ew[[1]][53:65] = 0
ewL = rep(1, nrow(y2)*(ncol(y2)-1) + ncol(y2)*(nrow(y2)-1))
ewL[((nrow(y2)-1)*(ncol(y2)-1)+1):(ncol(y2)*(nrow(y2)-1))] = 0

emat = cbind(1:nrow(y2) + (ncol(y2)-2)*nrow(y2),
              1:nrow(y2) + (ncol(y2)-1)*nrow(y2))
E = Matrix::spMatrix(nrow(emat), length(y2), rep(1:nrow(emat),2), emat,
             c(rep(1,nrow(emat)),rep(-1,nrow(emat))))

lams = fp$lambda
for (thisLam in lams[-length(lams)]) {
  fl1 = fusedLattice(y2, w, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="prox", edgeWeights=ewL, E=E)
  fl2 = fusedLattice(y, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="prox")
  fg1 = fusedGraph(y2, w, edges=el, lambda=thisLam, eps=1e-08,
                      maxIter=250, method="prox", edgeWeights=ew, , E=E)

  spath = matrix(coef(fp, lambda=thisLam)$beta, nrow=nrow(y), ncol=ncol(y))
  fl1_beta = matrix(fl1$beta, nrow=nrow(y2), ncol=ncol(y2))[,-ncol(y2)]
  fl2_beta = matrix(fl2$beta, nrow=nrow(y), ncol=ncol(y))
  fg1_beta = matrix(fg1$beta, nrow=nrow(y2), ncol=ncol(y2))[,-ncol(y2)]

  expect_true(max(abs(spath - fl1_beta)) < EPS*10)
  expect_true(max(abs(spath - fl2_beta)) < EPS*10)
  expect_true(max(abs(spath - fg1_beta)) < EPS*10)
}

# (11) Run fused lasso with genlasso over a simple K3 graph;
#      compare to a tetris-like variant over a 2D lattice:
#
#          1 1 2 2
#          1 1 1 2
#          3 3 3 3
#
#      Weights and constraints are hand-coded.
set.seed(2)
nb_mat = matrix(c(1,1,2,2,3,3),ncol=2)
gr = graph.edgelist(nb_mat, directed=FALSE)

y1 = runif(3)
index = c(1,1,3,1,1,3,2,1,3,2,2,3)
y2 = matrix(y1[index],3,4)
w = c(1/5,1/3,1/4)[index]

ew = c(NA,1/3,NA,1/3,1/3,1/3,NA,1,NA,1/3,NA,NA,NA,1/3,NA,NA,NA)
ew[is.na(ew)] = 1

E = Matrix::spMatrix(9, 12, rep(1:9,each=2), c(1,2,2,4,4,5,5,8,3,6,6,9,9,12,7,10,10,11), rep(c(1,-1),9))

fp = fusedlasso(y1, graph=gr)

lams = fp$lambda
for (thisLam in seq(lams[1], 0, length.out=100)) {
  fl1 = fusedLattice(y2, weights=w, edgeWeights=ew,
                      lambda=thisLam, eps=1e-8,
                      maxIter=250, E=E)
  spath = as.numeric(coef(fp, lambda=thisLam)$beta)
  fl1_beta = tapply(fl1$beta, c(1,1,3,1,1,3,2,1,3,2,2,3), median)
  expect_true(max(abs(spath - fl1_beta)) < EPS)
}
