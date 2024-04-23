setwd("./") # set to the path of ConjugateNNGP

rm(list = ls())
load("./nngp_10.RData")
source("./functions.R")
library(foreach)
library(doParallel)
library(parallel)
library(rbenchmark)
library(rbenchmark)
library(ggplot2)

phi <- 17.64706; deltasq <- 10
delta <- sqrt(deltasq)
a = as; b = bs;
ind_x <-c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1) : N), each = M)), 1:N)
ind_y <- c(c(t(NN.matrix$NN_ind))[which(c(t(NN.matrix$NN_ind)) > 0)], 1:N)
X.ord <- X[NN.matrix$ord, ]
Y.ord <- Y[NN.matrix$ord]
AD <- getADstan(neardist = NN.matrix$NN_dist,  neardistM = NN.matrix$NN_distM, 
                N = N, M = M, phi = phi) 
D <- AD[M + 1, ]
ind_x_X <- rep(1:N, P); ind_y_X <- rep(1:P, each = N)
ind_x_X_up <- c(ind_x_X, 1:N)
ind_y_X_up <- c(ind_y_X, (P + 1):(N + P))
X_star_up <- 
  sparseMatrix(ind_x_X_up, ind_y_X_up, x = c(as.vector(X.ord) / delta,
                                             rep(1 / delta, N)))
X_star_down <- 
  sparseMatrix(i = ind_x, j = (ind_y + P), 
               x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, N))) / 
  sqrt(D) 
X_star <- rbind(X_star_up, X_star_down)
Y_star <- c(Y.ord / delta, rep(0, N))
XTY_star <- crossprod(X_star, Y_star)
XTX_star <- crossprod(X_star, X_star)

#compare solving time between CG, symbolic symbolic decomposition and default solve for 1000 iteration#
solver <- new(CholeskySolver)
benchmark(solve(XTX_star, XTY_star[1:(N + P)]),
          cgsparse(XTX_star, XTY_star[1:(N + P)]),
          solver$solve(XTX_star, XTY_star[1:(N + P)]),replications = 1000)

#this parts shows that the result of CG, symbolic decomposition and default solve are the same#
gamma_hat <- cgsparse(XTX_star, XTY_star[1:(N + P)])
xtt=solve(XTX_star, XTY_star[1:(N + P)])
summary(gamma_hat-xtt)
all.equal(gamma_hat, xtt)
gamma_hat_dense <- as.matrix(gamma_hat)  
xtt_dense <- as.matrix(xtt)  
summary(gamma_hat_dense-xtt_dense)
all.equal(gamma_hat_dense, xttt_dense)
xttt=solver$solve(XTX_star, XTY_star[1:(N + P)])
summary(gamma_hat-xttt)
all.equal(gamma_hat, xttt)
