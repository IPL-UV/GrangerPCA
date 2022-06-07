source("functions.R")
library(ggplot2)
library(reshape2)
library(xtable)
library(fastICA)
library(ropls)

################# Example 1 ####################################################
######## Simulate data from the following VAR model
########
######## Y(t)   = sin(t / 10) + e
######## Z1(t)  = cos(3*t/10) + f
######## Z2(t)  = cos(5*t/100) + g
######## X1(t)  = 0.5 * X1(t - 2) + 0.1*Y(t - 1) + 0.2*t(t - 2) + 0.05*Y(t - 3) +
########          0.1*Z1(t-1) + 0.1*Z1(t-2) + 0.5*X9(t - 3) + e1
######## X2(t)  = 0.5 * X2(t - 2) + 0.1*Y(t - 1) + 0.2*t(t - 2) + 0.05*Y(t - 3) +
########          0.1*Z1(t-1) + 0.1*Z1(t-2) + 0.5*X10(t - 3) + e2
######## X3(t)  = 0.5*Z1(t-1) + e3
######## X4(t)  = 0.5*Z1(t-1) + e4
######## X5(t)  = 0.2*Z1(t-2) - 0.9*X1(t-1) + e5
######## X6(t)  = 0.2*Z1(t-2) - 0.9*X2(t-1) + e6
######## X7(t)  = -0.2*Z2(t-3) + 0.6 * X5(t - 2) + e7
######## X8(t)  = -0.2*Z2(t-3) + 0.6 * X6(t - 2) + e8
######## X9(t)  = 2*Z1(t-2) + 3*Z2(t-1) + e9
######## X10(t) = 2*Z1(t-2) + 3*Z2(t-1) + e10
######## where e, f, g, e1, ..., e10 are iid standard normal noises.
########
########       Y------>(X1, X2) ----------------> (X5, X6)
########                  ^                       ^    \
########                  |       _______________/      \
########                  |      /                       \
########  ----Z1--------------------> (X3, X4)            \
########  |                                                V
########  |   Z2--------------------------------------> (X7, X8)
########  |    |
########  V    V
########  (X9, X10)
################################################################################
set.seed(2022)
#dx <- 10  # low-dimension
dx <- 500   # high-dimension
dxc <- dx / 5
N <- 1000

Y <- matrix(rnorm(N, mean = sin(0.2 * (1:N)) + cos(0.1 * (1:N))), nrow = N)
Z1 <- matrix(rnorm(N, mean = cos(0.4 * (1:N))), nrow = N)
Z2 <- matrix(rnorm(N, mean = cos(0.05 * (1:N))), nrow = N)

X <- matrix(rnorm(dx * N, sd = 1), nrow = N, ncol = dx)
f1 <- function(xx, yy, zz) 0.2*cos(7*xx) + 0.3*sin(4*xx)*cos(6*yy) + 0.5*sin(3*zz)*sin(20*xx)
f2 <- function(xx) -0.3*sin(10*xx)*cos(5*xx)
f3 <- function(xx, yy, zz) -0.3*sin(xx)*cos(3**yy) + 0.5*sin(5*yy)*sin(7*zz) 
f4 <- function(xx, yy, zz) 0.4*sin(8*xx) - 0.2*sin(7*zz)*cos(zz)*cos(8*yy)

for (i in 4:nrow(X)){
  X[i, 1:dxc] <- X[i, 1:dxc] + #A
     0.2 * X[i - 2, 1:dxc] + f1(Y[i-1,],Y[i-2,],Y[i-3]) + 
    0.2 * Z1[i - 1, ] + 0.2 * Z1[i - 2, ]
  X[i, (dxc+1):(2*dxc)] <- X[i, (dxc+1):(2*dxc)] + #B
    f2(X[i - 1, (dxc+1):(2*dxc)]) + 0.5 * Z1[i - 1, ]
  X[i, (2*dxc+1):(3*dxc)] <- X[i, (2*dxc+1):(3*dxc)] + #C
     0.8 * Z1[i - 2, ] + 
      f3(X[i - 2, (2*dxc+1):(3*dxc)], mean(X[i - 1, 1:dxc]), mean(X[i - 2, 1:dxc]))
  X[i, (3*dxc+1):(4*dxc)] <- #D
    X[i, (3*dxc+1):(4*dxc)]  - 0.7 * Z2[i - 3, ]  +
    f4(X[i - 2, (3*dxc+1):(4*dxc)], Z2[i - 3, ], mean(X[i - 2, (2*dxc+1):(3*dxc)]))
  X[i, (4*dxc+1):(5*dxc)] <- #E
    X[i, (4*dxc+1):(5*dxc)] + 0.1 * X[i - 1, (4*dxc+1):(5*dxc)] + 0.8*Z1[i - 2, ] + 0.6*Z2[i - 1, ] + 0.2*Z2[i - 2, ]
}


### apply all methods
maxlag <- 5
PC <- prcomp(X, retx = TRUE, scale. = TRUE, center = TRUE,rank. = 5)
kPC <- kpca(X, features = 5, kpar = list(sigma = 1.5 / dx))

GPC <- gpca(PC$x[,1:5], Y, maxlag = maxlag, scale. = FALSE, center = FALSE)
GPC$rotation <- PC$rotation[,1:5] %*% GPC$rotation

kGPCA <- gpca(kPC@rotated, Y, maxlag = maxlag, scale. = FALSE, center = FALSE)

system.time(TLCC <- tlcancor(X, Y, maxlag))

system.time(TLOPLS <- tlopls(X, Y, maxlag, rank = 10))

resall <- list(PCA = PC, 
               #GC = GC, 
               "GPCA" = GPC,
               "kGPCA" = kGPCA,
               TLCC = TLCC,
               TLPLS = TLOPLS#,
               #ICA = ICA
               )

################  compute pvals for granger test
#pvals <- data.frame(lapply(resall, function(res){
#  apply(res$x[,1:5], MARGIN = 2, function(xx) lmtest::grangertest(Y, xx, order = 3)$`Pr(>F)`[2])
#}))

pvals <- data.frame(lapply(resall, function(res){
  apply(res$x[,1:5], MARGIN = 2, function(xx) granger_test(Y, xx, res$x[,1:5],3)$`Pr(>F)`[2])
}))

data.frame(lapply(resall, function(res){
  apply(res$x[,1:5], MARGIN = 2, var)
}))

DS <- lapply(resall, function(res){
  D <- melt(scale(res$rotation[,1:5], center = FALSE))
  D$Var2 <- as.numeric(D$Var2)
  D
})

DD <- melt(DS, measure.vars = "value")
DD$Var1 <- as.factor(DD$Var1)
DD$Var2 <- as.factor(DD$Var2)

cplot <- ggplot(DD) + 
  geom_tile(aes(x = Var2, y = Var1, fill = value))  + 
  scale_y_discrete(breaks = seq(from = dxc / 2, to = dx, by =dxc), labels = c("A", "B", "C", "D", "E"))+ 
  ylab("original coordinates") + xlab("components") + 
  scale_fill_gradient2(low = "darkred", mid =  "white", high = "darkblue", limit = c(-2,2), oob = scales::squish) + 
  facet_grid(rows = vars(L1)) + 
  theme_bw() + theme(legend.position = "none", legend.title = element_blank())

cplot
ggsave(cplot, file = paste0(dx,"components_example1_nonlin.pdf"), width = 3.5, height = 4, units = "in")

########### time-series plot
tpl <- 150
DD2 <- melt(data.frame(
  #PC1 = -scale(PC$x)[1:tpl,1],
  #PC2 = scale(PC$x)[1:tpl,2],
  "GPCA" = sign(cor(GPC$x[,1], X[,1])) * scale(GPC$x)[1:tpl,1],
  #"GC" = -scale(GC$x)[1:tpl,1],
  "TLCC" = sign(cor(TLCC$x[,1], X[,1])) *  scale(TLCC$x)[1:tpl,1],
  "TLPLS" = sign(cor(TLOPLS$x[,1], X[,1])) * scale(TLOPLS$x)[1:tpl,1],
  #"ICA1" = -ICA$x[1:100,1],
  time = 1:tpl, 
  "X_m" = rowMeans(scale(X)[1:tpl, 1:dxc])),
  id.vars = c("time"))

signalplot <- ggplot(DD2) + 
  geom_line(aes(y = value, x = time, color = variable )) + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.box.spacing = unit(0,'points'),
        axis.title.y = element_blank())
signalplot
ggsave(signalplot, file = paste0(dx,"signals_example1_nonlin.pdf"), width = 3.5, height = 2, units = "in")

### save p-vals table
rownames(pvals) <- paste0("C", 1:5)

print(xtable(pvals, caption = "Granger causality test p-values for 
             each component of the different methods.", label = "tab:pvals",
             display = c("s", rep("e",length(resall)))), booktabs = TRUE, file = paste0(dx,"pvals_ex1_nonlin.tex"))
