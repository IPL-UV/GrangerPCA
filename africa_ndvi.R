library(ggplot2)
library(ggspatial)
library(reshape2)
library(raster)
library(ncdf4)
library(fastICA)
library(kernlab)
source("functions.R")

path <- "databases/data/africa_ndvi/"
filename_template <- paste0(path,"ndvi_TMSTMP.nc")

years <- 2007:2017
days <- formatC(seq(from =5, to = 365, by = 8), width = 3, flag = "0")
timestamps <- apply(expand.grid(days, years)[,2:1], 1, paste0, collapse = "")

ndvi.list <- lapply(timestamps, function(tmstmp) raster(sub("TMSTMP", tmstmp, filename_template)))

ndvi <- stack(ndvi.list)  

LC <- raster("../LatentGranger/databases/data/LC_hd_global_2012_full.tif")
LC <- crop(LC, extent(-2*20 + 360, 2*55 + 360, -2*37 + 180, 2*38 + 180))
#LC <- aggregate(LC, fact = 3, fun = median)

extent(LC) <- extent(ndvi)
ndvi <- raster::mask(ndvi, LC, maskvalue = 0)
#ndvi <- raster::mask(ndvi, LC, maskvalue = 16)

ndvi.array <- as.array(ndvi)

X <- apply(ndvi.array, c(1,2), c)
X <- matrix(c(X), nrow = dim(X)[1])
X <- t(na.omit(t(X))) ##omit columns all NA


enso <- read.table("../LatentGranger/databases/data/inino34_daily.dat")
enso$smoothed <- filter(x = enso$V2, filter =  rep(1/8, 8), method = "convolution", sides = 1)

enso$V1 <- paste0(enso$V1)
Y <- c()
for (year in years){
  whr <- stringr::str_starts(enso$V1, paste0(year))
  dd <- (1:365)[as.numeric(days)]
  Y <- c(Y, enso$smoothed[whr][dd])
}


Y <- matrix(Y, nrow = length(Y)) 

dim(Y)
dim(X)

maxlag <- 10
rank <- 10


#X <- X %*% diag( 1 / sqrt(colSums(X^2)))

PCA <- prcomp(X , scale. = TRUE, center = TRUE, rank = rank)

#VM <- list(rotation = varimax(PCA$rotation)$loadings)
#VM$x <- scale(X) %*% VM$rotation

D <- dist(scale(X))^2
bw_median <- sqrt(median(D)/2)
sigma_median <- sqrt(1/(2 * median(D)))
kPCA <- kpca(scale(X), kernel = "rbfdot",
            features = rank, kpar = list(sigma = 0.05* sigma_median))

GPCA <- gpca(PCA$x, Y, maxlag = maxlag, scale. = FALSE, center = FALSE, 
            intercept = TRUE)


kGPCA <- gpca(kPCA@rotated, Y, maxlag = maxlag, scale. = FALSE, center = FALSE, 
            intercept = TRUE)

mask <- !is.na(ndvi.array[,,1])

ALLMTHDS <- list(
  "PCA" = PCA$x,
  "GPCA" = GPCA$x,
  "kPCA" = kPCA@rotated,
  "kGPCA" = kGPCA$x
)

pvals <- data.frame(lapply(ALLMTHDS, function(x){
  apply(x[,1:10], MARGIN = 2, function(xx) granger_test(Y, xx, x[,1:10],maxlag)$`Pr(>F)`[2])
})) 

#pvals <- data.frame(lapply(ALLMTHDS, function(x){
#  apply(x[,1:10], MARGIN = 2, function(xx) lmtest::grangertest(Y, xx, order = maxlag)$`Pr(>F)`[2])
#})) 

## build all maps 
allmaps <- lapply(ALLMTHDS, function(MTH) sapply(1:10, function(i){
  xx <- MTH[,i]
  # compute the sign correction
  sC <- c(sign(cor(Y, xx)))
  AA <- matrix(nrow = nrow(mask), ncol = ncol(mask), data = NA)
  AA[mask] <-  sC*coefficients(lm(X ~ MTH))[1+i,]
  #M <- max(abs(AA), na.rm = TRUE)
  #raster(x = 2 * ((AA + M)/ (2*M) - 0.5), template = ndvi)
  
  raster(x = AA / (max(AA, na.rm = TRUE) + max(-AA, na.rm = TRUE) ) , template = ndvi)
}))




maps <- lapply(allmaps, stack)

dfs <- lapply(maps, df_spatial)

df <- melt(dfs, id.vars = c('x', "y"))

#thr <- 5e-3
#df$value[df$value > thr] <- thr
#df$value[df$value < -thr] <- -thr
mapplot <- ggplot() + geom_raster(df, mapping = aes(x = x, y = y, fill = value)) + 
  scale_fill_gradient2(low = "red3", mid = "white", 
                       high = "blue3", na.value = "white", 
                       limits = c(-0.3,0.3), oob = scales::oob_squish) + 
  facet_grid(cols = vars(variable), rows = vars(L1)) + theme(axis.ticks = element_blank(), 
                                            axis.text = element_blank(), 
                                            axis.title = element_blank(), 
                                            legend.position = "none") +  
  ggtitle("")

mapplot
ggsave("maps_plot_africa.pdf", plot = mapplot, width = 6, height = 3, units = 'in')


tpl <- 506
DD2 <- melt(data.frame(
   PC3 = scale(PCA$x)[1:tpl,6],
  #"GPCA1" = sign(cor(GPCA$x[,1], Y[,1])) * scale(GPCA$x)[1:tpl,1],
  time = 1:tpl, 
  "ENSO3.4" = Y[1:tpl,1]),
  id.vars = c("time"))

signalplot <- ggplot(DD2) + 
  geom_line(aes(y = value, x = time, color = variable ), size = 1) + theme_bw() + 
  scale_color_manual(values = c("blue", "black")) + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.box.spacing = unit(0,'points'),
        axis.title.y = element_blank())
signalplot
ggsave(signalplot, file = paste0(dx,"signals_africa.pdf"), width = 3.5, height = 2, units = "in")



### save p-vals table
rownames(pvals) <- paste0("C", 1:10)

print(xtable(pvals, caption = "Granger causality test p-values for 
             each component of the different methods.", label = "tab:pvals",
             display = c("s", rep("e",length(ALLMTHDS)))), booktabs = TRUE, file = paste0("pvals_africa.tex"))
