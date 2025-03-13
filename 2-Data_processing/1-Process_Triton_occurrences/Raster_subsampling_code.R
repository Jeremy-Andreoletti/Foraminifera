library(raster)

# code for extracting water depths
cutoff.rs <- brick(paste0("Outputs/Ma_gcm_interp_Ceno/scot_extrap_Ceno_", age, ".tif"))
cutoff.rs[cutoff.rs < 200] <- NA # exclude those shallower than 200m
plot(cutoff.rs)

# code for subsampling (from Erin)
p <- data.frame(rasterToPoints(cutoff.rs))
#assuming that you want to get 1000 sample points
sampled <- p[sample(nrow(p), 1000),]
plot(sampled$x, sampled$y)


