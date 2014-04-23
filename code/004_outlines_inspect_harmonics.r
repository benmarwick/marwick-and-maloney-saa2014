# determine how many harmonics to use in
# the elliptical fourier analysis
# now just working through http://www.vincentbonhomme.fr/Momocs/vignettes/A-graphical-introduction-to-Momocs.pdf

nb.h = 50 # over estimate harmonics to see what minimum number is ok...

# one method...
# windows() # pops up a new window to hold plots
# set up a grid of plots...
# par(mfrow=c(ceiling(length(fac$nms)/2),ceiling(length(fac$nms)/2)-1)) 
# they may not all fit, if so, skip on to the next approach...
# lapply(coords_center@coo, function(i) hpow(Coo(i), nb.h = nb.h))
# dev.off() # may need to repeat this a few times to free up the graphics device

# another approach... watch the colours flash by...
# lapply(coords_center@coo, function(i) hqual(Coo(i), harm.range = seq(1, nb.h, 10),  plot.method = c("panel")[1]))
# lapply(coords_center@coo, function(i) hqual(Coo(i), harm.range = seq(1, nb.h, 10),  plot.method = c("stack")[1]))

# yet another approach, probably the best one since it's more objective and repeatable 
# hpow(coords_scale, nb.h = nb.h,
#     title="eFourier with extrema and mean dev.")
hpow(coords_scale, probs=c(0.25, 0.5, 0.75), drop=FALSE, legend=TRUE, nb.h = nb.h,
     title="eFourier three quartiles")
# hpow(coords_scale, method="rfourier",
#     title="rFourier")
# hpow(coords_scale, method="tfourier",
#     title="tFourier")

# inspect harmonics plots and output to see minimum number needed to get 0.9999 of total 
# harmonic power.

nb.h = 30 # after looking at the plots
