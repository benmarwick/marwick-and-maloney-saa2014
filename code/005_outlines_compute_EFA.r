#
# Perform Fourier analysis
coords_F <- eFourier(coords_center, nb.h=nb.h) # elliptical - main one
coords_R <- rFourier(coords_center, nb.h=nb.h) # radii variation
coords_T <- tFourier(coords_center, nb.h=nb.h) # tangent

# put the names and classifications on again
slot(coords_F, 'fac') <- fac
slot(coords_F, 'names') <- as.character(fac$nms)

# explore contribution of harmonics
hcontrib(coords_F)
hcontrib(coords_R)
hcontrib(coords_T)
