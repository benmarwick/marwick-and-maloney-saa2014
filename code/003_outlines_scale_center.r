# Now for statistical analysis of the outlines...
# I use Momocs2, which is a copy of the package that I've made 
# some changes to (it's only on my hard drive). You can use Momocs
# from CRAN (the usual repository), it should work the same for 
# these functions...
library(Momocs) # more details: http://www.vincentbonhomme.fr/Momocs/vignettes/A-graphical-introduction-to-Momocs.pdf
# many excellent uses here: http://biologyr.com/tag/momocs/ (some adapted here below)

# convert outlines to 'Coo' format
coords <- Coo(lapply(outlines_500, function(i) (simplify2array(i))))
# center the coords
coords_center <- Coo(lapply(coords@coo, coo.center))

# scale the coords
coords_scale <- Coo(lapply(coords_center@coo, function(i) coo.scale(i, 1)))

# give the images their file names and types
fac <- fac[-nm, ] # exclude names of shapes with <500 points
slot(coords_scale, 'fac') <- fac
slot(coords_scale, 'names') <- as.character(fac$nms)

# quick plot
panel(coords_scale,  borders="black", cols="grey90") # no names
panel(coords_scale,  borders="black", names=TRUE, cols="grey90") # with names
