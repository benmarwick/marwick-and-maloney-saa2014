
# compute Principal Component Analysis
coords_D <- pca(coords_F)
dimnames(coords_D$coe)[[1]] <- as.character(fac_d$nms)
row.names(coords_D$li) <- as.character(fac_d$nms)
# coords_D$fac <- as.factor(fac$class)

# how much variation explained by each axes?
summary(coords_D)
screeplot(coords_D) #  represents the amount of inertia (usually variance) associated to each dimension.

# plot
dudi.plot(coords_D)
# many other types of plot

# observe how to change the point labelling
dudi.plot(coords_D, title="coords_D with no class but with ellipses")
dudi.plot(coords_D, 1, title="coords_D with no class but with ellipses")

## Now explore the _many_ kinds of plot that can be made with the PCA output...

dudi.plot(coords_D, 1, ellipses=FALSE, neighbors=TRUE, 
          shapes=FALSE, star=FALSE, col.nei="black",
          title="PCA with Gabriel's neighboring graph")

# this one looks good, I think:
dudi.plot(coords_D, labels = TRUE, points=FALSE, boxes=FALSE, shapes=TRUE, pos.shp="li", clabel = 0.75, title="PCA with labels and reconstructed shapes")

dudi.plot(coords_D, 1, points=FALSE, labels=TRUE,
          boxes=FALSE, shapes=FALSE,
          title="PCA with labels and ellipse")

dudi.plot(coords_D, 1, arrows=TRUE, dratio.arrow=0.2, shapes=FALSE,
          title="PCA with harmonic correlations")

# even more plots...
morpho.space(coords_D)
title("Default")

morpho.space(coords_D, nr.shp=3, nc.shp=3,
             col.shp="#1A1A1A22", border.shp=NA, pch.pts=5)
title("Custom nb of shapes")

morpho.space(coords_D, pos.shp="li")
title("Using the PC coordinates")

morpho.space(coords_D, pos.shp="circle")
title("A circle of shapes")

morpho.space(coords_D, xlim=c(-0.25, 0.25), ylim=c(-0.2, 0.2),
             pos.shp=expand.grid(seq(-0.2, 0.2, 0.1) , seq(-0.2, 0.2, 0.1)),
             border.shp="black", col.shp="#00000011", col.pts="firebrick")
title("Custom shape positions")

morpho.space(coords_D, xax=3, yax=4)
title("Morpho space on PC3 and PC4")

morpho.space(coords_D, amp.shp=3)
title("3 times magnif. differences")

morpho.space(pca(coords_R), rotate.shp=pi/6)
title("Using rFourier")

# display the contribution of Principal Component to shape description.
PC.contrib(coords_D)
PC.contrib(coords_D, PC.r=1:3, sd=3,
           cols=paste0(col.sari(3), "55"), borders=col.sari(3))


## stop here, below is experimental, incomplete or not working... ##
# 
# # http://biologyr.com/2013/06/30/como-colorear-por-grupo-un-morfoespacio-en-momocs-how-to-color-by-group-a-morphospace-in-momocs/
fix(morpho.space) # update function
# inser this after line 63
for (i in 1:dim(shapes)[3]) {
  coo.draw(shapes[,,i], points=FALSE,
           border=border.shp, col=col.shp[i], first.point=first.point)}
#
rgb(green=1,0,0, alpha=0.5)->green_trans 
rgb(blue=1,0,0, alpha=0.5)->blue_trans
rgb(red=1,0,0, alpha=0.5)->red_trans
# 
# # just proof of concept, I don't really want to group those last three together...
q <- character(length(coords_D$fac$class))
q[coords_D$fac$class == "flake"] <- green_trans
q[coords_D$fac$class == "retouchedflake"] <- red_trans
# q[coords_D$fac$class== "blade"] <- blue_trans
# q[coords_D$fac$class %in% c("unifacial", "retouched", "RT")] <- blue_trans
# 
# # colour by artefact class
morpho.space(coords_D, pos.shp = c("li"), rotate.shp=3.3, col.pts=0,scale.shp=0.06, 
              col.shp=q, border=0)
# 
# # http://biologyr.com/2013/09/14/una-forma-facil-de-ponerle-transparencia-a-los-colores-de-un-morfoespacio-en-r-an-easy-way-to-put-transparency-to-the-colors-of-a-morphospace-in-r/
# # We put transparency to the colors we want.
# # The closer to 0 is the second argument in the following lines ,
# # More transparent the color you get ( which comes from the first argument) 
library(scales) # for alpha function
alpha ( " orange" , 0.25) -> orange_trans
alpha ( "red" , 0.25) -> red_trans
alpha ("green " , 0.25) -> green_trans

# # roughly, create a vector of colors by assigning a color to each specimen
q <- character(length(coords_D$fac$class))
q[coords_D$fac$class == "flake"] <- red_trans
q[coords_D$fac$class == "retouchedflake"] <- green_trans
# q[coords_D$fac$class== "blade"] <- green_trans
# q[coords_D$fac$class %in% c("unifacial", "retouched", "RT")] <- green_trans
# 
# # Get the analysis of Fourier descriptors and make them componetes analysis
# eFourier ( coords_D , norm = T, nb.h = 50, smooth.it = 100) -> efou_dorsal
# pca_efou_dorsal <- pca ( efou_dorsal )
# 
# # Now plot the morphospace with some transparent colors
# # Color vector giving the " q " we created above
morpho.space( coords_D , col.pts = 0 , pos.shp = c ( "li" ) ,
               rotate.shp = 1.6, scale.shp = 0.05 , col.shp = q, border = " white")
# # works a treat