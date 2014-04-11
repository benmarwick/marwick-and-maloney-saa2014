# from https://gist.githubusercontent.com/benmarwick/6260541/raw/175ea076e29fd0064603fa00f5cce46fedc3b3e5/artefact-morpho.R

# This script is a workflow for analysing 2D artefact outlines from 3D 
# scan objects captured by NextEngine and ScanStudio. Part of the process
# occurs in ScanStudio and GIMP and the quantative analysis of the outlines
# is done in R. In GIMP we process the images into B&W silhouettes ready for R. 
# In R we do elliptical fourier analysis to summarise the image outlines 
# then PCA and MANOVA to discriminate between shape variations
# and test for differences

########### Using GIMP to preparing ScanStudio images for input into R ###########
## Get screenshot from scanstudio into GIMP
# Edit > Paste Into  (image > check canvas size)
## Crop away background
# Click on the rectangle select tool, select area around artefact
# leave a bit of space between the arteafact and the crop rectangle, then
# Image > Crop to selection
## convert to B&W by ... (http://www.gimp.org/tutorials/Color2BW/)
# (other possible methods: http://biologyr.com/2013/04/22/como-binarizar-imagenes-a-blanco-y-negro-con-el-programa-gratis-gimp-how-to-binarize-images-to-black-and-white-with-the-free-program-gimp/)
# Colours > Threshold... then increase min value until infill is solid black 
# (typically over 200, a few bits of black in the background are ok)
## Clean up background
# Click on the paintbrush tool, select white colour (HTML notation: ffffff) 
# and manually paint over bits of black in the background and 
# any fuzzy bits near the edge of the artefact. Take care not
# to alter the outline of the artefact as this is what we're 
# analysing - the shape of the artefact. 
# 
# save as xcf (use original scanstudio filename but with no spaces in filename, use underscores instead)
# export to jpg (use original scanstudio filename but no spaces in filename, use underscores instead)
############ finish ########################################

############ now working in RStudio to analyse images ###### 
# setup: make sure you have
# installed http://www.imagemagick.org/script/binary-releases.php#windows
# I am using ImageMagick-6.8.8-6-Q16-x64-dll.exe (with all install options)
# and restarted RStudio

# Much of this is based on code in Claude J. 2008. Morphometrics with R. Springer, New York.
# But he hasn't made a package, so we'll just source his functions
# get functions from Claude's website (you can download this txt file
# and source from your disk instead)
source("F:/My Documents/My Papers/conferences/SAA2014/Marwick and Malony/code/claude_2008.r")
# source("http://www.isem.cnrs.fr/IMG/txt/Rfunctions.txt") # author's copy
# book errata: http://www.isem.univ-montp2.fr/recherche/teams/developmental-biology-and-evolution/staff/claudejulien/?lang=en
# source("C:\\Users\\marwick\\Downloads\\Rfunctions.txt") # my local copy
# library for greyscale stuff
# you may need to do this: install.packages("pixmap")
library(pixmap)

# set working dir (where the jpgs are)
# make sure that all the jpgs in this folder have
# filenames that are consistently patterned
# withe underscores separating the terms, exactly like this:
# "type_subtype_region_id.jpg 
# we're going to use the 'type' bit as a grouping factor
# setwd("C:\\Users\\marwick\\UW-Dropbox\\Dropbox\\Malony 's lithic scan data/japi scans/")

# make a list of jpgs to work on
# list_imgs <- list.files(getwd(), pattern = ".jpg$") # not really a list, but a chr vector...
# inspect list to make sure nothing dodgy is there
list_imgs
# classify the list for later analysis, we'll just get the
# word before the first underscore since that's a 'type' word
class <- unname(sapply(list_imgs, function(i) unlist(strsplit(i, split="_"))[1]))
# get vector of artefact IDs
nms <- gsub(".jpg", "", list_imgs)
# make a data frame of classifications and names for later use
fac <- data.frame(class = class, nms = nms)
# inspect to ensure it's going as planned...
fac

# exlcude classes with less than two members because they make mean values
# impossible to compute and stop the thin plate spline analysis from working
library(dplyr)
# find number of members of each class
tbl <- summarise(group_by(fac, class), count = length(class))
# if there are classes with <2 members, do this, otherwise skip
ifelse(length(tbl$class < 2)
# any classes have <2 members?
dropme <- filter(tbl, count < 2)$class
if(length(dropme) != 0){
# drop that/those class(es) from the table 
fac <- filter(fac, class != dropme)
# and from the list of images (by -ve indexing)
list_imgs <- list_imgs[-which(fac$class %in% dropme)]
} else {
  temp <- NULL  # do nothing
}


# Now loop over each jpg in the list to trace its outline
# and extract coordinates of the outline
# This is a bit slow, but plots pop up and down to show
# that's working (can be commented out to speed up loop)

# create an object to store the results
outlines <- vector("list", length = length(list_imgs))

# here's the loop
for(i in 1:length(list_imgs)) {
  
  # from Claude p. 40, convert image to ppm
  # could work with any image format since
  # imagemagick can handle >100 formats
  img <- nms[i] 
  # this is a call to the windows command line
  shell(shQuote(paste0("convert ", img, ".jpg ", img, ".ppm")))
  M<- read.pnm(paste0(img, ".ppm"))
  plot(M)
  
  # Convert the RGB image to a gray-scale image. 
  M<- as(M, "pixmapGrey")
  
  # Binarize the image (Claude p. 48)
  
  M@grey[which(M@grey>=0.9)]<-1
  M@grey[which(M@grey<0.9)]<-0.7
  plot(M)
  
  # auto-trace outline of shape
  # start<-locator(1) # manual option: now click on shape to pick starting point
  # or make a guess of where to start
  start <- list(x = 100, y = 100) # first guess
  t0 <- try(Rc<-Conte(c(round(start$x),round(start$y)),M@grey))
  if("try-error" %in% class(t0)) {  # if first guess is no good, make a second guess
    start <- list(x = 200, y = 200)  
    Rc<-Conte(c(round(start$x),round(start$y)),M@grey)
  } else {
    Rc<-Conte(c(round(start$x),round(start$y)),M@grey)
  }
  lines(Rc$X, Rc$Y, lwd=4)
  
  # this is the outline, we can inspect it...
  str(Rc) # inspect data, just for reassurance 
  # draw plot
  plot(Rc$X, Rc$Y, type = "l", asp = TRUE)
  
  # store result in list for further analysis
  outlines[[i]] <- Rc
  
  # for a long loop, it's nice to see what step you're up to
  print(i)
}

# quick look to see what we got...
str(outlines)
# did each image get an outline?
length(outlines) == nrow(fac) # should be TRUE

# Now for statistical analysis of the outlines...
# I use Momocs2, which is a copy of the package that I've made 
# some changes to (it's only on my hard drive). You can use Momocs
# from CRAN (the usual repository), it should work the same for 
# these functions...
library(Momocs) # more details: http://www.vincentbonhomme.fr/Momocs/vignettes/A-graphical-introduction-to-Momocs.pdf
# many excellent uses here: http://biologyr.com/tag/momocs/ (some adapted here below)

# convert outlines to 'Coo' format
coords <- Coo(lapply(outlines, function(i) (simplify2array(i))))
# center the coords
coords_center <- Coo(lapply(coords@coo, coo.center))

# scale the coords
coords_scale <- Coo(lapply(coords_center@coo, function(i) coo.scale(i, 1)))

# give the images their file names and types 
slot(coords_scale, 'fac') <- fac
slot(coords_scale, 'names') <- as.character(fac$nms)

# quick plot
panel(coords_scale,  borders="black", cols="grey90") # no names
panel(coords_scale,  borders="black", names=TRUE, cols="grey90") # with names

# determine how many harmonics to use in
# the elliptical fourier analysis
# now just working through http://www.vincentbonhomme.fr/Momocs/vignettes/A-graphical-introduction-to-Momocs.pdf

nb.h = 50 # over estimate harmonics to see what minimum number is ok...

# one method...
windows() # pops up a new window to hold plots
# set up a grid of plots...
par(mfrow=c(ceiling(length(fac$nms)/2),ceiling(length(fac$nms)/2)-1)) 
# they may not all fit, if so, skip on to the next approach...
lapply(coords_center@coo, function(i) hpow(Coo(i), nb.h = nb.h))
dev.off() # may need to repeat this a few times to free up the graphics device

# another approach... watch the colours flash by...
lapply(coords_center@coo, function(i) hqual(Coo(i), harm.range = seq(1, nb.h, 10),  plot.method = c("panel")[1]))
lapply(coords_center@coo, function(i) hqual(Coo(i), harm.range = seq(1, nb.h, 10),  plot.method = c("stack")[1]))

# yet another approach, probably the best one since it's more objective and repeatable 
hpow(coords_scale, nb.h = 50,
     title="eFourier with extrema and mean dev.")
hpow(coords_scale, probs=c(0.25, 0.5, 0.75), drop=FALSE, legend=TRUE, nb.h = 50,
     title="eFourier three quartiles")
hpow(coords_scale, method="rfourier",
     title="rFourier")
hpow(coords_scale, method="tfourier",
     title="tFourier")

# inspect harmonics plots and output to see minimum number needed to get 0.9999 of total 
# harmonic power.

nb.h = 43 # for my batch of images

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

# compute Principal Component Analysis
coords_D <- pca(coords_F)
dimnames(coords_D$coe)[[1]] <- as.character(fac$nms)
# plot
dudi.plot(coords_D)
# many other types of plot

# observe how to change the point labelling
dudi.plot(coords_D, title="coords_D with no class but with ellipses")
dudi.plot(coords_D, 1, title="coords_D with no class but with ellipses")
dudi.plot(coords_D, 2, title="coords_D with no class but with ellipses")

## Now explore the _many_ kinds of plot that can be made with the PCA output...

dudi.plot(coords_D, 1, ellipses=FALSE, neighbors=TRUE, 
          shapes=FALSE, star=FALSE, col.nei="black",
          title="coords_D with Gabriel's neighboring graph")

# this one looks good, I think:
dudi.plot(coords_D, 1, labels=FALSE, points=FALSE, boxes=FALSE,
          shapes=TRUE, pos.shp="li",  
          title="coords_D with labels and reconstructed shapes")

dudi.plot(coords_D, 1, points=FALSE, labels=TRUE,
          boxes=FALSE, shapes=FALSE,
          title="coords_D with labels and ellipse")

dudi.plot(coords_D, 1, arrows=TRUE, dratio.arrow=0.2, shapes=FALSE,
          title="coords_D with harmonic correlations")

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


# Statistical test differences between groups 
# through a Multivariate Analysis of Variance.
(mnv <- manova.Coe(coords_F, "class", retain = 5))


# inspect thin plate splines
x <- meanShapes(coords_F, 'class')
str(x) # what we got
stack(Coo(x), borders=col.gallus(2))
# compare a few classes
tps.grid(x[[3]],  x[[4]]) # have a look to ?"["
tps.grid(x$flake, x$retouched, grid.outside=2,
         grid.size=80, amp=2, plot.full=FALSE)


## stop here, below is experimental, incomplete or not working... ##

# http://biologyr.com/2013/06/30/como-colorear-por-grupo-un-morfoespacio-en-momocs-how-to-color-by-group-a-morphospace-in-momocs/
fix(morpho.space) # update function
# inser this after line 63
for (i in 1:dim(shapes)[3]) {
  coo.draw(shapes[,,i], points=FALSE,
           border=border.shp, col=col.shp[i], first.point=first.point)}
#
rgb(green=1,0,0, alpha=0.5)->green_trans 
rgb(blue=1,0,0, alpha=0.5)->blue_trans
rgb(red=1,0,0, alpha=0.5)->red_trans

# just proof of concept, I don't really want to group those last three together...
q <- character(length(coords_D$fac$class))
q[coords_D$fac$class == "biface"] <- green_trans
q[coords_D$fac$class == "flake"] <- red_trans
q[coords_D$fac$class== "blade"] <- blue_trans
q[coords_D$fac$class %in% c("unifacial", "retouched", "RT")] <- blue_trans

# colour by artefact class
morpho.space(coords_D, pos.shp = c("li"), rotate.shp=3.3, col.pts=0,scale.shp=0.06, 
             col.shp=q, border=0)

# http://biologyr.com/2013/09/14/una-forma-facil-de-ponerle-transparencia-a-los-colores-de-un-morfoespacio-en-r-an-easy-way-to-put-transparency-to-the-colors-of-a-morphospace-in-r/
# We put transparency to the colors we want.
# The closer to 0 is the second argument in the following lines ,
# More transparent the color you get ( which comes from the first argument) 
library(scales) # for alpha function
alpha ( " orange" , 0.25) -> orange_trans
alpha ( "red" , 0.25) -> red_trans
alpha ("green " , 0.25) -> green_trans

# roughly, create a vector of colors by assigning a color to each specimen
q <- character(length(coords_D$fac$class))
q[coords_D$fac$class == "biface"] <- orange_trans
q[coords_D$fac$class == "flake"] <- red_trans
q[coords_D$fac$class== "blade"] <- green_trans
q[coords_D$fac$class %in% c("unifacial", "retouched", "RT")] <- green_trans

# Get the analysis of Fourier descriptors and make them componetes analysis
eFourier ( coo_dorsal , norm = T, nb.h = 50, smooth.it = 100) -> efou_dorsal
pca_efou_dorsal <- pca ( efou_dorsal )

# Now plot the morphospace with some transparent colors
# Color vector giving the " q " we created above
morpho.space( coords_D , col.pts = 0 , pos.shp = c ( "li" ) ,
              rotate.shp = 1.6, scale.shp = 0.05 , col.shp = q, border = " white")
# works a treat


# cluster analysis and dendrograms
library(MASS)
library(shapes)

# interpolate point so each shape has the same number
interp <- Coo(lapply(outlines, function(i) coo.sample.int(simplify2array(i), 3000) ))

gorf<-gorf.dat
panf<-panf.dat[c(5,1,2:4,6:8),,]
pongof<-pongof.dat[c(5,1,2:4,6:8),,]
APE<-array(c(panf, gorf, pongof),dim=c(8,2,80))
APE<-array(c(panf, gorf),dim=c(8,2,56))
AP<-orp(pgpa(APE)$rotated)
m<-t(matrix(AP, 16, 56))
par(mar=c(0.5,2,1,1))
layout(matrix(c(1,2),2,1))
plot(hclust(dist(m), method="average"),main="UPGMA"
     ,labels=c(rep("P",26),rep("G",30)),cex=0.7)
plot(hclust(dist(m), method="complete"),main="COMPLETE"
     ,labels=c(rep("P",26),rep("G",30)),cex=0.7)



# from package vignette...
botFg <- meanShapes(botF)
str(botFg)
stack(Coo(botFg), borders=col.gallus(2))
tps.grid(botFg[[1]], botFg[[2]]) # have a look to ?"["
tps.grid(botFg$beer, botFg$whisky, grid.outside=2,
         grid.size=80, amp=3, plot.full=FALSE)

set.seed(007)
my_list<- lapply(1:10, function(i) matrix(data = runif(sample(seq(2,10,2), 1)), ncol = 2))

my_array<-Map(function(x,y,z) array(n,c(x,y,z)), myrow, mycol,myz) # but this creates the list of 1 array
