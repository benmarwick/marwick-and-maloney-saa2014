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









## stop here, below is experimental, incomplete or not working... ##
# 
# # http://biologyr.com/2013/06/30/como-colorear-por-grupo-un-morfoespacio-en-momocs-how-to-color-by-group-a-morphospace-in-momocs/
# fix(morpho.space) # update function
# # inser this after line 63
# for (i in 1:dim(shapes)[3]) {
#   coo.draw(shapes[,,i], points=FALSE,
#            border=border.shp, col=col.shp[i], first.point=first.point)}
# #
# rgb(green=1,0,0, alpha=0.5)->green_trans 
# rgb(blue=1,0,0, alpha=0.5)->blue_trans
# rgb(red=1,0,0, alpha=0.5)->red_trans
# 
# # just proof of concept, I don't really want to group those last three together...
# q <- character(length(coords_D$fac$class))
# q[coords_D$fac$class == "biface"] <- green_trans
# q[coords_D$fac$class == "flake"] <- red_trans
# q[coords_D$fac$class== "blade"] <- blue_trans
# q[coords_D$fac$class %in% c("unifacial", "retouched", "RT")] <- blue_trans
# 
# # colour by artefact class
# morpho.space(coords_D, pos.shp = c("li"), rotate.shp=3.3, col.pts=0,scale.shp=0.06, 
#              col.shp=q, border=0)
# 
# # http://biologyr.com/2013/09/14/una-forma-facil-de-ponerle-transparencia-a-los-colores-de-un-morfoespacio-en-r-an-easy-way-to-put-transparency-to-the-colors-of-a-morphospace-in-r/
# # We put transparency to the colors we want.
# # The closer to 0 is the second argument in the following lines ,
# # More transparent the color you get ( which comes from the first argument) 
# library(scales) # for alpha function
# alpha ( " orange" , 0.25) -> orange_trans
# alpha ( "red" , 0.25) -> red_trans
# alpha ("green " , 0.25) -> green_trans
# 
# # roughly, create a vector of colors by assigning a color to each specimen
# q <- character(length(coords_D$fac$class))
# q[coords_D$fac$class == "biface"] <- orange_trans
# q[coords_D$fac$class == "flake"] <- red_trans
# q[coords_D$fac$class== "blade"] <- green_trans
# q[coords_D$fac$class %in% c("unifacial", "retouched", "RT")] <- green_trans
# 
# # Get the analysis of Fourier descriptors and make them componetes analysis
# eFourier ( coo_dorsal , norm = T, nb.h = 50, smooth.it = 100) -> efou_dorsal
# pca_efou_dorsal <- pca ( efou_dorsal )
# 
# # Now plot the morphospace with some transparent colors
# # Color vector giving the " q " we created above
# morpho.space( coords_D , col.pts = 0 , pos.shp = c ( "li" ) ,
#               rotate.shp = 1.6, scale.shp = 0.05 , col.shp = q, border = " white")
# # works a treat
# 
# 
# # cluster analysis and dendrograms
# library(MASS)
# library(shapes)
# 
# # interpolate point so each shape has the same number
# interp <- Coo(lapply(outlines, function(i) coo.sample.int(simplify2array(i), 3000) ))
# 
# gorf<-gorf.dat
# panf<-panf.dat[c(5,1,2:4,6:8),,]
# pongof<-pongof.dat[c(5,1,2:4,6:8),,]
# APE<-array(c(panf, gorf, pongof),dim=c(8,2,80))
# APE<-array(c(panf, gorf),dim=c(8,2,56))
# AP<-orp(pgpa(APE)$rotated)
# m<-t(matrix(AP, 16, 56))
# par(mar=c(0.5,2,1,1))
# layout(matrix(c(1,2),2,1))
# plot(hclust(dist(m), method="average"),main="UPGMA"
#      ,labels=c(rep("P",26),rep("G",30)),cex=0.7)
# plot(hclust(dist(m), method="complete"),main="COMPLETE"
#      ,labels=c(rep("P",26),rep("G",30)),cex=0.7)
# 
# 
# 
# # from package vignette...
# botFg <- meanShapes(botF)
# str(botFg)
# stack(Coo(botFg), borders=col.gallus(2))
# tps.grid(botFg[[1]], botFg[[2]]) # have a look to ?"["
# tps.grid(botFg$beer, botFg$whisky, grid.outside=2,
#          grid.size=80, amp=3, plot.full=FALSE)
# 
# set.seed(007)
# my_list<- lapply(1:10, function(i) matrix(data = runif(sample(seq(2,10,2), 1)), ncol = 2))
# 
# my_array<-Map(function(x,y,z) array(n,c(x,y,z)), myrow, mycol,myz) # but this creates the list of 1 array
