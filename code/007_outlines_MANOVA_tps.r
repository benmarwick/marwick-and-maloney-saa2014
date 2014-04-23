# Statistical test differences between groups 
# through a Multivariate Analysis of Variance.
print(mnv <- manova.Coe(coords_F, "class", retain = 5))

#  hierarchical clustering that hinges on dist and
# hclust for calculation, and phylo.plot from the 
# ape package for graphical output
# clust(coords_F) # default for momocs but misses names

cl <- hclust(d = dist(coords_F@coe, method = 'euclidean'))
plot(cl, labels = coords_F@names)


# inspect thin plate splines
x <- meanShapes(coords_F, 'class')
str(x) # what we got
stack(Coo(x), borders=col.gallus(2))
# compare a few classes
# tps.grid(x[[3]],  x[[4]]) # have a look to ?"["
tps.grid(x$flake, x$retouched, grid.outside=2,
         grid.size=80, amp=1, plot.full=FALSE)
