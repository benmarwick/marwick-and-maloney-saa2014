# 
# check to see how many coords the shapes have, we want to exclude 
# any odd ones with very few
# get index of shape with less than 500 coords
idx <- lapply(outlines, function(i) length(i[[1]])) < 500
nm <- match(TRUE, idx) # which one(s) have less than 500 points?
# subset to get only those with >500 points
outlines_500 <- outlines[-nm]
