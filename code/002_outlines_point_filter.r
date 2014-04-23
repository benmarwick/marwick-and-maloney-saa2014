# 
# check to see how many coords the shapes have, we want to exclude 
# any odd ones with very few
n <- 300
# get index of shape with less than n coords
idx <- lapply(outlines, function(i) length(i[[1]])) < n
mtch <- match(TRUE, idx) 
nm <- ifelse(is.na(match(TRUE, idx)), 0, match(TRUE, idx) )  # which one(s) have less than n points?
# subset to get only those with >500 points
outlines_n <- if(nm == 0)  outlines else outlines[-nm]
# update other data
fac_d <- fac_d[!idx, ]

# report
message(paste0("Started with ", length(idx), " outlines, and removed ", ifelse(nm == 0, 0, length(nm)), " leaving ", length(outlines_n), " outlines for further analysis." ))
