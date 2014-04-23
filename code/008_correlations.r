# correlation of PC axis with retouch index

library(gdata)

dat <- read.xls("C:\\Users\\marwick\\UW-Dropbox\\Dropbox\\Malonys-lithic-scan-data\\spreadsheet ME3.xls")
str(dat)

retouch <- na.omit(dat[dat$II > 0, ])

str(coords_D)
# first PC
coords_D$li[,1]
# second PC
# first PC
coords_D$li[,2]
# check
plot(coords_D$li[,1], coords_D$li[,2])

PCs <- data.frame(name = row.names(coords_D$coe), 
           PC1 = coords_D$li[,1], 
           PC2 = coords_D$li[,2])

retouch_subset <- retouch[retouch$Image_Name %in% intersect(retouch$Image_Name, PCs$name), ]
PCs_subset <- PCs[PCs$name%in% intersect(retouch$Image_Name, PCs$name), ]
retouch_subset[,5:ncol(retouch_subset)] <- sapply(retouch_subset[,5:ncol(retouch_subset)], as.numeric, as.character)
cors <- lapply(5:ncol(retouch_subset), function(i) cor.test(retouch_subset[,i], PCs_subset$PC1 ))
names(cors) <- names(retouch_subset[,5:ncol(retouch_subset)] )
# inspect output
cors

# plot
join <- merge(retouch_subset, PCs_subset, by.x = "Image_Name", by.y = "name")
library(ggplot2)
ggplot(join, aes(PC1, GIUR)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm") +
  theme_minimal(base_size = 16)

# Kimberley point correlation of PC1 with serration

dat <- read.xls("C:\\Users\\marwick\\UW-Dropbox\\Dropbox\\Malonys-lithic-scan-data\\spreadsheet_allsites_EFA.xls")
kimb <- dat[grepl("kimberleypoint", dat$Image_Name), ]

str(coords_D)
# first PC
coords_D$li[,1]
# second PC
# first PC
coords_D$li[,2]
# check
plot(coords_D$li[,1], coords_D$li[,2])

PCs <- data.frame(name = row.names(coords_D$coe), 
                  PC1 = coords_D$li[,1], 
                  PC2 = coords_D$li[,2])

kimb_subset <- kimb[kimb$Image_Name %in% intersect(kimb$Image_Name, PCs$name), ]
PCs_subset <- PCs[PCs$name%in% intersect(kimb$Image_Name, PCs$name), ]
kimb_subset[,5:ncol(kimb_subset)] <- sapply(kimb_subset[,5:ncol(kimb_subset)], as.numeric, as.character)
cor.test(kimb_subset$Medial_W_th, PCs_subset$PC1, use = "complete.obs")
plot(kimb_subset$Medial_W_th, PCs_subset$PC1)
cor.test(kimb_subset$n_serr_per_10mm, PCs_subset$PC2, use = "complete.obs")

perct_perimter_serr
n_serr_per_10mm
names(cors) <- names(kimb_subset[,5:ncol(kimb_subset)] )
# inspect output
cors

# plot
join <- merge(kimb_subset, PCs_subset, by.x = "Image_Name", by.y = "name")
library(ggplot2)
ggplot(join, aes(PC1, n_serr)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm") +
  theme_minimal(base_size = 16)
