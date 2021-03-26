library(scales)

#import data
coords = readRDS(file = "data/magdif/max_pca_num.rds") 
colors = readRDS(file = "data/magdif/colors.rds")
mb_abs = readRDS(file = "data/magdif/abs_magdiff.rds") 
mb_rel = readRDS(file = "data/magdif/rel_magdiff.rds")
arrows = readRDS(file = "data/aging_vectors/arrows.rds")

young_r = readRDS(file = "data/aging_vectors/new_coords/y_coords_r.rds")[["arrows4"]]
old_r = readRDS(file = "data/aging_vectors/new_coords/o_coords_r.rds")[["arrows4"]]

y_nms = grep("young",rownames(arrows),value = TRUE)
o_nms = grep("old",rownames(arrows),value = TRUE)

#rownames(young_g) = y_nms
rownames(young_r) = y_nms
#rownames(old_g) = o_nms
rownames(old_r) = o_nms

arrows = rbind(young_r,old_r)
arrows = arrows[rownames(coords),]


c1 = names(colors)[colors=="#4882C3"]
c2 = names(colors)[colors=="#F26A6A"]
c3 = names(colors)[colors=="#13751B"]
c4 = names(colors)[colors=="#FF6A00"]
c5 = names(colors)[colors=="#E2CF00"]

colors[c(c1,c2,c3)] = "#4882C3"
c1 = c(c1,c2,c3)


#draw magdif 
#par(mar = c(3,3,3,3),mfrow = c(3,3))
cex = .7

exp = mb_abs
#exp_sign = ((exp >=0)*2)-1 
#exp = log2(abs(exp * 50)) * exp_sign

length = 100
val_range = seq(-max(abs(exp[grep("old",rownames(coords))])),
                max(abs(exp[grep("old",rownames(coords))])),
                length.out = length)
marker = colorRampPalette(c("#0060C3","white","#ff4d4d"))(100)

clrs = marker[sapply(exp,function(v){which.min(abs(val_range-v))})]

###################
###  FIGURE 3B  ###
###################
plot(coords, col = colors, cex = cex, pch = 19,xaxt='n',yaxt ="n",main = "Principal Component Analysis",xlab = "",ylab = "")
points(coords[grep("old",rownames(coords)),],pch = 5, cex = 1.3)
title(main = "",xlab = "PC1 (71.87%)", ylab="PC2 (6.99%)", line=1, cex.lab=1.2)



#plot all together
cex = .7
y_cls = grep("young",rownames(arrows),value = TRUE)
o_cls = grep("old",rownames(arrows),value = TRUE)

###################
###  FIGURE 3D  ###
###################
plot(coords, col = "white",xaxt='n',yaxt ="n", main = "",xlab = "",ylab = "")
points(coords[y_cls,], col = colors[y_cls], cex = cex, pch = 19)
points(coords[o_cls,], col = alpha(colors[o_cls],.3), cex = cex, pch = 19)
arrows(arrows[y_cls,1],
       arrows[y_cls,2],
       arrows[y_cls,3],
       arrows[y_cls,4],
       length = .05, angle = 30, code = 2,lwd = 1, col = "black")
title(main = "2 month old",xlab = "PC1", ylab="PC2", line=1, cex.lab=1.2)


###################
###  FIGURE 3E  ###
###################
plot(coords, col = "white",xaxt='n',yaxt ="n",main = "",xlab = "",ylab = "")
points(coords[o_cls,], col = colors[o_cls], cex = cex, pch = 19)
points(coords[y_cls,], col = alpha(colors[y_cls],.3), cex = cex, pch = 19)
arrows(arrows[o_cls,1],
       arrows[o_cls,2],
       arrows[o_cls,3],
       arrows[o_cls,4],
       length = .05, angle = 30, code = 2,lwd = 1, col = "black")
title(main = "4.5 month old",xlab = "PC1", ylab="PC2", line=1, cex.lab=1.2)



###################
###  FIGURE 3F  ###
###################
lwd = 4
cex = .3
scale = 5

og = list(c1= c1,c2 = c2,c3 = c3, c4=c4,c5=c5)
new = list(c1= c1,c4=c4,c5=c5)

grp_lst = new

group_origins = grp_avg_origin(arrows[,1:2],grp_lst)

all_cells_transposed = get_transposed_vectors(arrows,grp_lst,scale = scale)

group_avg_cells_transposed = get_average_group_vector(arrows,
                                                      group_lst = grp_lst,
                                                      comparison_lst = list(young = y_cls,old = o_cls),
                                                      scale = scale)


#individual cells
plot(coords, col = colors, cex = cex, pch = 19,xaxt='n',yaxt ="n", main = "State changes")
points(group_origins)
arrows(all_cells_transposed[y_cls,1],
       all_cells_transposed[y_cls,2],
       all_cells_transposed[y_cls,3],
       all_cells_transposed[y_cls,4],
       length = .05, angle = 30, code = 2,lwd = 1, col = alpha("#0060C3",.3))
arrows(all_cells_transposed[o_cls,1],
       all_cells_transposed[o_cls,2],
       all_cells_transposed[o_cls,3],
       all_cells_transposed[o_cls,4],
       length = .05, angle = 30, code = 2,lwd = 1, col = alpha("#ff4d4d",.3))

arrows(group_avg_cells_transposed$young[,1],
       group_avg_cells_transposed$young[,2],
       group_avg_cells_transposed$young[,3],
       group_avg_cells_transposed$young[,4],
       length = .1, angle = 35, code = 2,lwd = lwd, col = "#0060C3")

arrows(group_avg_cells_transposed$old[,1],
       group_avg_cells_transposed$old[,2],
       group_avg_cells_transposed$old[,3],
       group_avg_cells_transposed$old[,4],
       length = .1, angle = 35, code = 2,lwd = lwd, col = "#ff4d4d")

