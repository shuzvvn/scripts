#!/usr/bin/env Rscript

# combine_heatmaps.4.R

# Shu-Ting Cho <vivianlily6@hotmail.com>
# combine four heatmaps with same row and col order
# invert color if needed 
# v1 2021/01/13
# v2 2021/01/15 show value in cell
#               use percent for Roary
# v3 2021/01/21 auto calculate output size by input matrix size or set w/h
#               allow not input all four mt
#               heatmap mt reorder by the first mt's dend
# v4 2021/02/09 input as mt1, mt2, mt3, mt4, then assign names later
#               heatmap mt reorder by the MT1's dend
#               convert decimal to percent for matrix in a given list
#               invert color scale for matrix in a given list
#               assign max of color scale to specific matrix

# Usage: 
# combine_heatmaps.4.R \
# --MT1=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_SNPs.csv \
# --MT2=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_cgMLST.csv \
# --MT3=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_SKA_SNPs.csv \
# --MT4=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_MLVA.csv \
# --names='SNPs,cgMLST,SKA,MLVA' \
# --invert='MT2' --percent='MT4' --max1=10 --max4=100 \
# --outfile=/mnt/c/Users/vivia/project/edshat02/edshat02.03/CD_heatmaps.pdf


################### Require Packages ####################
# 1. optparse                                           #
# install.packages("optparse")                          #
#                                                       #
# 2. ComplexHeatmap                                     #
# if (!requireNamespace("BiocManager", quietly = TRUE)) #
#     install.packages("BiocManager")                   #
# BiocManager::install("ComplexHeatmap")                #
#                                                       #
# 3. scales                                             #
# install.packages("scales")                            #
#                                                       #
# 4. circlize                                           #
#########################################################

# load libraries
library(optparse)

######################## get options ########################
# specify options in a list
option_list = list(
	make_option(c("-a", "--MT1"), type="character", default="NA", 
		help="MT1 file name (Required)", metavar="filename"),
	make_option(c("-b", "--MT2"), type="character", default="NA", 
		help="MT2 file name (Optional)", metavar="filename"),
	make_option(c("-c", "--MT3"), type="character", default="NA", 
		help="MT3 file name (Optional)", metavar="filename"),
	make_option(c("-d", "--MT4"), type="character", default="NA", 
		help="MT4 file name (Optional)", metavar="filename"),
    make_option(c("-n", "--names"), type="character", default="MT_NAME_1,MT_NAME_2,MT_NAME_3,MT_NAME_4", 
		help="heatmap names [default= %default]", metavar="data1,data2,data3,data4"),
    make_option(c("-o", "--outfile"), type="character", default="out_heatmaps.pdf", 
    	help="output.pdf file name [default= %default]", metavar="filename"),

    make_option(c("-e", "--max1"), type="double", default=0, 
		help="Maximum number for color scale of MT1 [default= auto]", metavar="number"),
    make_option(c("-f", "--max2"), type="double", default=0, 
		help="Maximum number for color scale of MT2 [default= auto]", metavar="number"),
    make_option(c("-g", "--max3"), type="double", default=0, 
		help="Maximum number for color scale of MT3 [default= auto]", metavar="number"),
    make_option(c("-i", "--max4"), type="double", default=0, 
		help="Maximum number for color scale of MT4 [default= auto]", metavar="number"),

    make_option(c("-r", "--invert"), type="character", default="", 
		help="list of data which color should be inverted", metavar="MT1,MT3"),
    make_option(c("-p", "--percent"), type="character", default="", 
		help="list of data which value should be present as percentage number", metavar="MT1,MT3"),

    make_option(c("-x", "--col_low"), type="character", default="yellow", 
		help="heatmap color low [default= %default]", metavar="color"),
    make_option(c("-y", "--col_high"), type="character", default="blue", 
		help="heatmap color high [default= %default]", metavar="color"),
    make_option(c("-w", "--width"), type="double", default=0, 
		help="output file width [default= auto]", metavar="number"),
    make_option(c("-t", "--height"), type="double", default=0, 
		help="output file height [default= auto]", metavar="number"),
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
    	help="Print extra output [default= %default]")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--MT1=FILENAME_1] [--MT2=FILENAME_2] [--MT3=FILENAME_3] [--MT4=FILENAME_4] [--names=MT_NAME_1,MT_NAME_2,MT_NAME_3,MT_NAME_4] [--outfile=FILENAME]", option_list=option_list)) 




######################## read inputs ########################

### heatmap names, default: MT_NAME_1, MT_NAME_2, MT_NAME_3, MT_NAME_4
# if item in names not in input, remove from names
heatmap.names <- strsplit(opt$names, ",")[[1]]
# heatmap.names <- heatmap.names[heatmap.names %in% names(matrix.list)]

if ( opt$verbose ) {
	message <- cat("heatmap names:", heatmap.names)
	write(message, stderr())
}


### list of data which should be shown as percentage number
data.perc <- strsplit(opt$percent, ",")[[1]]
if ( opt$verbose ) {
	message <- cat("shown as percentage:", data.perc)
	write(message, stderr())
}


### list of data which color should be inverted
data.rev <- strsplit(opt$invert, ",")[[1]]
if ( opt$verbose ) {
	message <- cat("invert color:", data.rev)
	write(message, stderr())
}







### Reading input matrix
if ( opt$verbose ) { 
    write("Reading input matrix...", stderr())
}

library(scales)

matrix.list <- list()

# read tables
# data frame to matrix
# combine matrix into one list 
if ( opt$MT1 != "NA" ) {
	MT1.df <- read.csv(opt$MT1, sep=",", header = T, row.names = 1)
	MT1.mt <- as.matrix(MT1.df)
	if (is.element("MT1", data.perc)) {
		MT1.mt[] <- as.numeric(percent(MT1.mt, accuracy=.1, suffix = ""))		
	} 
	matrix.list[["MT1"]] <- MT1.mt
} 
if ( opt$MT2 != "NA" ) {
	MT2.df <- read.csv(opt$MT2, sep=",", header = T, row.names = 1)
	MT2.mt <- as.matrix(MT2.df)
	if (is.element("MT2", data.perc)) {
		MT2.mt[] <- as.numeric(percent(MT2.mt, accuracy=.1, suffix = ""))		
	} 
	matrix.list[["MT2"]] <- MT2.mt
} 
if ( opt$MT3 != "NA" ) {
	MT3.df <- read.csv(opt$MT3, sep=",", header = T, row.names = 1)
	MT3.mt <- as.matrix(MT3.df)
	if (is.element("MT3", data.perc)) {
		MT3.mt[] <- as.numeric(percent(MT3.mt, accuracy=.1, suffix = ""))		
	} 
	matrix.list[["MT3"]] <- MT3.mt
} 
if ( opt$MT4 != "NA" ) {
	MT4.df <- read.csv(opt$MT4, sep=",", header = T, row.names = 1)
	MT4.mt <- as.matrix(MT4.df)
	if (is.element("MT4", data.perc)) {
		MT4.mt[] <- as.numeric(percent(MT4.mt, accuracy=.1, suffix = ""))		
	} 
	matrix.list[["MT4"]] <- MT4.mt
}


### output pdf file
outfile <- opt$outfile




# # heatmap ordered by heatmap.names
# matrix.list <- matrix.list[heatmap.names]

if ( opt$verbose ) {
	message <- cat("input matrix size:", nrow(matrix.list[[1]]), "x", ncol(matrix.list[[1]]))
	write(message, stderr())
}


### heatmap color, default: low -> high : yellow -> blue
col_low <- opt$col_low
col_high <- opt$col_high
if ( opt$verbose ) {
	message <- cat("color:", col_low, "->", col_high)
	write(message, stderr())
}


### output size, default: 0 (auto)
width = opt$width
height = opt$height


### max for each color scale, default: 0 (auto)
for (i in c(1:4))
{
	var <- paste("max", i, sep = "")
	assign(var, opt[[var]])
}


######################## main ########################
library(ComplexHeatmap)
library(circlize)

# row and col ordered by the dendrograms of the first heatmap
out <- Heatmap(matrix.list[[1]])
gene.order <- rownames(matrix.list[[1]])[row_order(out)]


# reorder row and col in matrix.list
for ( mt in names(matrix.list) )
{
	matrix.list[[mt]] <- matrix.list[[mt]][gene.order, gene.order]
}




# get heatmaps
i = 0
for (data in names(matrix.list))
{
	i = i + 1
	htp <- paste("ht", i, sep = "")

	# check if color should be perc
	if (is.element(data, data.perc)) {
		fontsize = 7
	} else {
		fontsize = 10
	}

	# check if color should be invert
	if (is.element(data, data.rev)) {
		cols <- c(col_high, col_low)
		text.col.h <- "black"
		text.col.l <- "white"
	} else {
		cols <- c(col_low, col_high)
		text.col.h <- "white"
		text.col.l <- "black"
	}

	# draw heatmap, because cell_fun is weird, have to generate new data everytime
	if (data == "MT1") {
		mat1 <- matrix.list[[data]]
		fontsize1 <- fontsize
		text.col.h1 <- text.col.h
		text.col.l1 <- text.col.l
		# check if need to change max for color scale
		if (max1 != 0) {
			col_fun = colorRamp2(c(0, max1), cols)
		} else {
			col_fun = colorRamp2(c(0, max(mat1)), cols)
			max1 <- max(mat1)
		}
		# draw heatmap
		assign(htp, 
			Heatmap(mat1, name = heatmap.names[1], col = col_fun,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat1[i, j] < max1/2) {
        			grid.text(sprintf("%s", mat1[i, j]), x, y, gp = gpar(fontsize = fontsize1, col = text.col.l1))
        		} else {
        			grid.text(sprintf("%s", mat1[i, j]), x, y, gp = gpar(fontsize = fontsize1, col = text.col.h1))
        		}
   			 }))
	} else if (data == "MT2") {
		mat2 <- matrix.list[[data]]
		fontsize2 <- fontsize
		text.col.h2 <- text.col.h
		text.col.l2 <- text.col.l
		# check if need to change max for color scale
		if (max2 != 0) {
			col_fun = colorRamp2(c(0, max2), cols)
		} else {
			col_fun = colorRamp2(c(0, max(mat2)), cols)
			max2 <- max(mat2)
		}
		assign(htp, 
			Heatmap(mat2, name = heatmap.names[2], col = col_fun,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat2[i, j] < max2/2) {
        			grid.text(sprintf("%s", mat2[i, j]), x, y, gp = gpar(fontsize = fontsize2, col = text.col.l2))
        		} else {
        			grid.text(sprintf("%s", mat2[i, j]), x, y, gp = gpar(fontsize = fontsize2, col = text.col.h2))
        		}
   			 }))
	} else if (data == "MT3") {
		mat3 <- matrix.list[[data]]
		fontsize3 <- fontsize
		text.col.h3 <- text.col.h
		text.col.l3 <- text.col.l
		# check if need to change max for color scale
		if (max3 != 0) {
			col_fun = colorRamp2(c(0, max3), cols)
		} else {
			col_fun = colorRamp2(c(0, max(mat3)), cols)
			max3 <- max(mat3)
		}
		assign(htp, 
			Heatmap(mat3, name = heatmap.names[3], col = col_fun,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat3[i, j] < max3/2) {
        			grid.text(sprintf("%s", mat3[i, j]), x, y, gp = gpar(fontsize = fontsize3, col = text.col.l3))
        		} else {
        			grid.text(sprintf("%s", mat3[i, j]), x, y, gp = gpar(fontsize = fontsize3, col = text.col.h3))
        		}
   			 }))
	} else {
		mat4 <- matrix.list[[data]]
		fontsize4 <- fontsize
		text.col.h4 <- text.col.h
		text.col.l4 <- text.col.l
		# check if need to change max for color scale
		if (max4 != 0) {
			col_fun = colorRamp2(c(0, max4), cols)
		} else {
			col_fun = colorRamp2(c(0, max(mat4)), cols)
			max4 <- max(mat4)
		}
		assign(htp, 
			Heatmap(mat4, name = heatmap.names[4], col = col_fun,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat4[i, j] < max4/2) {
        			grid.text(sprintf("%s", mat4[i, j]), x, y, gp = gpar(fontsize = fontsize4, col = text.col.l4))
        		} else {
        			grid.text(sprintf("%s", mat4[i, j]), x, y, gp = gpar(fontsize = fontsize4, col = text.col.h4))
        		}
   			 }))
	}
}



# calculate output w and h, write output pdf

if (length(matrix.list) == 4) {
	if ( width == 0 && height == 0 ) {
		height = nrow(matrix.list[[1]]) * 0.307 + 1.3632
		width = height * 4 - 5
	} else if ( width == 0 && height != 0 ) {
		width = height * 4 - 5
	} else if ( width != 0 && height == 0 ) {
		height = (width + 5) / 4
	}

	if ( opt$verbose ) {
		message <- cat("output size (W x H):", width, "x", height, "\nwriting output...")
		write(message, stderr())
	}

	pdf(file = outfile, width=width, height=height)
	# combine four heatmaps
	draw(ht1 + ht2 + ht3 + ht4, heatmap_legend_side = "bottom")
	dev.off()

} else if (length(matrix.list) == 3) {
	if ( width == 0 && height == 0 ) {
		height = nrow(matrix.list[[1]]) * 0.307 + 1.3632
		width = height * 3 - 3.3
	} else if ( width == 0 && height != 0 ) {
		width = height * 3 - 3.3
	} else if ( width != 0 && height == 0 ) {
		height = (width + 3.3) / 3
	}

	if ( opt$verbose ) {
		message <- cat("output size (W x H):", width, "x", height, "\nwriting output...")
		write(message, stderr())
	}

	pdf(file = outfile, width=width, height=height)
	# combine four heatmaps
	draw(ht1 + ht2 + ht3, heatmap_legend_side = "bottom")
	dev.off() 

} else if (length(matrix.list) == 2) {
	if ( width == 0 && height == 0 ) {
		height = nrow(matrix.list[[1]]) * 0.307 + 1.3632
		width = height * 1.9896 - 1.9107
	} else if ( width == 0 && height != 0 ) {
		width = height * 1.9896 - 1.9107
	} else if ( width != 0 && height == 0 ) {
		height = (width + 1.9107) / 1.9896
	}

	if ( opt$verbose ) {
		message <- cat("output size (W x H):", width, "x", height, "\nwriting output...")
		write(message, stderr())
	}

	pdf(file = outfile, width=width, height=height)
# combine four heatmaps
draw(ht1 + ht2, heatmap_legend_side = "bottom")
dev.off() 
} else {
	if ( width == 0 && height == 0 ) {
		height = nrow(matrix.list[[1]]) * 0.307 + 1.3632
		width = height * 0.9977 - 0.5045
	} else if ( width == 0 && height != 0 ) {
		width = height * 0.9977 - 0.5045
	} else if ( width != 0 && height == 0 ) {
		height = (width + 0.5045) / 0.9977
	}

	if ( opt$verbose ) {
		message <- cat("output size (W x H):", width, "x", height, "\nwriting output...")
		write(message, stderr())
	}

	pdf(file = outfile, width=width, height=height)
# combine four heatmaps
draw(ht1, heatmap_legend_side = "bottom")
dev.off() 
}
