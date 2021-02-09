#!/usr/bin/env Rscript

# combine_heatmaps.3.R

# Shu-Ting Cho <vivianlily6@hotmail.com>
# combine four heatmaps with same row and col order
# invert color if needed 
# v1 2021/01/13
# v2 2021/01/15 show value in cell
#               use percent for Roary
# v3 2021/01/21 auto calculate output size by input matrix size or set w/h
#               allow not input all four mt
#               heatmap mt reorder by the first mt's dend

# Usage: 
# combine_heatmaps.3.R \
# --SNPs=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_SNPs.csv \
# --cgMLST=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_cgMLST.csv \
# --SKA=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_SKA_SNPs.csv \
# --Roary=/mnt/c/Users/vivia/project/edshat02/rawData/CD_ST1_A_roary_corr.csv \
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
#########################################################

# load libraries
library(optparse)

######################## get options ########################
# specify options in a list
option_list = list(
	make_option(c("-a", "--SNPs"), type="character", default="NA", 
		help="SNPs.csv file name (Required)", metavar="filename"),
	make_option(c("-b", "--cgMLST"), type="character", default="NA", 
		help="cgMLST.csv file name (Required)", metavar="filename"),
	make_option(c("-c", "--SKA"), type="character", default="NA", 
		help="SKA.csv file name (Required)", metavar="filename"),
	make_option(c("-d", "--Roary"), type="character", default="NA", 
		help="Roary.csv file name (Required)", metavar="filename"),
    make_option(c("-o", "--outfile"), type="character", default="out_heatmaps.pdf", 
    	help="output.pdf file name [default= %default]", metavar="filename"),
    make_option(c("-l", "--order"), type="character", default="SNPs,cgMLST,SKA,Roary", 
		help="heatmap order [default= %default]", metavar="data1,data2,data3,data4"),
    make_option(c("-r", "--invert"), type="character", default="", 
		help="list of data which color should be inverted", metavar="data1,data3"),
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
opt <- parse_args(OptionParser(usage = "%prog [options] [--SNPs=FILENAME] [--cgMLST=FILENAME] [--SKA=FILENAME] [--Roary=FILENAME] [--outfile=FILENAME]", option_list=option_list)) 


######################## read inputs ########################
# print some progress messages to stderr if "quietly" wasn't requested
if ( opt$verbose ) { 
    write("Reading input matrix...", stderr())
}

matrix.list <- list()

# read tables
# data frame to matrix
# combine matrix into one list 
if ( opt$SNPs != "NA" ) {
	SNPs.df <- read.csv(opt$SNPs, sep=",", header = T, row.names = 1)
	SNPs.mt <- as.matrix(SNPs.df)
	matrix.list[["SNPs"]] <- SNPs.mt
} 
if ( opt$cgMLST != "NA" ) {
	cgMLST.df <- read.csv(opt$cgMLST, sep=",", header = T, row.names = 1)
	cgMLST.mt <- as.matrix(cgMLST.df)
	matrix.list[["cgMLST"]] <- cgMLST.mt
} 
if ( opt$SKA != "NA" ) {
	SKA.df <- read.csv(opt$SKA, sep=",", header = T, row.names = 1)
	SKA.mt <- as.matrix(SKA.df)
	matrix.list[["SKA"]] <- SKA.mt
} 
if ( opt$Roary != "NA" ) {
	Roary.df <- read.csv(opt$Roary, sep=",", header = T, row.names = 1)
	Roary.mt <- as.matrix(Roary.df)
	# Roary is decimal, convert to percent but not show "%"
	library(scales)
	Roary.mt[] <- as.numeric(percent(Roary.mt, accuracy=.1, suffix = ""))
	matrix.list[["Roary"]] <- Roary.mt
}


# output pdf file
outfile <- opt$outfile


# heatmap order, default: SNPs, cgMLST, SKA, Roary
# if item in order not in input, remove from order
heatmap.order <- strsplit(opt$order, ",")[[1]]
heatmap.order <- heatmap.order[heatmap.order %in% names(matrix.list)]

if ( opt$verbose ) {
	message <- cat("heatmap order:", heatmap.order)
	write(message, stderr())
}

# heatmap ordered by heatmap.order
matrix.list <- matrix.list[heatmap.order]


message <- cat("input matrix size:", nrow(matrix.list[[1]]), "x", ncol(matrix.list[[1]]))
write(message, stderr())

# list of data which color should be inverted
data.rev <- strsplit(opt$invert, ",")[[1]]
if ( opt$verbose ) {
	message <- cat("invert color:", data.rev)
	write(message, stderr())
}

# heatmap color, default: low -> high : yellow -> blue
col_low <- opt$col_low
col_high <- opt$col_high
if ( opt$verbose ) {
	message <- cat("color:", col_low, "->", col_high)
	write(message, stderr())
}

# output size, default: 0 (auto)
width = opt$width
height = opt$height

######################## main ########################
library(ComplexHeatmap)


# row and col ordered by the dendrograms of the first heatmap
out <- Heatmap(matrix.list[[1]])
gene.order <- rownames(matrix.list[[1]])[row_order(out)]


# reorder row and col in matrix.list
for ( mt in names(matrix.list) )
{
	matrix.list[[mt]] <- matrix.list[[mt]][gene.order, gene.order]
}


# color scale, normal or invert
col.f <- colorRampPalette(c(col_low, col_high))(256)
col.r <- rev(col.f)


# get heatmaps
i = 0
for (data in names(matrix.list))
{
	i = i + 1
	htp <- paste("ht", i, sep = "")
	# check if color should be invert
	if (is.element(data, data.rev == TRUE)) {
		col <- col.r
	} else {
		col <- col.f
	}
	# draw heatmap, because cell_fun is weird, have to generate new data everytime
	if (data == "Roary") {
		mat <- matrix.list[[data]]
		assign(htp, 
			Heatmap(mat, name = data, col = col,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat[i, j] < max(mat)/2) {
        			grid.text(sprintf("%s", mat[i, j]), x, y, gp = gpar(fontsize = 7))
        		} else {
        			grid.text(sprintf("%s", mat[i, j]), x, y, gp = gpar(fontsize = 7, col = "white"))
        		}
   			 }))
	} else if (data == "SNPs") {
		mat2 <- matrix.list[[data]]
		assign(htp, 
			Heatmap(mat2, name = data, col = col,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat2[i, j] < max(mat2)/2) {
        			grid.text(sprintf("%i", mat2[i, j]), x, y, gp = gpar(fontsize = 10))
        		} else {
        			grid.text(sprintf("%i", mat2[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
        		}
   			 }))
	} else if (data == "SKA") {
		mat3 <- matrix.list[[data]]
		assign(htp, 
			Heatmap(mat3, name = data, col = col,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat3[i, j] < max(mat3)/2) {
        			grid.text(sprintf("%i", mat3[i, j]), x, y, gp = gpar(fontsize = 10))
        		} else {
        			grid.text(sprintf("%i", mat3[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
        		}
   			 }))
	} else {
		mat4 <- matrix.list[[data]]
		assign(htp, 
			Heatmap(mat4, name = data, col = col,
				cluster_rows = FALSE, cluster_columns = FALSE, # no dendrogram
        		row_names_side = "left", column_names_side = "top",
        		heatmap_legend_param = list(direction = "horizontal"),
        		cell_fun = function(j, i, x, y, width, height, fill) {
        		if (mat4[i, j] < max(mat4)/2) {
        			grid.text(sprintf("%i", mat4[i, j]), x, y, gp = gpar(fontsize = 10))
        		} else {
        			grid.text(sprintf("%i", mat4[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
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
