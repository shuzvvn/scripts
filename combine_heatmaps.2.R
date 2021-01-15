#!/usr/bin/env Rscript

# combine_heatmaps.2.R

# Shu-Ting Cho <vivianlily6@hotmail.com>
# combine four heatmaps with same row and col order
# invert color if needed 
# v1 2021/01/13
# v2 2021/01/15 show value in cell, use percent for Roary

# Usage: 
# combine_heatmaps.2.R \
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
library(ComplexHeatmap)
library(scales)

######################## get options ########################
# specify options in a list
option_list = list(
	make_option(c("-a", "--SNPs"), type="character", default=NULL, 
		help="SNPs.csv file name (Required)", metavar="filename"),
	make_option(c("-b", "--cgMLST"), type="character", default=NULL, 
		help="cgMLST.csv file name (Required)", metavar="filename"),
	make_option(c("-c", "--SKA"), type="character", default=NULL, 
		help="SKA.csv file name (Required)", metavar="filename"),
	make_option(c("-d", "--Roary"), type="character", default=NULL, 
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
    make_option(c("-w", "--width"), type="integer", default=20, 
		help="output file width [default= %default]", metavar="number"),
    make_option(c("-t", "--height"), type="integer", default=6, 
		help="output file height [default= %default]", metavar="number"),
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
    	help="Print extra output [default]")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--SNPs=FILENAME] [--cgMLST=FILENAME] [--SKA=FILENAME] [--Roary=FILENAME] [--outfile=FILENAME]", option_list=option_list)) 


######################## read inputs ########################
# print some progress messages to stderr if "quietly" wasn't requested
if ( opt$verbose ) { 
    write("Reading input matrix...", stderr())
}

# read tables
SNPs.df <- read.csv(opt$SNPs, sep=",", header = T, row.names = 1)
cgMLST.df <- read.csv(opt$cgMLST, sep=",", header = T, row.names = 1)
SKA.df <- read.csv(opt$SKA, sep=",", header = T, row.names = 1)
Roary.df <- read.csv(opt$Roary, sep=",", header = T, row.names = 1)

# output pdf file
outfile <- opt$outfile


# heatmap order, default: SNPs, cgMLST, SKA, Roary
heatmap.order <- strsplit(opt$order, ",")[[1]]
if ( opt$verbose ) {
	message <- cat("heatmap order:", heatmap.order)
	write(message, stderr())
}

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

# output size, default: 10 x 4
width = opt$width
height = opt$height
if ( opt$verbose ) {
	message <- cat("output size:", width, "x", height)
	write(message, stderr())
}



######################## main ########################
# data frame to matrix
SNPs.mt <- as.matrix(SNPs.df)
cgMLST.mt <- as.matrix(cgMLST.df)
SKA.mt <- as.matrix(SKA.df)
Roary.mt <- as.matrix(Roary.df)

# Roary is decimal, convert to percent but not show "%"
# Roary.mt.perc <- Roary.mt
Roary.mt[] <- as.numeric(percent(Roary.mt, accuracy=.1, suffix = ""))


# row and col ordered by SNPs dendrograms 
out <- Heatmap(SNPs.mt)
gene.order <- rownames(SNPs.mt)[row_order(out)]

# combine matrix into one list and reorder row and col
matrix.list <- list(SNPs = SNPs.mt[gene.order, gene.order], 
	                cgMLST = cgMLST.mt[gene.order, gene.order], 
	                SKA = SKA.mt[gene.order, gene.order], 
	                Roary = Roary.mt[gene.order, gene.order])

#Roary.mt.perc <- Roary.mt.perc[gene.order, gene.order]

# heatmap ordered by heatmap.order
matrix.list <- matrix.list[heatmap.order]


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



# write output pdf
pdf(file = outfile, width=width, height=height)
# combine four heatmaps
draw(ht1 + ht2 + ht3 + ht4, heatmap_legend_side = "bottom")
dev.off() 