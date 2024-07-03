rm(list = ls())  
if(!require("gwasrapidd")) install.packages("gwasrapidd", update = FALSE, ask = FALSE)
if(!require("stringr")) install.packages("stringr", update = FALSE, ask = FALSE)
if(!require("data.table")) install.packages("data.table", update = FALSE, ask = FALSE)
if (!require("devtools")) { install.packages("devtools") } else {}
# devtools::install_github("rondolab/MR-PRESSO", force = TRUE)
library(TwoSampleMR)
library(MRPRESSO)
library(circlize)
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)

# Load data
result_all <- read.csv("C:/Users/gary/Desktop/2.1400 to PD/2.1400 to PD/1400_PD_res.csv")
result_all2 <- read.csv("C:/Users/gary/Desktop/1.338 to PD/1.338 to PD/338_PD_res.csv")

# Function to process the data
process_data <- function(data, sample_type) {
  methods <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Simple mode", "Weighted mode")
  col_names <- c("Metabolites", "Beta", "P-Value", "Odd Ratio")
  result_list <- list()
  
  for (method in methods) {
    subset <- data[data$method == method,][, c(4, 7, 9, if (method == "Inverse variance weighted") 12 else NULL)]
    colnames(subset) <- c("Metabolites", paste0(method, " Beta"), paste0(method, " P-Value"), if (method == "Inverse variance weighted") "Odd Ratio" else NULL)
    subset$Metabolites <- paste0(sample_type, " ", subset$Metabolites)
    result_list[[method]] <- subset
  }
  
  return(result_list)
}

# Process both datasets
data1 <- process_data(result_all, "blood")
data2 <- process_data(result_all2, "csf")

# Combine and merge the data
combined_data <- list()
for (method in names(data1)) {
  combined_data[[method]] <- rbind(data1[[method]], data2[[method]])
}

# Merge all dataframes by "Metabolites"
data <- Reduce(function(x, y) merge(x, y, by = "Metabolites"), combined_data)

# Filter the data based on IVW P-Value < 0.05 and consistent sign of Beta values
data <- data[data[,"Inverse variance weighted P-Value"] < 0.05,]

# Ensure all Beta values have the same sign
consistent_sign <- sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"MR Egger Beta"]) &
  sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"Weighted median Beta"]) &
  sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"Simple mode Beta"]) &
  sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"Weighted mode Beta"])

data <- data[consistent_sign,]

# Set row names and remove the Metabolites column
rownames(data) <- data[,1]

# List of metabolites to be removed
metabolites_to_remove <- c("csf Glycerophosphoinositol levels",
                           "blood Adenosine 5'-diphosphate (ADP) to 2'-deoxyuridine ratio",
                           "blood Glycosyl ceramide (d18:1/20:0, d16:1/22:0) levels")

# Remove specified rows from the data
data <- subset(data, !(rownames(data) %in% metabolites_to_remove))

# Remove the Metabolites column
data <- data[,-1]

# Convert to matrix
data <- as.matrix(data)

# Display the filtered data
print(data)

# Prepare parameters
pdf(file = "5.pdf", width = 15, height = 15)
# Set color
pal <- brewer.pal(11,'Spectral')
col_pval = colorRamp2(c(0, 0.05, 1), c("#9E0142","#FFFFBF","#5E4FA2"))

# Group
split_bac <- factor(str_split(rownames(data)," ", simplify = TRUE)[,1], levels = unique(str_split(rownames(data)," ", simplify = TRUE)[,1]))

circos.clear()
circos.par(start.degree = 90, gap.degree = 2, gap.after = c(2,25), track.height = 0.1, clock.wise = FALSE, circle.margin = 0.75)
circos.heatmap.initialize(data, split = split_bac, cluster = FALSE)

# IVW column P-values
circos.heatmap(data[,2], col = col_pval, track.height = 0.05)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(data)[2]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), 
                CELL_META$ycenter, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside", niceFacing = FALSE)
  }
}, bg.border = NA)

# Remove row names
# Assume data and CELL_META are correctly defined
# Use the modified row names in circos.track
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  # Extract row names, filter and sort by CELL_META$subset and CELL_META$row_order
  original_names <- rownames(data)[CELL_META$subset]
  sorted_names <- original_names[CELL_META$row_order]
  
  # Remove the first word from each row name
  # Match everything from the start to the first space (including space) and replace with an empty string
  trimmed_names <- str_replace(sorted_names, "^\\w+\\s+", "")
  
  # Draw text
  circos.text(seq_along(trimmed_names)-0.5,  
              CELL_META$cell.ylim[2]+(convert_y(strwidth(trimmed_names, units = "inches"), "inches"))/4+convert_y(2,"mm"),  
              trimmed_names,  
              cex = 0.5,  
              facing = "reverse.clockwise")  
}, cell.padding = c(0.02, 0, 0.02, 0))

# MR Egger column P-values
circos.heatmap(data[,5], col = col_pval, bg.border = "black", bg.lwd = 1, bg.lty = 1, track.height = 0.05)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(data)[5]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), 
                CELL_META$ycenter, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside", niceFacing = FALSE)
  }
}, bg.border = NA)

# WM column P-values
circos.heatmap(data[,7], col = col_pval, bg.border = "black", bg.lwd = 1, bg.lty = 1, track.height = 0.05)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(data)[7]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), 
                CELL_META$ycenter, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside", niceFacing = FALSE)
  }
}, bg.border = NA)

# IVW OR values column
row_or = data[,3]
circos.track(ylim = range(row_or), panel.fun = function(x, y) {
  y = row_or[CELL_META$subset]
  y = y[CELL_META$row_order]
  circos.lines(CELL_META$cell.xlim, c(1,1), lty = 2, col = "#a5a7ab")
  circos.points(seq_along(y)-0.5, y, col = ifelse(y < 1, "#005f81", "#b03d26"), cex = 0.5, pch = 16)
}, cell.padding = c(0.02, 0, 0.02, 0))
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(data)[3]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), 
                CELL_META$ycenter, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside", niceFacing = FALSE)
  }
}, bg.border = NA)
circos.yaxis(side = "left",
             