rm(list = ls())  # Clear workspace

# Install required packages if not already installed
if(!require("gwasrapidd")) install.packages("gwasrapidd", update = FALSE, ask = FALSE)
if(!require("stringr")) install.packages("stringr", update = FALSE, ask = FALSE)
if(!require("data.table")) install.packages("data.table", update = FALSE, ask = FALSE)
if(!require("devtools")) install.packages("devtools")

# Load libraries
library(TwoSampleMR)
library(MRPRESSO)
library(circlize)
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)

# Load result datasets
result_all <- read.csv("C:/Users/王嘉立/Desktop/2.1400 to PD/1400_PD_res.csv")
result_all2 <- read.csv("C:/Users/王嘉立/Desktop/1.338 to PD/338_PD_res.csv")

# Function to process MR result data
process_data <- function(data, sample_type) {
  methods <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Simple mode", "Weighted mode")
  result_list <- list()
  
  for (method in methods) {
    subset <- data[data$method == method,][, c(4, 7, 9, if (method == "Inverse variance weighted") 12 else NULL)]
    colnames(subset) <- c("Metabolites", paste0(method, " Beta"), paste0(method, " P-Value"), if (method == "Inverse variance weighted") "Odd Ratio" else NULL)
    subset$Metabolites <- paste0(sample_type, " ", subset$Metabolites)
    result_list[[method]] <- subset
  }
  return(result_list)
}

# Function to draw metabolite labels
draw_metabolite_labels <- function(data, cex = 0.75, extra_mm = 2) {
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    original_names <- rownames(data)[CELL_META$subset]
    sorted_names <- original_names[CELL_META$row_order]
    trimmed_names <- stringr::str_replace(sorted_names, "^\\w+\\s+", "")
    
    text_height <- convert_y(strwidth(trimmed_names, units = "inches") * cex, "inches") / 2
    y_positions <- CELL_META$cell.ylim[2] + text_height + convert_y(extra_mm, "mm")
    
    circos.text(seq_along(trimmed_names) - 0.5, y_positions, trimmed_names, cex = cex, facing = "reverse.clockwise", niceFacing = TRUE)
  }, cell.padding = c(0.02, 0, 0.02, 0))
}

# Process both datasets
data1 <- process_data(result_all, "blood")
data2 <- process_data(result_all2, "csf")

# Combine data by method
combined_data <- list()
for (method in names(data1)) {
  combined_data[[method]] <- rbind(data1[[method]], data2[[method]])
}

# Merge dataframes by 'Metabolites'
data <- Reduce(function(x, y) merge(x, y, by = "Metabolites"), combined_data)

# Filter by IVW P-value < 0.05 and consistent Beta sign
consistent_sign <- sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"MR Egger Beta"]) &
  sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"Weighted median Beta"]) &
  sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"Simple mode Beta"]) &
  sign(data[,"Inverse variance weighted Beta"]) == sign(data[,"Weighted mode Beta"])

# Subset filtered data
data <- data[data[,"Inverse variance weighted P-Value"] < 0.05 & consistent_sign,]
rownames(data) <- data[,1]  # Set metabolite names as row names

# Remove specific metabolites from data
metabolites_to_remove <- c("csf Glycerophosphoinositol levels",
                           "blood Adenosine 5'-diphosphate (ADP) to 2'-deoxyuridine ratio",
                           "blood Glycosyl ceramide (d18:1/20:0, d16:1/22:0) levels")
data <- subset(data, !(rownames(data) %in% metabolites_to_remove))

# Remove 'Metabolites' column
data <- data[,-1]

# Convert to matrix
data <- as.matrix(data)

# Define heatmap colors
pal <- brewer.pal(11,'Spectral')
col_pval <- colorRamp2(c(0, 0.05, 1), c("#9E0142","#FFFFBF","#5E4FA2"))

# Define groups for sectors
split_bac <- factor(str_split(rownames(data)," ",simplify = T)[,1], levels = unique(str_split(rownames(data)," ",simplify = T)[,1]))

# Start PDF output
pdf(file = "8.pdf", width = 20, height = 20)

# Initialize Circos heatmap
circos.clear()
circos.par(start.degree = 90, gap.degree = 2, gap.after = c(2,25), track.height = 0.1, clock.wise = F, circle.margin = 1.5)
circos.heatmap.initialize(data, split = split_bac, cluster = F)

colnames(data)[2] <- "IVW P-Value"
colnames(data)[5] <- "MR Egger P-Value"
colnames(data)[7] <- "WM P-Value"
colnames(data)[3] <- "OR"

# Add heatmap for IVW P-values
circos.heatmap(data[,2], col = col_pval, track.height = 0.05)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) {
    cn <- colnames(data)[2]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), CELL_META$ycenter, cn, cex = 0.75, adj = c(0,0.5), facing = "inside")
  }
}, bg.border = NA)

# Draw metabolite labels
draw_metabolite_labels(data, cex = 0.8, extra_mm = 4)

# Add MR-Egger P-value heatmap
circos.heatmap(data[,5], col = col_pval, bg.border = "black", bg.lwd = 1, bg.lty = 1, track.height = 0.05)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) {
    cn <- colnames(data)[5]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), CELL_META$ycenter, cn, cex = 0.75, adj = c(0,0.5), facing = "inside")
  }
}, bg.border = NA)

# Add Weighted Median P-value heatmap
circos.heatmap(data[,7], col = col_pval, bg.border = "black", bg.lwd = 1, bg.lty = 1, track.height = 0.05)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) {
    cn <- colnames(data)[7]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), CELL_META$ycenter, cn, cex = 0.75, adj = c(0,0.5), facing = "inside")
  }
}, bg.border = NA)

# Add Odd Ratio scatter track
row_or <- data[,3]
circos.track(ylim = range(row_or), panel.fun = function(x, y) {
  y <- row_or[CELL_META$subset]
  y <- y[CELL_META$row_order]
  circos.lines(CELL_META$cell.xlim, c(1,1), lty = 2, col = "#a5a7ab")
  circos.points(seq_along(y)-0.5, y, col = ifelse(y < 1, "#005f81", "#b03d26"), cex = 0.5, pch = 16)
}, cell.padding = c(0.02, 0, 0.02, 0))
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) {
    cn <- colnames(data)[3]
    circos.text(CELL_META$cell.xlim[2] + convert_x(2, "mm"), CELL_META$ycenter, cn, cex = 0.75, adj = c(0,0.5), facing = "inside")
  }
}, bg.border = NA)
circos.yaxis(side = "left", labels.cex = 0.5, at = c(0.8, 1, 1.3))

# Add sector labels
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - convert_y(3, "mm"), CELL_META$sector.index, facing = "bending.inside", cex = 0.5, adj = c(0.5, 0))
}, bg.border = NA)

# Add Legends
pval <- Legend(title = "P.Value", col_fun = col_pval, at = c(0, 0.05, 0.5, 1))
or <- Legend(labels = c("Risk factor", "Protective factor"), title = "Odd Ratio", legend_gp = gpar(fill = c("#b03d26", "#005f81")))
all <- packLegend(pval, or)
grid.draw(all)

# Close PDF
dev.off()
