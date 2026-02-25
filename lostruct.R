#############################################################################################
#
# Run lostruct on full set of imputed + maf-filtered snps.
#
# Usage
# -----
# cd <output_directory>
# Rscript lostruct.R imputed_maf_filt_snp_file
# 
# Parameters
# ----------
# imputed_snp_file : Path
#    path to full snp data set (imputed, maf-filtered, but not thinned)
#    tab-separated datatable (.txt or .lfmm) with rows for individuals and columns for SNPs
#
#############################################################################################
library(lostruct)

# GET INPUT ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Usage: lostruct.R <imputed_maf_filt_snp_file>")
}
snp_file <- args[1]
name <- sub("\\.(lfmm|txt)$", "", basename(snp_file), perl = TRUE)

## read in data
cat(sprintf('\nINFO [%s] reading in SNPs\n', Sys.time()))
snps <- t(as.matrix(read.table(snp_file, header=T)))  # transpose so inds are columns and snps are rows

# RUN LOSTRUCT
cat(sprintf('\nINFO [%s] calculating eigen windows\n', Sys.time()))
ew <- eigen_windows(snps, win = 1000, k = 2)

cat(sprintf('\nINFO [%s] calculating distance matrix\n', Sys.time()))
dist_mat <- pc_dist(ew)

cat(sprintf('\nINFO [%s] running multidimensional scaling\n', Sys.time()))
fit2d <- cmdscale(dist_mat, eig = TRUE, k = 2)

mds_fig_file <- sprintf('%s_mds_plot.pdf', name)
pdf(mds_fig_file, height=4, width=8)
plot(fit2d$points, xlab = "Coordinate 1", ylab = "Coordinate 2", col=rainbow(1.2 * nrow(dist_mat)) )
dev.off()
cat(sprintf('\nINFO [%s] saved MDS plot to: %s\n', Sys.time(), mds_fig_file))

mds_file <- sprintf('%s_mds_points.txt', name)
write.table(fit2d$points, mds_file, row.names=FALSE, col.names=FALSE)
cat(sprintf('\nINFO [%s] saved MDS data to: %s\n', Sys.time(), mds_file))

print('DONE!')
