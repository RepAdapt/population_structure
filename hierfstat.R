##########################################################
#
# Run hierfstat on thinned SNP set when there are discrete populations.
#
# Usage
# -----
# cd <output_directory>
# Rscript hierfstat.R thinned_snp_file samp_to_pop_file
#
# Parameters
# ----------
# thinned_snp_file : Path
#     tab-separated datatable with rows for individuals and columns for SNPs (named SNPs and inds, not indices)
# samp_to_pop_file : Path
#     tab-separated datatable with two columns (arbitrarily named) 
#         - first column is sample ID - this must be the same order as rownames of thinned_snp_file
#         - second column is population ID (must be numeric ≥ 1)
#
# Notes
# -----
# - if run manually, direct output to a file : Rscript hierfstat.R thinned_snp_file samp_to_pop_file > output
#     - the output text should be pretty mimal and most stats are printed
#
# Outline
# -------
# 1. GET INPUT ARGUMENTS
# 2. ADD IN SNP AND SAMPLE NAMES
# 3. CALCULATE BASIC STATS
# 4. CALCULATE PAIRWISE FST
# 5. CALCULATE AVERAGED STATS
#
# TODO
# ----
# make sure samps from samptopop are coerced to be same order as how I'm pulling them - allow for subset of samps
##########################################################

library(hierfstat)

# 1. GET INPUT ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: hierfstat.R <thinned_snp_file> <samp_to_pop_file>")
}
thinned_snp_file <- args[1]
samp_to_pop_file <- args[2]

## read in data
thinned_snps <- read.table(thinned_snp_file)

samp_to_pop <- read.table(samp_to_pop_file, header=T)
stopifnot(all(samp_to_pop[,1] == rownames(thinned_snps)))  # assert sample order

## metadata
name <- sub("\\.(txt)$", "", basename(thinned_snp_file), perl = TRUE)


# 2. ADD IN SNP AND SAMPLE NAMES
thinned_snps[, 'pop'] <- samp_to_pop[,2]
snp_data <- data.frame(pop=samp_to_pop[,2], thinned_snps, check.names=FALSE)


# 3. CALCULATE BASIC STATS
cat(sprintf('\nINFO [%s] calculating hierfstat basic stats', Sys.time()))
stats <- basic.stats(snp_data)

cat(sprintf('\nINFO [%s] overall stats \n', Sys.time()))
print( data.frame(stats$overall))

stats_file <- sprintf('%s_hierfstat_overall_stats.txt', name)
write.table(data.frame(stats$overall), stats_file)
cat(sprintf('\nINFO [%s] basic hierfstat stats saved to: %s\n', Sys.time(), stats_file))


# 4. CALCULATE PAIRWISE FST
cat(sprintf('\nINFO [%s] running hierfstat pairwise FST\n', Sys.time()))
pw_fst <- pairwise.WCfst(snp_data)
cat(sprintf('\nINFO [%s] pairwise FST:\n', Sys.time()))
print(pw_fst)

pw_file <- sprintf('%s_hierfstat_pw_fst.txt', name)
write.table(pw_fst, pw_file, sep='\t', row.names=TRUE, col.names=TRUE)
cat(sprintf('\nINFO [%s] saved pairwise FST to : %s\n', Sys.time(), pw_file))


# 5. CALCULATE AVERAGED STATS
cat(sprintf('\nINFO [%s] calculating hierfstat averaged stats\n', Sys.time()))
mean_Ho  <- colMeans(stats$Ho,  na.rm = TRUE)
mean_Hs  <- colMeans(stats$Hs,  na.rm = TRUE)
mean_Fis <- colMeans(stats$Fis, na.rm = TRUE)

cat(sprintf('\nINFO [%s] mean_Ho :\n', Sys.time()))
print(mean_Ho)

cat(sprintf('\nINFO [%s] mean_Hs\n', Sys.time()))
print(mean_Hs)

cat(sprintf('\nINFO [%s] mean_Fis\n', Sys.time()))
print(mean_Fis)

mean_data <- cbind(mean_Ho, mean_Hs, mean_Fis)
mean_file <- sprintf('%s_hierfstat_mean_Ho_Hs_Fis.txt', name)
write.table(mean_data, mean_file, sep='\t', row.names=TRUE, col.names=TRUE)
cat(sprintf('\nINFO [%s] wrote averaged stats to: %s\n', Sys.time(), mean_file))

print('DONE!')
