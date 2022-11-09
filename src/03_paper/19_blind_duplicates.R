library(data.table)
library(ggplot2)

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Filter to blind duplicate samples
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (blind_duplicate)]
bd <- dat[unique(sinfo[,.(eid_7439, visit)]), on = .(eid_7439, visit), nomatch=0]

# Use the sample kept for UKB release as the reference sample
bd[, blind_duplicate := "duplicate"]
bd[sinfo[(in_ukb_raw)], on = .(eid_7439, visit, sample_id), blind_duplicate := "reference"]
bd <- bd[order(sample_id)]
bd[blind_duplicate == "duplicate", blind_duplicate := paste0("duplicate.", 1:.N), by=.(eid_7439, visit, variable)]

# extract original data in format for computing CV% and R2
raw_bd <- dcast(bd, eid_7439 + visit + variable ~ blind_duplicate, value.var="raw")
raw_bd <- melt(raw_bd, measure.vars=c("duplicate.1", "duplicate.2"), na.rm=TRUE, variable.name="dup.number", value.name="duplicate")
raw_bd <- raw_bd[!is.na(reference) & !is.na(duplicate)]

# extract post-qc data in format for computing CV% and R2
adj4_bd <- dcast(bd, eid_7439 + visit + variable ~ blind_duplicate, value.var="adj4")
adj4_bd <- melt(adj4_bd, measure.vars=c("duplicate.1", "duplicate.2"), na.rm=TRUE, variable.name="dup.number", value.name="duplicate")
adj4_bd <- adj4_bd[!is.na(reference) & !is.na(duplicate)]

# extract post-qc data with outlier plates excluded in format for computing CV% and R2
adj5_bd <- dcast(bd, eid_7439 + visit + variable ~ blind_duplicate, value.var="adj5")
adj5_bd <- melt(adj5_bd, measure.vars=c("duplicate.1", "duplicate.2"), na.rm=TRUE, variable.name="dup.number", value.name="duplicate")
adj5_bd <- adj5_bd[!is.na(reference) & !is.na(duplicate)]

# Function to compute R2 and CV% when comparing repeated measurements, following
#
# Bland, J.M., and Altman, D.G. (1996). Measurement error proportional to the mean. BMJ 313, 106..
#
# https://www-users.york.ac.uk/~mb55/meas/cv.htm#:~:text=In%20the%20study%20of%20measurement,within%2Dsubject%20coefficient%20of%20variation
replicability <- function(x1, x2) {
  # keep only pairwise complete samples
  non_missing <- !is.na(x1) & !is.na(x2)
  x1 <- x1[non_missing]
  x2 <- x2[non_missing]

  # Add small offset for 0 values 
  if (any(x1 == 0) || any(x2 == 0)) {
    offset <- pmin(min(x1[x1 != 0]), min(x2[x2 != 0]))/2
    x1 <- x1 + offset
    x2 <- x2 + offset
  } 

  # Remove outliers more than 4 x IQR from median
  xmed <- median(c(x1, x2))
  lqr <- quantile(c(x1, x2), 0.25)
  uqr <- quantile(c(x1, x2), 0.75)
  xlim <- c(xmed - (xmed - lqr)*4, xmed + (uqr - xmed)*4)
  outliers <- x1 < xlim[1] | x1 > xlim[2] | x2 < xlim[1] | x2 > xlim[2]
  x1 <- x1[!(outliers)]
  x2 <- x2[!(outliers)]

  # Compute different versions of CV%
   
  # Root square mean approach
  s2 <- ((x1 - x2)**2)/2
  m <- (x1 + x2)/2
  s2m2 <- s2/(m**2)
  cv_sqrt <- sqrt(mean(s2m2))

  # Log method
  lx1 <- log(x1)
  lx2 <- log(x2)
  s2l <- ((lx1 - lx2)**2)/2
  cv_log <- exp(sqrt(mean(s2l)))-1

  # Mean CV (not recommended, but I suspect this is what Nightingale uses as the numbers lien up best)
  cv_mean <- mean(sqrt(s2)/m)

  # Whole dataset comparison rather than per subject
  cv_dataset <- mean(s2) / ( (mean(x1) + mean(x2))/2 )

  # compute different versions of R2
  r2 <- summary(lm(x2 ~ x1))$r.squared
  r2_log <- summary(lm(lx2 ~ lx1))$r.squared
  
  list(cv_sqrt=cv_sqrt, cv_log=cv_log, cv_mean=cv_mean, cv_dataset=cv_dataset, r2=r2, r2_log=r2_log)
}

# Compute repeatability statistics
bd_stats <- rbind(idcol="dataset",
  raw=raw_bd[,replicability(reference, duplicate),by=variable],
  adj4=adj4_bd[,replicability(reference, duplicate),by=variable],
  adj5=adj5_bd[,replicability(reference, duplicate),by=variable]
)
bd_stats <- melt(bd_stats, id.vars=c("dataset", "variable"), variable.name="metric")
bd_stats <- dcast(bd_stats, variable + metric ~ dataset)

# Write out
fwrite(bd_stats, sep="\t", quote=FALSE, file="paper_output/blind_duplicate_cv_r2.txt")

# Plot comparison for non-derived biomarkers
g <- ggplot(bd_stats[metric %in% c("cv_sqrt", "r2") & variable %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]) +
  aes(x=raw, y=adj5) +
  geom_abline(intercept=0, slope=1, linetype=2) + 
  geom_point(shape=19, alpha=0.6, size=0.6) +
  facet_wrap(~ metric, scales="free") +
  scale_x_continuous("Original dataset", oob=scales::squish) +
  scale_y_continuous("Post-QC", oob=scales::squish) +
  theme_bw() +
  theme(axis.title=element_text(size=7), axis.text=element_text(size=6),
        strip.text=element_text(size=7))

ggsave(g, height=1.5, width=2.7, file="paper_output/blind_duplicate_cv_r2.pdf")
















