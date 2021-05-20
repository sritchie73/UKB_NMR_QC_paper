library(data.table)
library(ukbnmr)
library(MASS)
library(ggplot2)
source("src/utilities/logit.R")

# Create output directory
if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Apply the QC pipeline on all biomarkers, not just the non-derived ones, so we can
# find those most different.
new <- dat[, .(sample_id, visit, sample_degredation, plate_row, plate_column, plate_id, plate_measure_bin, spectrometer, variable, raw)]

# Percentages will be logit transformed instead of log transfomred
new <- new[!is.na(raw)]
new[variable %in% ukbnmr::nmr_info[Units == "%", Biomarker], raw := raw/100]

# Log/logit offset for 0 values
lower_lim <- new[!is.na(raw),.(min=min(raw), non_zero=min(raw[raw != 0])),by=variable]
lower_lim[, log_offset := ifelse(min == 0, non_zero / 2, 0)]
new[lower_lim, on = .(variable), raw := raw + log_offset]

# Log / logit transformation
new[variable %in% nmr_info[Units == "%", Biomarker], log_raw := logit(raw)]
new[variable %in% nmr_info[Units != "%", Biomarker], log_raw := log(raw)]

# Function for converting a set of numbers to a factor, ordering levels by group size
factor_by_size <- function(x) {
  group_sizes <- as.data.table(table(x))
  group_sizes <- group_sizes[order(-N)]
  factor(x, levels=group_sizes$x)
}

# Do the four adjustments
new[, adj := rlm(log_raw ~ sample_degredation)$residuals, by=.(variable)]
new[, adj := rlm(adj ~ factor_by_size(plate_row))$residuals, by=.(variable)]
new[, adj := rlm(adj ~ factor_by_size(plate_column))$residuals, by=.(variable)]
new[, adj := rlm(adj ~ factor_by_size(plate_measure_bin))$residuals, by=.(variable, spectrometer)]

# Rescale the adjusted residuals back to absolute units
rescale <- function(log_adj, log_raw) {
  # The residuals from any regression are defined as the difference between the
  # observed independent variable (e.g. biomarker concentrations) and the parameter
  # estimated by the regression. A resulting key property is that the distribution
  # of the residuals is centred on 0 for this estimated parameter. In the case of
  # robust linear regression, this parameter is an estimate of the mean that is
  # robust to outliers. As a consequence of the way residuals are defined, their
  # distribution can be scaled to match the distribution of the independent variable
  # by giving it the same estimated mean. For robust linear regression, residuals
  # can be returned to the same scale as the observed independent variable by
  # estimating the mean of the observed independent variable using robust linear
  # regression and adding it to the residuals.
  adj <- log_adj + as.vector(coef(rlm(log_raw ~ 1)))[1]
  return(adj)
}
new[, adj := rescale(adj, log_raw), by=.(variable)]
new[variable %in% nmr_info[Units == "%", Biomarker], adj := invlogit(adj)]
new[variable %in% nmr_info[Units != "%", Biomarker], adj := exp(adj)]
new[lower_lim, on = .(variable), adj := adj - log_offset]

# Add offset so that there are no negative values
shift <- new[,.(right_shift=-pmin(0, min(adj))), by=variable]
lower_lim <- lower_lim[shift, on = .(variable)]
new[lower_lim, on = .(variable),  adj := adj + right_shift]

# Percentages back to 0-100 instead of 0 to 1
new[variable %in% nmr_info[Units == "%", Biomarker], adj := adj * 100]

# Add to multi-step qc table
dat[new, on = .(sample_id, visit, variable), adj4_no_rederive := i.adj]

# Write out updated 
fwrite(dat ,sep="\t", quote=FALSE, file="data/tech_qc/multistep_adjusted_values.txt")

# Show difference between methods:
if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")

gg_dt <- dat[variable %in% nmr_info[Type != "Non-derived" & (Nightingale), Biomarker]]
gg_dt <- gg_dt[!is.na(adj4) & !is.na(adj4_no_rederive)]
gg_dt <- gg_dt[order(tolower(variable))]
gg_dt[, variable := factor(variable, levels=unique(variable))]

g <- ggplot(gg_dt) +
  aes(x=adj4, y=adj4_no_rederive) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_x_continuous(name="Derived and composite biomarkers computed from adjusted parts") +
  scale_y_continuous(name="Biomarker adjusted for technical covariates") +
  facet_wrap(~ variable, ncol=12, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=16.8, height=16.8, units="in", file="diagnostic_plots/adjusted_vs_from_parts.png")

