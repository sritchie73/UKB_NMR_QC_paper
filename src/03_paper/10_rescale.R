library(data.table)
library(MASS)
library(ggplot2)
library(ggrastr)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step and extract Clinical_LDL_C
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")
dat <- dat[variable == "Clinical_LDL_C" & !is.na(raw)]

# Show how distributions change at each step.

g1 <- ggplot(dat) +
  aes(x = raw) + 
  geom_density(trim=TRUE, color="#e08214", size=0.4) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Clinical LDL cholesterol (mmol/L)") +
  theme_bw() + 
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
			axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
      panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  ) 

g2 <- ggplot(dat) +
  aes(x = log_raw) + 
  geom_density(trim=TRUE, color="#e08214", size=0.4) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="After log transformation") +
  theme_bw() + 
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
			axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
      panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  ) 

g3 <- ggplot(dat) +
  aes(x = log_adj4) + 
  geom_density(trim=TRUE, color="#542788", size=0.4) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Residuals") +
  theme_bw() + 
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
			axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
      panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  ) 

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
dat[, log_adj4_rescale := rescale(log_adj4, log_raw)]

g4 <- ggplot(dat) +
  geom_density(aes(x=log_raw), trim=TRUE, color="#e08214", size=0.4) +
  geom_density(aes(x=log_adj4_rescale), trim=TRUE, color="#542788", size=0.4) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Rescaled") +
  theme_bw() + 
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
			axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
      panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  ) 

g5 <- ggplot(dat) +
  geom_density(aes(x=raw), trim=TRUE, color="#e08214", size=0.4) +
  geom_density(aes(x=postqc), trim=TRUE, color="#542788", size=0.4) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="postqc") +
  theme_bw() + 
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
			axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
      panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  ) 

g6 <- ggplot(dat) +
  aes(x = raw, y = postqc) +
  rasterize(geom_point(shape=21, size=0.45, stroke=0.2, color="black", fill="white")) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_x_continuous(name="Raw (mmol/L)", limits=c(0, 8.5)) +
  scale_y_continuous(name="Post-QC (mmol/L)", limits=c(0, 8.5)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  )

g <- plot_grid(g1, g2, g3, g4, g5, g6, nrow=2, align="hv") 
ggsave(g, width=7.2, height=4.8, file="paper_output/rescale.pdf")

