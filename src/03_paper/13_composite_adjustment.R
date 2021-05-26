library(data.table)
library(ukbnmr)
library(MASS)
library(foreach)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(rcartocolor)
library(ggnewscale)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in
# UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Show example where biomarker concentrations differ after adjustment
gg_dt <- dat[variable == "XXL_VLDL_PL_pct" & !is.na(adj4) & !is.na(adj4_no_rederive) & !is.na(adj5)]

g1 <- ggplot(gg_dt) +
  aes(x=adj4, y=adj4_no_rederive) +
  rasterize(geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2)) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.2) +
  scale_x_continuous(name="Post-QC value") +
  scale_y_continuous(name="Adjusted directly") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# Show effects of sample degradation time on XL_HDL_L
effects <- foreach(var = c("XL_HDL_FC", "XL_HDL_CE", "XL_HDL_PL", "XL_HDL_TG"), .combine=rbind) %do% {
  l1 <- rlm(scale(log(raw)) ~ sample_degredation, data=dat[variable == var])
  cf <- coef(summary(l1))
  ci <- confint.default(l1)
  data.table(variable = var, beta = cf[2,1], l95 = ci[2,1], u95 = ci[2,2])
}

g2 <- ggplot(effects) +
  aes(y=variable, x=beta, xmin=l95, xmax=u95) +
  geom_errorbarh(height=0, size=0.2) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_vline(xintercept=0, linetype=2, color="red", size=0.2) +
  scale_x_continuous("SD change in log-biomarker", limits=c(-0.001, 0.0035), breaks=c(-0.001, 0, 0.003)) + 
  ylab("") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# Show differences as a function of sample degredation time
gg_dt <- dat[variable == "XL_HDL_L" & !is.na(adj4) & !is.na(adj4_no_rederive) & !is.na(adj5)]
gg_dt <- gg_dt[order(sample_degredation)]

pal <- colorRampPalette(carto_pal(name="ag_Sunset"))(255)
g3 <- ggplot(gg_dt) +
  aes(x=adj4, y=adj4 - adj4_no_rederive, fill=sample_degredation) +
  rasterize(geom_point(shape=21, color="black", size=0.45, stroke=0.2)) +
  geom_hline(yintercept=0, linetype=2, color="red", size=0.2) +
  scale_x_continuous(name="Post-QC value", limits=c(0,1), breaks=c(0, 0.5, 1)) +
  scale_y_continuous(name="Adj - postQC", limits=c(-0.005, 0.015), breaks=c(-0.005, 0, 0.005, 0.01, 0.015)) +
  scale_fill_gradientn(name="Sample degradation", colors=pal, limits=c(0, 151)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=7), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none"
  ) 


# Show more dramatic example using age, sex, and bmi.
asb <- fread("data/age_sex_bmi_adj/adjusted_vs_recomputed.txt")
asb[, adj_no_rederive := as.numeric(adj_no_rederive)]

gg_dt <- asb[variable == "HDL_PL" & !is.na(adj) & !is.na(adj_no_rederive)]
 
g4 <- ggplot(gg_dt) + 
  aes(x = adj, y = adj_no_rederive) +
  rasterize(geom_point(color="black", fill="white", size=0.45, stroke=0.2, shape=21)) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.2) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
   )

# Show effects on component parts
dat[unique(asb, by=c("eid_7439", "sample_id", "visit")), on = .(eid_7439, sample_id, visit), c("age", "sex", "bmi") := .(age, sex, bmi)]

effects <- foreach(var = c("S_HDL_PL", "M_HDL_PL", "L_HDL_PL", "XL_HDL_PL"), .combine=rbind) %do% {
  l1 <- rlm(scale(adj5) ~ scale(age) + factor(sex, levels=c("Female", "Male")) + scale(log(bmi)), data=dat[variable == var])
  cf <- coef(summary(l1))
  ci <- confint.default(l1)
  data.table(variable = var, covariate = c("age", "sex (male)", "bmi"), beta = cf[-1,1], l95 = ci[-1,1], u95 = ci[-1,2])
}

g5 <- ggplot(effects) +
  aes(y=variable, x=beta, xmin=l95, xmax=u95, color=covariate, fill=covariate) +
  geom_errorbarh(height=0, size=0.4, position=position_dodge(width=0.6)) +
  geom_point(shape=21, color="black", size=0.45, stroke=0.2, position=position_dodge(width=0.6)) +
  scale_fill_manual(name="", values=c("sex (male)"="#D86020", "age"="#B1A866", "bmi"="#572530")) +
  scale_color_manual(name="", values=c("sex (male)"="#D86020", "age"="#B1A866", "bmi"="#572530")) +
  geom_vline(xintercept=0, linetype=2, color="red", size=0.2) +
  scale_x_continuous("SD change in log-biomarker", limits=c(-0.81, 0.2), breaks=c(-0.5, -0.2, 0, 0.2)) + 
  ylab("") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none")

g <- plot_grid(g1, g2, g3, g4, g5, nrow=1)

ggsave(g, width=7.2, height=1.5, file="paper_output/composite_adjustment.pdf")

# Also plot version of g3 with legend so we can extract for plot
gg_dt <- dat[variable == "XL_HDL_L" & !is.na(adj4) & !is.na(adj4_no_rederive) & !is.na(adj5)]
gg_dt <- gg_dt[order(sample_degredation)]

pal <- colorRampPalette(carto_pal(name="ag_Sunset"))(255)
g3 <- ggplot(gg_dt) +
  aes(x=adj4, y=adj4 - adj4_no_rederive, fill=sample_degredation) +
  rasterize(geom_point(shape=21, color="black", size=0.45, stroke=0.2)) +
  geom_hline(yintercept=0, linetype=2, color="red", size=0.2) +
  scale_x_continuous(name="Post-QC value", limits=c(0,1), breaks=c(0, 0.5, 1)) +
  scale_y_continuous(name="Adj - postQC", limits=c(-0.0051, 0.015), breaks=c(-0.005, 0, 0.005, 0.010, 0.015)) +
  scale_fill_gradientn(name="Sample degradation", colors=pal, limits=c(0, 151)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=7), 
        legend.text=element_text(size=6), legend.title=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom"
  ) + guides(fill = guide_colorbar(title.position="top"))

ggsave(g3, width=1.8, height=1.5, file="paper_output/composite_adjustment_additional_legend.pdf")

