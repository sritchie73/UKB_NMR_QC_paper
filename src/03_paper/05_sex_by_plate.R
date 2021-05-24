library(data.table)
library(ggplot2)

sinfo <- fread("data/tech_qc/sample_information.txt")

# Compute % female samples on each plate
sinfo <- sinfo[!(sample_removed)]
plate_info <- sinfo[, .(pct_female = sum(Gender == "Female")/.N*100), by=.(plate_id)]

# Into bins of 5%
plate_info[, pct_bin := sprintf("%s - %s%%", floor(pct_female/5) * 5, floor(pct_female/5) * 5 + 5)]
plate_info[, bin_loc := floor(pct_female/5) * 5 + 2.5]
plate_info <- plate_info[, .N, by=.(pct_bin, bin_loc)]

g <- ggplot(plate_info) +
  aes(x = bin_loc, y=N) +
  geom_col(fill="#184068") + 
  scale_x_continuous(name="% samples from female participants on each plate", limits=c(30,75), expand=c(0,0)) +
  scale_y_log10(name="Number of plates") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7))
ggsave(g, width=3.6, height=3.6, file="paper_output/sex_by_plate.pdf")
