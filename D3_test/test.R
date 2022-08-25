library(ggplot2)
library(scales)
library(cowplot)
library(ggbeeswarm)

#remove(list = ls())

locus_length=500
#D3_gamma_bootstrap_results <- read.csv(file = '/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/D3_res.csv', header = TRUE, check.names = TRUE)
D3_gamma_bootstrap_results <- read.csv(file = '/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/tools/D3_introgression/D3_res_gamma_0.csv', header = TRUE, check.names = TRUE)

### P value calculcation for individual sims ### 

D3_gamma_bootstrap_results<-D3_gamma_bootstrap_results[D3_gamma_bootstrap_results$locus_length == locus_length & D3_gamma_bootstrap_results$A != "Z_0" & D3_gamma_bootstrap_results$B != "Z_0" & D3_gamma_bootstrap_results$C != "Z_0", ]
get_pval <- function(estimate, stdev){ #gets pval from statistic value and stdev
  
  z <- estimate / stdev #z score
  p <- pnorm(estimate, mean = 0, sd = stdev, lower.tail = FALSE)
  #p <- pnorm(estimate, mean = 0, sd = stdev, lower.tail = TRUE) #p value
  #p <- 2 *pnorm(abs(estimate), mean = 0, sd = stdev, lower.tail = FALSE)
  #p<-pnorm(z)
  return(p)
}


D3_gamma_bootstrap_results$D3_pvals <- mapply(get_pval, #calculate D3 p-values
                                              estimate = D3_gamma_bootstrap_results$D3,
                                              stdev = D3_gamma_bootstrap_results$D3_stdev)

### Plot for D3 ###

D3_plot <- ggplot(D3_gamma_bootstrap_results, #setup plot object
                  aes(x = as.factor(scale), y = D3))

#Add violin
D3_plot <- D3_plot + geom_violin(trim=FALSE) 
#+ geom_boxplot(width=0.05)+stat_summary(fun.data=mean_sdl, mult=0.1, geom="pointrange", color="red")

#Add points
#D3_plot <- D3_plot + geom_jitter(position = position_jitter(0.01))
#Remove grid background and add borders

D3_plot <- D3_plot + theme_bw(base_size = 9) + theme(panel.border = element_blank(), 
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(), 
                                                      axis.line = element_line(colour = "black"))

#Add dashed line at 0

D3_plot <- D3_plot + geom_hline(aes(yintercept = 0), 
                                linetype = "dashed")

#Adjust axis labels

D3_plot <- D3_plot + labs(x = expression(paste(
  'Scale')),
  y = expression(italic('D')[3]))

ggsave(paste0("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/D3_",locus_length,".pdf"), plot=D3_plot, units = "mm", height = 40, width = 50)

### Calculate rejection rates ### 
#
 D3_short_power <- sum(subset(D3_gamma_bootstrap_results, scale == "short")$D3_pvals < 0.05)/100*100
 D3_medium_power <- sum(subset(D3_gamma_bootstrap_results,scale == "medium")$D3_pvals < 0.05)/100*100
 D3_long_power <- sum(subset(D3_gamma_bootstrap_results, scale == "long")$D3_pvals < 0.05)/100*100


D3_hist_short<-ggplot(D3_gamma_bootstrap_results[D3_gamma_bootstrap_results$scale == "short", ], aes(x=D3_pvals)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black")) +
 labs(x = expression('D3 p-values'), y = expression(paste('Percentage '))) +
  scale_x_continuous(limits = c(-0.05, 1.05))+ scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

D3_hist_medium<-ggplot(D3_gamma_bootstrap_results[D3_gamma_bootstrap_results$scale == "medium", ], aes(x=D3_pvals)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('D3 p-values'), y = expression(paste('Percentage '))) +
  scale_x_continuous(limits = c(-0.05, 1.05))+ scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

D3_hist_long<-ggplot(D3_gamma_bootstrap_results[D3_gamma_bootstrap_results$scale == "long", ], aes(x=D3_pvals)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('D3 p-values'), y = expression(paste('Percentage '))) +
  scale_x_continuous(limits = c(-0.05, 1.05))+ scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")


prow = plot_grid(D3_hist_short + theme(legend.position = "none"),
                 D3_hist_medium + theme(legend.position = "none"),
                 D3_hist_long + theme(legend.position = "none"),
                 labels = c("A", "B", "C"), nrow = 1, align = "hv", rel_heights = c(0.8), vjust = c(1.0, 1.0, 1.0))
legend = get_legend(D3_hist_long)

plot_grid(legend, prow, ncol = 1, rel_heights = c(0.07, 1))
ggsave(paste0("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/tools/D3_introgression/D3_res_gamma_0_pval.pdf"), units = "mm", height = 50, width = 200)
