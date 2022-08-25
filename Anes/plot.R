library(ggplot2)
library(scales)
library(cowplot)
library(ggbeeswarm)

remove(list = ls())

# locus_length=2000
quarnet_results <- read.csv(file = '/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/quarnetGoF_iq_outlierp.csv', header = TRUE, check.names = TRUE)
quarnet_results_true <- read.csv(file = '/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/quarnetGoF_true_outlierp.csv', header = TRUE, check.names = TRUE)

quarnet_results_500<-quarnet_results[quarnet_results$locus_length == 500 & quarnet_results$t1 != "Z" & quarnet_results$t2 != "Z" & quarnet_results$t3 != "Z" & quarnet_results$t4 != "Z", ]
quarnet_results_2000<-quarnet_results[quarnet_results$locus_length == 2000 & quarnet_results$t1 != "Z" & quarnet_results$t2 != "Z" & quarnet_results$t3 != "Z" & quarnet_results$t4 != "Z", ]


hist_short_true<-ggplot(quarnet_results_true[quarnet_results_true$scale == "short", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black")) +
 labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.1, 1.1))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_short_iq500<-ggplot(quarnet_results_500[quarnet_results_500$scale == "short", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_short_iq2000<-ggplot(quarnet_results_2000[quarnet_results_2000$scale == "short", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_medium_true<-ggplot(quarnet_results_true[quarnet_results_true$scale == "medium", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_medium_iq500<-ggplot(quarnet_results_500[quarnet_results_500$scale == "medium", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_medium_iq2000<-ggplot(quarnet_results_2000[quarnet_results_2000$scale == "medium", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_long_true<-ggplot(quarnet_results_true[quarnet_results_true$scale == "long", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_long_iq2000<-ggplot(quarnet_results_2000[quarnet_results_2000$scale == "long", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  # scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

hist_long_iq500<-ggplot(quarnet_results_500[quarnet_results_500$scale == "long", ], aes(x=p_value)) +
  geom_histogram(binwidth=0.05, color="#e9ecef", alpha=0.9, aes(y=(..count..)/sum(..count..))) +
  geom_vline(aes(xintercept=0.05), linetype="dashed", size = 0.3) +
  theme_bw(base_size = 19) + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black")) +
  labs(x = expression('p-values'), y = expression(paste('Percentage '))) +
  scale_x_continuous(limits = c(-0.05, 1.05))+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "y", minor = "none")

sig_pro = c()
sig_pro <- append(sig_pro, sum(subset(quarnet_results_true, scale == "short")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_2000, scale == "short")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_500, scale == "short")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_true, scale == "medium")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_2000, scale == "medium")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_500, scale == "medium")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_true, scale == "long")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_2000, scale == "long")$p_value < 0.05)/100*100)
sig_pro <- append(sig_pro, sum(subset(quarnet_results_500, scale == "long")$p_value < 0.05)/100*100)


prow = plot_grid(hist_short_true + theme(legend.position = "none"),
                 hist_short_iq2000 + theme(legend.position = "none"),
                 hist_short_iq500 + theme(legend.position = "none"),
                 hist_medium_true + theme(legend.position = "none"),
                 hist_medium_iq2000 + theme(legend.position = "none"),
                 hist_medium_iq500 + theme(legend.position = "none"),
                 hist_long_true + theme(legend.position = "none"),
                 hist_long_iq2000 + theme(legend.position = "none"),
                 hist_long_iq500 + theme(legend.position = "none"),
                 labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), nrow = 3, align = "hv", rel_heights = c(0.8), vjust = c(1.0, 1.0, 1.0))
legend = get_legend(hist_long_true)

p<-plot_grid(legend, prow, ncol = 1, rel_heights = c(0.07, 1))
ggsave(paste0("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/quarnet_outlier_pval",".pdf"), units = "mm", height = 200, width = 300)
