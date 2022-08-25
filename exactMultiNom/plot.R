library(ggplot2)
path<-"/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/statistics.csv"
df <- read.csv(file = path)
print(df)
# ggplot(df, aes(leafNum, Count, fill=dlRate, colour=dlRate)) +
#   labs(x="Size of gene family",colour="Duplication and loss level",fill="Duplication and loss level") +
#   geom_histogram(stat = "identity",alpha=0.6) +
#   scale_fill_grey(start = 0.8, end = 0.2) +
#   scale_colour_grey(start = 1, end = 0.5) +
#   geom_vline(aes(xintercept=16), linetype="dashed", size = 0.3) +
#   facet_grid(popSize~type,scales='free') +
#   theme(legend.position = "top") + coord_flip()

p<-ggplot(df, aes(x=chisq)) + geom_histogram(binwidth=1)
# plot(p)
# + geom_histogram(stat = "identity",alpha=0.6) +
#   scale_fill_grey(start = 0.8, end = 0.2) +
#   scale_colour_grey(start = 1, end = 0.5) +
#   geom_vline(aes(xintercept=16), linetype="dashed", size = 0.3) +
#   facet_grid(popSize~type,scales='free') +
#   theme(legend.position = "top") + coord_flip()

# df <- data.frame(
#   sex=factor(rep(c("F", "M"), each=200)),
#   weight=round(c(rnorm(200, mean=55, sd=5), rnorm(200, mean=65, sd=5)))
#   )
# head(df)
# p<-ggplot(df, aes(x=weight)) + geom_histogram()
