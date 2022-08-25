library(ExactMultinom)

num_loci <- 10000
exact_multi_nom_1rep <- function(path) {
    gt_obs_exp <- read.csv(file = path)
    # print(gt_obs_exp)
    # print(gt_obs_exp$obs*10000)
    pvalue <- multinom.test(
        round(gt_obs_exp$prob_obs*num_loci),
        gt_obs_exp$prob_exp,
        stat = "Prob",
        # method = "exact",
        method = "Monte-Carlo",
        theta = 1e-07,
        timelimit = 10,
        N = 10000
        )
    print(pvalue$pvals_mc[1])
    return(pvalue$pvals_mc[1])
}

# res <- exact_multi_nom_1rep(
#     "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/1/heter/gt_prob_true_obs_exp.csv",
#     10000)
# print(res$pvals_mc[1])

root_path <- "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"

# prob_path <- c()
num_replicate<-10
res_df <- data.frame(scale=c(), locus_length=c(), replicate=c(), true_st=c(), true_gt=c(), exact_mn_mc=c())
scale_list<-c("short", "medium", "long")
llen_list <- c(500, 1000)
st_list <- c("True" = "truest", "False" = "mlst")
gt_list <- c("True" = "truegt", "False" = "iqgt")
for (scale in scale_list) {
    for (i in 1:num_replicate) {
        for (llen in llen_list) {
            print(i)
            for (st in c("True", "False")) {
                for (gt in c("True", "False")) {
                    csv_path <- paste("gt_prob_all",st_list[st], gt_list[gt],"0.csv", sep="_")
                    data_path <- paste(root_path, scale, i, "heter", llen, csv_path, sep = "/")
                    print(data_path)
                    pvalue <- exact_multi_nom_1rep(data_path)
                    print(pvalue)
                    # prob_path <- c(prob_path, data_path)
                    res_df <- rbind(res_df, c(scale, llen, i, st, gt, pvalue))
                }
            }
        }
    }
}

write.csv(res_df, paste(root_path, "exact_mn.csv"), sep="/")
# res <- lapply(prob_path, exact_multi_nom_1rep)
# print(res)

# x<-unlist(res)
# plot=qplot(x, geom="histogram", xlim = c(0,1), binwidth=0.05, xlab = "p-value",col=I("blue")) 
# # ggsave(paste(root_path, "/exact_multinom_true_gt_st.pdf", sep = ""), plot=plot, width = 12, height = 10, dpi = 300)
# ggsave(paste(root_path, "/exact_multinom_true_gt_MLst.pdf", sep = ""), plot=plot, width = 12, height = 10, dpi = 300)
