library(ExactMultinom)

exact_multi_nom_1rep <- function(gt_obs_exp, num_loci, method) {
    # gt_obs_exp <- read.csv(file = path)
    pvalue <- multinom.test(
        round(gt_obs_exp$prob_obs*num_loci),
        gt_obs_exp$prob_exp,
        stat = "Prob",
        # method = "exact",
        method = method,
        theta = 1e-07,
        timelimit = 10,
        N = 10000
        )
    # print(pvalue$pvals_mc[1])
    if (method == "Monte-Carlo"){
        return(pvalue$pvals_mc[1])
    }
    else if(method == "exact"){
        return(pvalue$pvals_ex[1])
    }
}