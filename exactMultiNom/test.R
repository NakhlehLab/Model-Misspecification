# install.packages("EMT")
library(EMT)
observed =c(A = 5, B = 2, C = 2)
# prob =rep(1/length(observed),3)
prob=c(0.7,0.15,0.15)
out <- multinomial.test(observed, prob, useChisq = FALSE, MonteCarlo = FALSE)


library(ExactMultinom)
out2 <- multinom.test(
        observed,
        prob,
        stat = "Prob",
        method = "exact",
        # theta = 1e-05,
        # timelimit = 10,
        )
print(out)
print(out2)

# Test fairness of a die (that is, whether each side has the same probability)
# p_fair = rep(1/6,6) # Hypothesized probabilities for each side
# x = c(16,17,12,15,15,25) # Observed number of times each side appeared on 100 throws
# # Exact multinomial test (using probability ordering by default):
# multinom.test(x,p_fair)
# # Exact multinomial test using log-likelihood ratio:
# multinom.test(x,p_fair,stat = "LLR")
# # Classical chi-square test (using asymptotics to estimate p-value and Pearson's chi-square):
# multinom.test(x,p_fair,stat = "Chisq",method = "asymptotic")
# # Using Monte-Carlo approach and probability ordering
# multinom.test(x,p_fair,stat = "Prob",method = "Monte-Carlo")