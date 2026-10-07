# Small shared fixtures. Sizes are kept tiny so that the boosting and stability
# selection tests run quickly.

# list of `n` simulated data sets (columns: y, X1..XP)
sim_data <- function(n = 2, N = 60, P = 12, p_ref = 3, tau = -1, sigma = 2,
                     rho = 0.6, seed = 1) {
    withr::with_seed(seed, lapply(seq_len(n), function(i) {
        get_simulated_data(N = N, P = P, p_ref = p_ref, tau = tau,
                           sigma = sigma, rho = rho)
    }))
}

# stratified train/test split of simulated data sets
cv_fixture <- function(..., seed = 2) {
    withr::with_seed(seed, lapply(sim_data(...), get_stratified_cv_data))
}

# list of `n` data sets with a continuous (age), a binary (sex) and a noise
# predictor, for gamboost
gam_fixture <- function(n = 2, N = 80, seed = 3) {
    set.seed(seed)
    lapply(seq_len(n), function(i) {
        age <- rnorm(N)
        sex <- rbinom(N, 1, 0.5)
        x1 <- rnorm(N)
        y <- rbinom(N, 1, plogis(1.5 * age + 0.8 * sex))
        # mboost expects centred covariates (intercept = FALSE)
        age <- age - mean(age)
        sex <- sex - mean(sex)
        x1 <- x1 - mean(x1)
        data.table(y = factor(y, levels = 0:1), age = age, sex = sex, x1 = x1)
    })
}

# data resembling the FaPaCa data: Status (factor), features and an ID column
fapaca_like_data <- function(n = 60, seed = 4) {
    set.seed(seed)
    data.table(
        ID     = seq_len(n),
        Status = factor(rep(c("C", "L"), times = c(2 * n / 3, n / 3))),
        age    = rnorm(n, 50, 10),
        sex    = rbinom(n, 1, 0.5),
        prot1  = rnorm(n),
        prot2  = rnorm(n)
    )
}
