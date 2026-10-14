#load data
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# standardise within each training fold during tuning
tune_res <- tune.rcc(X, Y, validation = "Mfold", scale = TRUE)

# use the same scaling choice for the final model
rcc_res <- rcc(X, Y, lambda1 = tune_res$opt.lambda1[1],
               lambda2 = tune_res$opt.lambda2[1], scale = TRUE)

# plot output
plot(tune_res)
