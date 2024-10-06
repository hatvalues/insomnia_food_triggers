library(lattice)
library(ggplot2)
library(dplyr)
library(survival)

dframe <- read.csv("data/insomnia-diary.csv")
month_name <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
dframe$month <- month_name[dframe$month_number]
dframe$quarter <- (dframe$month_number - 1) %/% 3 + 1
dframe[["combo_foods"]] <- factor(with(dframe, cheese * 1 + brassica * 2 + meat * 4 + spice * 8))

dframe[["trend"]] = 1:length(dframe$time_sleeping)
dframe[["ar_1"]] = c(dframe$time_sleeping[1], head(dframe$time_sleeping, -1))
dframe[["ar_2"]] = c(rep(dframe$time_sleeping[1], 2), head(dframe$time_sleeping, -2))
dframe[["ar_3"]] = c(rep(dframe$time_sleeping[1], 3), head(dframe$time_sleeping, -3))
dframe[["ar_4"]] = c(rep(dframe$time_sleeping[1], 4), head(dframe$time_sleeping, -4))


# EXPLORATORY
densityplot(~time_sleeping, data = dframe, groups = dframe$quarter, auto.key = TRUE)
densityplot(~time_sleeping, data = dframe, groups = dframe$combo_foods, auto.key = TRUE)
ggplot(data = dframe, mapping = aes(x = time_sleeping)) +
    geom_density(mapping = aes(colour = cheese)) +
    geom_density(mapping = aes(colour = brassica)) +
    geom_density(mapping = aes(colour = meat)) +
    geom_density(mapping = aes(colour = spice))
#    geom_density(mapping = aes(colour = month))

lm_model <- lm(time_sleeping~.+ I(quarter^2) - trend - month_number - month - combo_foods, data=dframe)

glm_model <- glm(time_sleeping~.+ I(quarter^2) - trend - month_number - month - combo_foods, data= dframe, family = Gamma)
summary(glm_model)
plot(glm_model)


weibull_model <- survreg(Surv(time_sleeping) ~ .+ I(quarter^2) - trend - month_number - month - combo_foods, data = dframe, dist = "weibull")

plot(1/predict(glm_model), predict(weibull_model))
ggplot(dframe, mapping = aes(x = trend)) +
  geom_point(mapping = aes(y = time_sleeping), colour = 1) +
  geom_point(mapping = aes(y = predict(weibull_model)), colour = 2) +
  geom_point(mapping = aes(y = predict(glm_model, type = "response")), colour = 3)

time_seq <- seq(0, max(dframe$time_sleeping), length.out = 100)

# Predict the linear predictor (eta) from the model for each observation
linear_predictor <- predict(weibull_model, type = "lp")

# Extract shape and scale parameters
shape <- 1 / weibull_model$scale  # Weibull shape parameter
scale <- exp(linear_predictor)    # Scale parameter for each observation

# Calculate survival probabilities for each observation across all time points
# Initialize a matrix to store survival probabilities for all observations
surv_probs <- matrix(NA, nrow = length(time_seq), ncol = length(scale))

# Loop over each time point and calculate the survival probabilities
for (i in 1:length(time_seq)) {
  surv_probs[i, ] <- exp(-(time_seq[i] / scale)^shape)
}

# Plotting the survival curves for a few observations (e.g., the first 5 observations)
matplot(time_seq, surv_probs[, 1:5], type = "l", lty = 1, col = 1:5,
        xlab = "Time", ylab = "Survival Probability",
        main = "Survival Curves for Selected Observations (Weibull Model)")
legend("topright", legend = paste("Obs", 1:5), col = 1:5, lty = 1)


# Calculate the log-likelihoods
glm_df <- glm_model$df.residual
weibull_df <- weibull_model$df.residual
loglik_glm <- logLik(glm_model)[1]  # Log-likelihood from the GLM model
loglik_weibull <- weibull_model$loglik[2]  # Log-likelihood from the Weibull model
deviance_glm <- -2 * loglik_glm
deviance_weibull <- -2 * loglik_weibull

# Calculate the likelihood ratio test statistic (subtracting GLM from Weibull since Weibull is better)
lrt_stat <- 2 * (loglik_weibull - loglik_glm)

# Degrees of freedom difference
df_diff <- abs(glm_df - weibull_df)

# Calculate the p-value
p_value <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)

# Display the results
cat("Likelihood Ratio Test Statistic:", lrt_stat, "\n")
cat("p-value:", p_value, "\n")


# Extract fitted values
fitted_values <- predict(weibull_model, type = "response")

# Extract residuals
residuals_weibull <- residuals(weibull_model, type = "deviance")

plot(fitted_values, residuals_weibull,
     xlab = "Fitted Values",
     ylab = "Deviance Residuals",
     main = "Residuals vs Fitted Values (Weibull Model)",
     pch = 20)
spline_fit <- smooth.spline(fitted_values, residuals_weibull)
# Add the fitted spline line to the plot
lines(spline_fit, col = "red", lwd = 1, lty=5)


qqnorm(residuals_weibull, main = "Q-Q Plot of Residuals (Weibull Model)")
qqline(residuals_weibull, col = "red")

sqrt_abs_residuals <- sqrt(abs(residuals_weibull))
plot(fitted_values, sqrt_abs_residuals,
     xlab = "Fitted Values",
     ylab = "Sqrt(|Residuals|)",
     main = "Scale-Location Plot (Weibull Model)",
     pch = 20)
spline_fit <- smooth.spline(fitted_values, sqrt_abs_residuals)
# Add the fitted spline line to the plot
lines(spline_fit, col = "red", lwd = 1, lty=5)

