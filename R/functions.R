# Packages and functions 

# ######################################## GLLVMs

coef.plot <- function(m_, X_coef) {
coefficients <- coef(m_)

dat_plot <-  confint(m_) %>%
  data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  filter(str_starts(ID, X_coef))  %>%
  cbind(as.vector(coefficients$Xcoef[,1])) %>%
  mutate(class = str_split(ID, pattern = "[:]") %>% map_chr(.,last)) %>%
  rename(
    lower = X2.5..,
    upper = X97.5..,
    mean = `as.vector(coefficients$Xcoef[, 1])`
  ) %>%
  dplyr::select(class, mean, lower, upper) %>%
  mutate(Sig = case_when((lower <= 0 & upper <= 0) | (lower > 0 & upper > 0) ~ 1 , TRUE ~ 0)) %>%
  arrange(desc(mean)) %>% 
  mutate(class = fct_inorder(class))  

 p <- ggplot(dat_plot, aes(x = mean, y = class, col = as.factor(Sig))) + 
   geom_point(size = 2, position=position_dodge(width=0.5)) + 
   geom_errorbar(aes(xmin = lower, xmax = upper),
                     position=position_dodge(width=0.5),
                 width = 0.2) + 
   theme_bw() + 
   xlab("Differences of predicted cover (in logit scale)") +
   ylab("Class") +
   geom_vline(xintercept = 0, linetype = "dashed") +
   scale_color_manual(values = c("1" = "black", "0" = "gray66")) +
   theme(legend.position = "none", axis.text.x = element_text(size=11),
         axis.text.y = element_text(size=11),axis.title.y=element_text(size=12),
         axis.title.x=element_text(size=12),
         strip.text = element_text(size = 11),
         strip.background = element_rect(fill = 'white'))

return(p)
}

plot_residuals <- function(m_, X, X_coef, name){

filename <- file.path("appendix/figures", paste0(name, ".png"))

png(filename, width = 1000, height = 400)
## Fitted values 
fitted_vals <- predict(m_, X, )

## residuals 
res <- residuals(m_, type = "response") 

par(mfrow = c(1, 3))
plot(fitted_vals, res$residuals, 
     xlab = "Fitted values", ylab = "Residuals", 
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)
hist(res$residuals, main = "Histogram of Residuals", xlab = "Residuals")
boxplot(res$residuals ~ X_coef, main = "Residuals by Xcoef", xlab="", ylab = "residuals")

# Close graphics device
dev.off()
  
# Return invisible result
invisible(filename)

# Capture and return the plot
#plot_obj <- recordPlot()
#return(plot_obj)
}


goodness_of_fit <- function(m_, y){

## Fitted values 
fitted_vals <- predict(m_)

df_true <- as.data.frame(y) %>%
  mutate(obs_id = row_number()) %>%
  pivot_longer(-obs_id, names_to = "comm", values_to = "observed")

df_fit <- as.data.frame(fitted_vals) %>%
  mutate(obs_id = row_number()) %>%
  pivot_longer(-obs_id, names_to = "comm", values_to = "fitted") %>%
  mutate(fitted = plogis(fitted))

# Join observed and fitted
df_plot <- left_join(df_true, df_fit, by = c("obs_id", "comm"))
return(df_plot)
}


plot_fit <- function(df_plot){
# Plot fitted vs observed
pfit <- ggplot(df_plot, aes(x = observed, y = fitted)) +
  geom_point(alpha = 0.4) +
 # geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  facet_wrap(~comm, scales = "free") +
  labs(
    title = "",
    x = "Observed",
    y = "Predicted"
  ) +
  theme_minimal()
return(pfit)
}
