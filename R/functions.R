# Functions 

check_brms <- function(model,             # brms model
                        seed = 10,
                        integer = FALSE,   # integer response? (TRUE/FALSE)
                        plot = TRUE,       # make plot?
                        ...                # further arguments for DHARMa::plotResiduals 
 ) {
   mdata <- brms::standata(model)
   if (!"Y" %in% names(mdata)) {
     stop("Cannot extract the required information from this brms model")
   }
   options(mc.cores = 1)
   on.exit(options(mc.cores = parallel::detectCores()))
   ndraws <- nrow(brms::as_draws_df(model))
   manual_preds_brms <- matrix(0, ndraws, nrow(model$data))
   random_terms <- insight::find_random(
     model, split_nested = TRUE, flatten = TRUE
   )
   new_data <- model$data |>
     dplyr::mutate(across(
       all_of(random_terms), \(x)paste0("NEW_", x) |> as.factor()
     ))
   set.seed(seed)
   brms_sims <- brms::posterior_predict(
   model, re_formula = NULL, newdata = new_data,
     allow_new_levels = TRUE, sample_new_levels = "gaussian"
   ) |>
     t()
   fitted_median_brms <- apply(brms_sims, 1, median) 
   dharma.obj <- DHARMa::createDHARMa(
     simulatedResponse = brms_sims,
     observedResponse = mdata$Y, 
     fittedPredictedResponse = fitted_median_brms,
     integerResponse = integer)
   if (plot) {
     plot(dharma.obj, ...)
   }
   invisible(dharma.obj)

 }
