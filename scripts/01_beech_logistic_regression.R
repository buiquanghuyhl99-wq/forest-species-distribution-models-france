# ==============================================================================
# Beech Species Distribution Modelling with Logistic Regression
# ==============================================================================
# Author: Huy Quang Bui
#
# Purpose:
#   Build, screen, validate, and export a logistic-regression species
#   distribution model for beech presence/absence data.
#
# Notes:
#   - Raw data are intentionally excluded from this public repository.
#   - Update the CONFIGURATION section before running the script.
#   - This refactored workflow is adapted from an earlier teaching script by
#     C. Piedallu (2024). See the repository README for acknowledgements.
# ==============================================================================


# 1. CONFIGURATION -------------------------------------------------------------

input_file <- file.path("data", "database_ed.csv")
output_dir <- "outputs"

# Binary response variable: presence / absence
response_variable <- "qupe"

# Columns before this position are treated as identifiers or metadata.
# Explanatory variables begin at this column index.
first_predictor_column <- 6

# Reproducible calibration-validation split
random_seed <- 2024
calibration_fraction <- 0.70

# Final model selected after forward screening.
# Terms may include transformations such as I(variable^2).
final_model_terms <- c(
  "Tmean_au",
  "I(Tmean_au^2)",
  "TW",
  "PW",
  "CN"
)

# Variables for exported partial-response curves.
response_curve_variables <- c("Tmean_au", "TW", "PW", "CN")


# 2. PACKAGE SETUP --------------------------------------------------------------

required_packages <- c("dplyr", "pROC")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Install the missing packages before running this script: ",
      paste(missing_packages, collapse = ", "),
      "\nRun: install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    )
  )
}

library(dplyr)

dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)


# 3. HELPER FUNCTIONS -----------------------------------------------------------

validate_input_data <- function(data, response_name, predictor_start) {
  if (!response_name %in% names(data)) {
    stop("Response variable not found: ", response_name)
  }

  if (predictor_start > ncol(data)) {
    stop(
      "first_predictor_column is larger than the number of columns in the dataset."
    )
  }

  response_values <- unique(stats::na.omit(data[[response_name]]))

  if (!all(response_values %in% c(0, 1, "0", "1"))) {
    stop(
      "The response variable must contain only binary values: 0 and 1."
    )
  }
}


coerce_response_to_numeric <- function(data, response_name) {
  data[[response_name]] <- as.numeric(as.character(data[[response_name]]))

  if (any(!stats::na.omit(data[[response_name]]) %in% c(0, 1))) {
    stop("Response conversion failed: expected values 0 and 1.")
  }

  data
}


calculate_deviance_explained <- function(model) {
  if (is.null(model$null.deviance) || model$null.deviance == 0) {
    return(NA_real_)
  }

  (model$null.deviance - model$deviance) / model$null.deviance
}


calculate_auc <- function(model, data, response_name) {
  probabilities <- stats::predict(model, newdata = data, type = "response")
  roc_object <- pROC::roc(
    response = data[[response_name]],
    predictor = probabilities,
    quiet = TRUE
  )

  as.numeric(pROC::auc(roc_object))
}


fit_candidate_model <- function(
  data,
  response_name,
  candidate_name,
  selected_terms = character(0),
  use_quadratic = FALSE
) {
  candidate_terms <- candidate_name

  if (use_quadratic && is.numeric(data[[candidate_name]])) {
    candidate_terms <- c(candidate_terms, paste0("I(", candidate_name, "^2)"))
  }

  model_terms <- unique(c(selected_terms, candidate_terms))
  model_formula <- stats::reformulate(model_terms, response = response_name)

  stats::glm(
    formula = model_formula,
    family = stats::binomial(),
    data = data,
    na.action = stats::na.exclude
  )
}


screen_predictors <- function(
  data,
  response_name,
  predictor_names,
  selected_terms = character(0),
  step_number = 1
) {
  baseline_formula <- if (length(selected_terms) == 0) {
    stats::reformulate("1", response = response_name)
  } else {
    stats::reformulate(selected_terms, response = response_name)
  }

  baseline_model <- stats::glm(
    formula = baseline_formula,
    family = stats::binomial(),
    data = data,
    na.action = stats::na.exclude
  )

  results <- lapply(predictor_names, function(candidate_name) {
    if (candidate_name %in% selected_terms) {
      return(NULL)
    }

    linear_model <- fit_candidate_model(
      data = data,
      response_name = response_name,
      candidate_name = candidate_name,
      selected_terms = selected_terms,
      use_quadratic = FALSE
    )

    quadratic_model <- if (is.numeric(data[[candidate_name]])) {
      fit_candidate_model(
        data = data,
        response_name = response_name,
        candidate_name = candidate_name,
        selected_terms = selected_terms,
        use_quadratic = TRUE
      )
    } else {
      linear_model
    }

    comparison <- stats::anova(
      baseline_model,
      quadratic_model,
      test = "Chisq"
    )

    tibble::tibble(
      step = step_number,
      candidate = candidate_name,
      variable_type = class(data[[candidate_name]])[1],
      deviance_explained_linear = calculate_deviance_explained(linear_model),
      deviance_explained_quadratic = calculate_deviance_explained(quadratic_model),
      auc_linear = calculate_auc(linear_model, data, response_name),
      auc_quadratic = calculate_auc(quadratic_model, data, response_name),
      p_value_vs_previous_model = comparison$`Pr(>Chi)`[2],
      rows_used_linear = stats::nobs(linear_model),
      rows_used_quadratic = stats::nobs(quadratic_model)
    )
  })

  dplyr::bind_rows(results) %>%
    arrange(desc(deviance_explained_quadratic), desc(auc_quadratic))
}


get_reference_value <- function(variable) {
  if (is.numeric(variable)) {
    return(stats::median(variable, na.rm = TRUE))
  }

  most_common <- names(sort(table(variable), decreasing = TRUE))[1]
  factor(most_common, levels = levels(variable))
}


export_partial_response_curve <- function(
  model,
  data,
  variable_name,
  figures_dir,
  tables_dir
) {
  if (!variable_name %in% names(data)) {
    warning("Response-curve variable not found: ", variable_name)
    return(invisible(NULL))
  }

  if (!is.numeric(data[[variable_name]])) {
    warning("Response curves are exported only for numeric variables: ", variable_name)
    return(invisible(NULL))
  }

  model_variables <- all.vars(stats::formula(model))
  predictor_variables <- setdiff(model_variables, response_variable)

  reference_data <- lapply(
    data[predictor_variables],
    get_reference_value
  ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  gradient <- seq(
    min(data[[variable_name]], na.rm = TRUE),
    max(data[[variable_name]], na.rm = TRUE),
    length.out = 150
  )

  prediction_data <- reference_data[rep(1, length(gradient)), , drop = FALSE]
  prediction_data[[variable_name]] <- gradient

  predicted_probability <- stats::predict(
    model,
    newdata = prediction_data,
    type = "response"
  )

  curve_table <- data.frame(
    variable = variable_name,
    value = gradient,
    predicted_probability = predicted_probability
  )

  utils::write.csv(
    curve_table,
    file.path(tables_dir, paste0("partial_response_", variable_name, ".csv")),
    row.names = FALSE
  )

  grDevices::png(
    filename = file.path(figures_dir, paste0("partial_response_", variable_name, ".png")),
    width = 1000,
    height = 700,
    res = 120
  )

  graphics::plot(
    gradient,
    predicted_probability,
    type = "l",
    xlab = variable_name,
    ylab = "Predicted probability of presence",
    main = paste("Partial response curve:", variable_name)
  )

  grDevices::dev.off()
}


# 4. LOAD AND PREPARE DATA ------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPlace your private CSV file in the data/ folder."
  )
}

raw_data <- utils::read.csv2(
  file = input_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

validate_input_data(
  data = raw_data,
  response_name = response_variable,
  predictor_start = first_predictor_column
)

model_data <- coerce_response_to_numeric(raw_data, response_variable)

predictor_names <- names(model_data)[first_predictor_column:ncol(model_data)]

# Convert character explanatory variables to factors.
model_data[predictor_names] <- lapply(
  model_data[predictor_names],
  function(variable) {
    if (is.character(variable)) factor(variable) else variable
  }
)

# Use complete rows across the response and explanatory variables so that
# candidate models are compared on the same observations.
analysis_columns <- unique(c(response_variable, predictor_names))
model_data <- model_data[stats::complete.cases(model_data[analysis_columns]), ]

if (nrow(model_data) == 0) {
  stop("No complete rows remain after filtering missing values.")
}

set.seed(random_seed)
calibration_indices <- sample(
  seq_len(nrow(model_data)),
  size = floor(calibration_fraction * nrow(model_data)),
  replace = FALSE
)

calibration_data <- model_data[calibration_indices, , drop = FALSE]
validation_data <- model_data[-calibration_indices, , drop = FALSE]

data_summary <- data.frame(
  metric = c(
    "rows_raw",
    "rows_complete",
    "rows_calibration",
    "rows_validation",
    "number_of_predictors"
  ),
  value = c(
    nrow(raw_data),
    nrow(model_data),
    nrow(calibration_data),
    nrow(validation_data),
    length(predictor_names)
  )
)

utils::write.csv(
  data_summary,
  file.path(output_dir, "tables", "data_summary.csv"),
  row.names = FALSE
)


# 5. EXPLORATORY CHECKS ---------------------------------------------------------

numeric_predictors <- predictor_names[
  vapply(model_data[predictor_names], is.numeric, FUN.VALUE = logical(1))
]

if (length(numeric_predictors) > 1) {
  correlation_matrix <- round(
    stats::cor(
      model_data[numeric_predictors],
      use = "pairwise.complete.obs",
      method = "pearson"
    ),
    digits = 3
  )

  utils::write.csv(
    correlation_matrix,
    file.path(output_dir, "tables", "predictor_correlation_matrix.csv")
  )
}


# 6. FORWARD-SCREENING TABLES ---------------------------------------------------

# These tables rank candidate predictors at each forward-selection stage.
# The final selected variables remain configurable near the top of this script.

screen_step_1 <- screen_predictors(
  data = calibration_data,
  response_name = response_variable,
  predictor_names = predictor_names,
  selected_terms = character(0),
  step_number = 1
)

utils::write.csv(
  screen_step_1,
  file.path(output_dir, "tables", "screening_step_1.csv"),
  row.names = FALSE
)

# Reproduce the intended selection sequence from the original analysis:
# Tmean_au (quadratic), TW, PW, and CN.
selection_sequence <- list(
  c("Tmean_au", "I(Tmean_au^2)"),
  c("Tmean_au", "I(Tmean_au^2)", "TW"),
  c("Tmean_au", "I(Tmean_au^2)", "TW", "PW")
)

for (step_index in seq_along(selection_sequence)) {
  screening_table <- screen_predictors(
    data = calibration_data,
    response_name = response_variable,
    predictor_names = predictor_names,
    selected_terms = selection_sequence[[step_index]],
    step_number = step_index + 2
  )

  utils::write.csv(
    screening_table,
    file.path(
      output_dir,
      "tables",
      paste0("screening_step_", step_index + 2, ".csv")
    ),
    row.names = FALSE
  )
}


# 7. FINAL MODEL ----------------------------------------------------------------

missing_final_variables <- setdiff(
  all.vars(stats::reformulate(final_model_terms)),
  names(calibration_data)
)

if (length(missing_final_variables) > 0) {
  stop(
    "Final-model variables are missing from the dataset: ",
    paste(missing_final_variables, collapse = ", ")
  )
}

final_formula <- stats::reformulate(
  final_model_terms,
  response = response_variable
)

final_model <- stats::glm(
  formula = final_formula,
  family = stats::binomial(),
  data = calibration_data,
  na.action = stats::na.exclude
)

capture.output(
  summary(final_model),
  file = file.path(output_dir, "tables", "final_model_summary.txt")
)

coefficient_table <- as.data.frame(summary(final_model)$coefficients)
coefficient_table$term <- rownames(coefficient_table)
rownames(coefficient_table) <- NULL
coefficient_table <- coefficient_table[
  ,
  c("term", setdiff(names(coefficient_table), "term"))
]

utils::write.csv(
  coefficient_table,
  file.path(output_dir, "tables", "final_model_coefficients.csv"),
  row.names = FALSE
)


# 8. VALIDATION -----------------------------------------------------------------

validation_probabilities <- stats::predict(
  final_model,
  newdata = validation_data,
  type = "response"
)

roc_curve <- pROC::roc(
  response = validation_data[[response_variable]],
  predictor = validation_probabilities,
  quiet = TRUE
)

validation_auc <- as.numeric(pROC::auc(roc_curve))

optimal_threshold <- pROC::coords(
  roc_curve,
  x = "best",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

validation_metrics <- data.frame(
  metric = c("auc", "threshold", "sensitivity", "specificity"),
  value = c(
    validation_auc,
    optimal_threshold$threshold,
    optimal_threshold$sensitivity,
    optimal_threshold$specificity
  )
)

utils::write.csv(
  validation_metrics,
  file.path(output_dir, "tables", "validation_metrics.csv"),
  row.names = FALSE
)

validation_export <- validation_data
validation_export$predicted_probability <- validation_probabilities

utils::write.csv(
  validation_export,
  file.path(output_dir, "tables", "validation_predictions.csv"),
  row.names = FALSE
)

grDevices::png(
  filename = file.path(output_dir, "figures", "validation_roc_curve.png"),
  width = 900,
  height = 700,
  res = 120
)

graphics::plot(
  roc_curve,
  main = paste0("Validation ROC curve — AUC: ", round(validation_auc, 3))
)

grDevices::dev.off()


# 9. PARTIAL-RESPONSE CURVES ----------------------------------------------------

for (variable_name in response_curve_variables) {
  export_partial_response_curve(
    model = final_model,
    data = calibration_data,
    variable_name = variable_name,
    figures_dir = file.path(output_dir, "figures"),
    tables_dir = file.path(output_dir, "tables")
  )
}


# 10. COMPLETION MESSAGE --------------------------------------------------------

message(
  "\nAnalysis complete.\n",
  "Author: Huy Quang Bui\n",
  "Calibration rows: ", nrow(calibration_data), "\n",
  "Validation rows: ", nrow(validation_data), "\n",
  "Validation AUC: ", round(validation_auc, 3), "\n",
  "Results saved in: ", output_dir
)
