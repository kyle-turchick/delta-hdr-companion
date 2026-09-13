# Exports the five included source-data CSVs for Figure 1b–d, Extended Data Figs. 1–2,
# and Extended Data Tables 1–2. Run separately from the repository root with:
#   Rscript --vanilla code/export_figure_table_source_data.R

find_source_data_project_root <- function(start = getwd(), marker = ".mi_hdr_project_root") {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, marker))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not find the repository root containing marker: ", marker, call. = FALSE)
}

source_data_filenames <- function() {
  c(
    observations = "observations.csv",
    fitted_mean_curves = "fitted_mean_curves.csv",
    model_comparisons = "model_comparisons.csv",
    model_terms = "model_terms.csv",
    figure1_illustrative_curves = "figure1_illustrative_curves.csv"
  )
}

source_data_analysis_metadata <- function() {
  tibble::tribble(
    ~analysis_id, ~analysis_display_order, ~dataset_id, ~analysis_role, ~adjustment_scenario, ~response_name, ~response_unit, ~delta_unit, ~information_criterion,
    "Foliage_delta_unadjusted_gaussian", 1L, "macarthur_1961", "primary", "unadjusted", "natural-log Shannon bird diversity", "dimensionless", "ft", "AICc",
    "Foliage_delta_mean_foliage_height_gaussian", 2L, "macarthur_1961", "sensitivity", "mean_foliage_height", "natural-log Shannon bird diversity", "dimensionless", "ft", "AICc",
    "Topography_delta_no_mean_covariates_gaussian", 3L, "allouche_2012", "primary", "no_mean_covariates", "breeding-bird species richness", "count", "m", "AIC",
    "Topography_delta_mean_elevation_linear_gaussian", 4L, "allouche_2012", "sensitivity", "mean_elevation_linear", "breeding-bird species richness", "count", "m", "AIC",
    "Topography_delta_mean_elevation_quadratic_gaussian", 5L, "allouche_2012", "sensitivity", "mean_elevation_quadratic", "breeding-bird species richness", "count", "m", "AIC"
  )
}

source_data_candidate_order <- function() {
  c(linear = 1L, logarithmic = 2L, quadratic = 3L)
}

source_data_candidate_sets <- function(hierarchy) {
  sets <- c(
    hierarchy$foliage_delta_models$sets,
    hierarchy$topography_delta_models$sets
  )
  required <- source_data_analysis_metadata()$analysis_id
  if (!all(required %in% names(sets))) {
    stop(
      "The analysis hierarchy is missing source-data candidate sets: ",
      paste(setdiff(required, names(sets)), collapse = ", "),
      call. = FALSE
    )
  }
  sets[required]
}

build_source_data_observations <- function(results, foliage_samples) {
  if (!identical(names(foliage_samples), LETTERS[1:13])) {
    stop("The foliage table must provide sites A-M in order.", call. = FALSE)
  }
  finite_counts <- vapply(
    foliage_samples,
    function(values) sum(is.finite(as.numeric(values))),
    integer(1)
  )
  if (any(finite_counts <= 0L)) {
    stop("Every foliage site must contain at least one finite retained value.", call. = FALSE)
  }

  macarthur <- results$macarthur$metrics %>%
    dplyr::arrange(match(site, LETTERS[1:13])) %>%
    dplyr::transmute(
      dataset_id = "macarthur_1961",
      observation_id = as.character(site),
      site_id = as.character(site),
      utm10 = NA_character_,
      source_fid = NA_character_,
      raster_filename = NA_character_,
      environmental_variable = "foliage height",
      environmental_mean = as.numeric(mu),
      environmental_mean_unit = "ft",
      environmental_sample_variance = as.numeric(sigma2),
      environmental_variance_unit = "ft^2",
      delta = as.numeric(delta),
      delta_unit = "ft",
      response_name = "natural-log Shannon bird diversity",
      response_value = as.numeric(bird_species_diversity),
      response_unit = "dimensionless",
      retained_value_count = as.integer(finite_counts[site])
    )

  catalonia <- results$allouche$allouche_model_df %>%
    dplyr::mutate(
      source_fid_sort = suppressWarnings(as.integer(sub("\\.tif$", "", as.character(fid))))
    )
  if (any(!is.finite(catalonia$source_fid_sort))) {
    stop("Every retained Catalonia row must have a numeric .tif lineage identifier.", call. = FALSE)
  }
  catalonia <- catalonia %>%
    dplyr::arrange(source_fid_sort) %>%
    dplyr::transmute(
      dataset_id = "allouche_2012",
      observation_id = as.character(UTM10),
      site_id = NA_character_,
      utm10 = as.character(UTM10),
      source_fid = as.character(fid),
      raster_filename = basename(as.character(fid)),
      environmental_variable = "elevation",
      environmental_mean = as.numeric(mu),
      environmental_mean_unit = "m",
      environmental_sample_variance = as.numeric(sigma2),
      environmental_variance_unit = "m^2",
      delta = as.numeric(delta),
      delta_unit = "m",
      response_name = "breeding-bird species richness",
      response_value = as.numeric(richness),
      response_unit = "count",
      retained_value_count = as.integer(n)
    )

  dplyr::bind_rows(macarthur, catalonia)
}

source_data_fixed_values <- function(candidate_set, adjustment_scenario) {
  has_mu <- "mu" %in% names(candidate_set$data)
  has_mu2 <- "mu2" %in% names(candidate_set$data)
  list(
    foliage = if (identical(adjustment_scenario, "mean_foliage_height") && has_mu) {
      mean(candidate_set$data$mu, na.rm = TRUE)
    } else {
      NA_real_
    },
    elevation = if (adjustment_scenario %in% c("mean_elevation_linear", "mean_elevation_quadratic") && has_mu) {
      mean(candidate_set$data$mu, na.rm = TRUE)
    } else {
      NA_real_
    },
    elevation_squared = if (identical(adjustment_scenario, "mean_elevation_quadratic") && has_mu2) {
      mean(candidate_set$data$mu2, na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

build_source_data_fitted_mean_curves <- function(candidate_sets) {
  metadata <- source_data_analysis_metadata()
  family_order <- source_data_candidate_order()
  primary_ids <- unname(figure1_primary_analysis_ids())

  dplyr::bind_rows(lapply(seq_len(nrow(metadata)), function(i) {
    meta <- metadata[i, , drop = FALSE]
    candidate_set <- candidate_sets[[meta$analysis_id]]
    selected_family <- figure1_select_plotted_family(candidate_set$summary)
    if (meta$analysis_id %in% primary_ids && !identical(selected_family, "logarithmic")) {
      stop("Both Figure 1 primary analyses must select the logarithmic family.", call. = FALSE)
    }
    fixed <- source_data_fixed_values(candidate_set, meta$adjustment_scenario)

    dplyr::bind_rows(lapply(names(family_order), function(model_family) {
      predictions <- candidate_set$predictions %>%
        dplyr::filter(as.character(model_family) == .env$model_family) %>%
        dplyr::arrange(x)
      if (nrow(predictions) != 200L) {
        stop("Each analysis-family prediction grid must contain 200 rows.", call. = FALSE)
      }

      include_interval <- meta$analysis_id %in% primary_ids && identical(model_family, selected_family)
      if (include_interval) {
        intervals <- figure1_prediction_confidence_band(candidate_set, model_family, level = 0.95)
        if (!identical(as.numeric(intervals$x), as.numeric(predictions$x))) {
          stop("Figure 1 confidence-band coordinates do not match stored predictions.", call. = FALSE)
        }
        confidence_lower <- as.numeric(intervals$conf_low)
        confidence_upper <- as.numeric(intervals$conf_high)
      } else {
        confidence_lower <- rep(NA_real_, nrow(predictions))
        confidence_upper <- rep(NA_real_, nrow(predictions))
      }

      tibble::tibble(
        analysis_id = meta$analysis_id,
        analysis_display_order = as.integer(meta$analysis_display_order),
        dataset_id = meta$dataset_id,
        analysis_role = meta$analysis_role,
        adjustment_scenario = meta$adjustment_scenario,
        model_family = model_family,
        candidate_display_order = as.integer(family_order[[model_family]]),
        prediction_index = seq_len(nrow(predictions)),
        delta = as.numeric(predictions$x),
        delta_unit = meta$delta_unit,
        fitted_mean_response = as.numeric(predictions$yhat),
        response_name = meta$response_name,
        response_unit = meta$response_unit,
        fixed_mean_foliage_height_ft = fixed$foliage,
        fixed_mean_elevation_m = fixed$elevation,
        fixed_mean_elevation_squared_m2 = fixed$elevation_squared,
        confidence_interval_included = rep(include_interval, nrow(predictions)),
        confidence_interval_type = if (include_interval) "pointwise_fitted_mean" else NA_character_,
        confidence_level = if (include_interval) 0.95 else NA_real_,
        confidence_lower = confidence_lower,
        confidence_upper = confidence_upper
      )
    }))
  }))
}

build_source_data_model_comparisons <- function(candidate_sets) {
  metadata <- source_data_analysis_metadata()
  family_order <- source_data_candidate_order()

  dplyr::bind_rows(lapply(seq_len(nrow(metadata)), function(i) {
    meta <- metadata[i, , drop = FALSE]
    summary_tbl <- candidate_sets[[meta$analysis_id]]$summary %>%
      dplyr::mutate(
        rank_value = as.integer(rank(IC, ties.method = "min")),
        best_value = IC == min(IC)
      )

    dplyr::bind_rows(lapply(names(family_order), function(model_family) {
      row <- summary_tbl %>%
        dplyr::filter(as.character(model_family) == .env$model_family)
      if (nrow(row) != 1L) {
        stop("Each analysis-family comparison must have exactly one row.", call. = FALSE)
      }
      is_quadratic <- identical(model_family, "quadratic")
      tibble::tibble(
        analysis_id = meta$analysis_id,
        analysis_display_order = as.integer(meta$analysis_display_order),
        dataset_id = meta$dataset_id,
        analysis_role = meta$analysis_role,
        adjustment_scenario = meta$adjustment_scenario,
        model_family = model_family,
        candidate_display_order = as.integer(family_order[[model_family]]),
        formula = as.character(row$formula),
        n = as.integer(row$n),
        k = as.integer(row$k),
        information_criterion = as.character(row$selection_criterion),
        information_criterion_value = as.numeric(row$IC),
        delta_ic = as.numeric(row$delta_IC),
        akaike_weight = as.numeric(row$Akaike_weight),
        rank = as.integer(row$rank_value),
        best_by_ic = as.logical(row$best_value),
        predictor_min = as.numeric(row$predictor_min),
        predictor_max = as.numeric(row$predictor_max),
        quadratic_vertex = if (is_quadratic) as.numeric(row$quadratic_vertex_x) else NA_real_,
        quadratic_vertex_inside_observed_range = if (is_quadratic) as.logical(row$quadratic_vertex_inside_observed_range) else NA,
        quadratic_high_side_declines = if (is_quadratic) as.logical(row$quadratic_high_side_declines) else NA,
        shape_supported_by_coefficients = as.logical(row$shape_supported_by_coefficients),
        shape_support_with_uncertainty = as.character(row$shape_support_with_uncertainty),
        shape_interpretation = as.character(row$shape_interpretation)
      )
    }))
  }))
}

build_source_data_model_terms <- function(candidate_sets) {
  metadata <- source_data_analysis_metadata()
  family_order <- source_data_candidate_order()

  dplyr::bind_rows(lapply(seq_len(nrow(metadata)), function(i) {
    meta <- metadata[i, , drop = FALSE]
    candidate_set <- candidate_sets[[meta$analysis_id]]

    dplyr::bind_rows(lapply(names(family_order), function(model_family) {
      term_rows <- candidate_set$terms %>%
        dplyr::filter(as.character(model_family) == .env$model_family)
      coefficient_names <- names(stats::coef(candidate_set$models[[model_family]]))
      if (!identical(as.character(term_rows$term), coefficient_names)) {
        stop("Exported model terms do not preserve the fitted model-matrix order.", call. = FALSE)
      }
      formula_row <- candidate_set$summary %>%
        dplyr::filter(as.character(model_family) == .env$model_family)

      tibble::tibble(
        analysis_id = meta$analysis_id,
        analysis_display_order = as.integer(meta$analysis_display_order),
        dataset_id = meta$dataset_id,
        analysis_role = meta$analysis_role,
        adjustment_scenario = meta$adjustment_scenario,
        model_family = model_family,
        candidate_display_order = as.integer(family_order[[model_family]]),
        formula = as.character(formula_row$formula),
        term_order = seq_len(nrow(term_rows)),
        term = as.character(term_rows$term),
        term_role = as.character(term_rows$term_role),
        estimate = as.numeric(term_rows$estimate),
        std_error = as.numeric(term_rows$std.error),
        statistic = as.numeric(term_rows$statistic),
        p_value = as.numeric(term_rows$p.value),
        conf_low = as.numeric(term_rows$conf.low),
        conf_high = as.numeric(term_rows$conf.high)
      )
    }))
  }))
}

source_data_expected_columns <- function() {
  list(
    observations = c(
      "dataset_id", "observation_id", "site_id", "utm10", "source_fid",
      "raster_filename", "environmental_variable", "environmental_mean",
      "environmental_mean_unit", "environmental_sample_variance",
      "environmental_variance_unit", "delta", "delta_unit", "response_name",
      "response_value", "response_unit", "retained_value_count"
    ),
    fitted_mean_curves = c(
      "analysis_id", "analysis_display_order", "dataset_id", "analysis_role",
      "adjustment_scenario", "model_family", "candidate_display_order",
      "prediction_index", "delta", "delta_unit", "fitted_mean_response",
      "response_name", "response_unit", "fixed_mean_foliage_height_ft",
      "fixed_mean_elevation_m", "fixed_mean_elevation_squared_m2",
      "confidence_interval_included", "confidence_interval_type",
      "confidence_level", "confidence_lower", "confidence_upper"
    ),
    model_comparisons = c(
      "analysis_id", "analysis_display_order", "dataset_id", "analysis_role",
      "adjustment_scenario", "model_family", "candidate_display_order", "formula",
      "n", "k", "information_criterion", "information_criterion_value", "delta_ic",
      "akaike_weight", "rank", "best_by_ic", "predictor_min", "predictor_max",
      "quadratic_vertex", "quadratic_vertex_inside_observed_range",
      "quadratic_high_side_declines", "shape_supported_by_coefficients",
      "shape_support_with_uncertainty", "shape_interpretation"
    ),
    model_terms = c(
      "analysis_id", "analysis_display_order", "dataset_id", "analysis_role",
      "adjustment_scenario", "model_family", "candidate_display_order", "formula",
      "term_order", "term", "term_role", "estimate", "std_error", "statistic",
      "p_value", "conf_low", "conf_high"
    ),
    figure1_illustrative_curves = c(
      "illustrative_family", "point_index", "plotting_x", "plotting_y"
    )
  )
}

source_data_column_types <- function() {
  list(
    observations = readr::cols(
      dataset_id = readr::col_character(), observation_id = readr::col_character(),
      site_id = readr::col_character(), utm10 = readr::col_character(),
      source_fid = readr::col_character(), raster_filename = readr::col_character(),
      environmental_variable = readr::col_character(), environmental_mean = readr::col_double(),
      environmental_mean_unit = readr::col_character(), environmental_sample_variance = readr::col_double(),
      environmental_variance_unit = readr::col_character(), delta = readr::col_double(),
      delta_unit = readr::col_character(), response_name = readr::col_character(),
      response_value = readr::col_double(), response_unit = readr::col_character(),
      retained_value_count = readr::col_integer()
    ),
    fitted_mean_curves = readr::cols(
      analysis_id = readr::col_character(), analysis_display_order = readr::col_integer(),
      dataset_id = readr::col_character(), analysis_role = readr::col_character(),
      adjustment_scenario = readr::col_character(), model_family = readr::col_character(),
      candidate_display_order = readr::col_integer(), prediction_index = readr::col_integer(),
      delta = readr::col_double(), delta_unit = readr::col_character(),
      fitted_mean_response = readr::col_double(), response_name = readr::col_character(),
      response_unit = readr::col_character(), fixed_mean_foliage_height_ft = readr::col_double(),
      fixed_mean_elevation_m = readr::col_double(), fixed_mean_elevation_squared_m2 = readr::col_double(),
      confidence_interval_included = readr::col_logical(), confidence_interval_type = readr::col_character(),
      confidence_level = readr::col_double(), confidence_lower = readr::col_double(),
      confidence_upper = readr::col_double()
    ),
    model_comparisons = readr::cols(
      analysis_id = readr::col_character(), analysis_display_order = readr::col_integer(),
      dataset_id = readr::col_character(), analysis_role = readr::col_character(),
      adjustment_scenario = readr::col_character(), model_family = readr::col_character(),
      candidate_display_order = readr::col_integer(), formula = readr::col_character(),
      n = readr::col_integer(), k = readr::col_integer(), information_criterion = readr::col_character(),
      information_criterion_value = readr::col_double(), delta_ic = readr::col_double(),
      akaike_weight = readr::col_double(), rank = readr::col_integer(),
      best_by_ic = readr::col_logical(), predictor_min = readr::col_double(),
      predictor_max = readr::col_double(), quadratic_vertex = readr::col_double(),
      quadratic_vertex_inside_observed_range = readr::col_logical(),
      quadratic_high_side_declines = readr::col_logical(),
      shape_supported_by_coefficients = readr::col_logical(),
      shape_support_with_uncertainty = readr::col_character(),
      shape_interpretation = readr::col_character()
    ),
    model_terms = readr::cols(
      analysis_id = readr::col_character(), analysis_display_order = readr::col_integer(),
      dataset_id = readr::col_character(), analysis_role = readr::col_character(),
      adjustment_scenario = readr::col_character(), model_family = readr::col_character(),
      candidate_display_order = readr::col_integer(), formula = readr::col_character(),
      term_order = readr::col_integer(), term = readr::col_character(),
      term_role = readr::col_character(), estimate = readr::col_double(),
      std_error = readr::col_double(), statistic = readr::col_double(),
      p_value = readr::col_double(), conf_low = readr::col_double(),
      conf_high = readr::col_double()
    ),
    figure1_illustrative_curves = readr::cols(
      illustrative_family = readr::col_character(), point_index = readr::col_integer(),
      plotting_x = readr::col_double(), plotting_y = readr::col_double()
    )
  )
}

source_data_assert_unique_key <- function(table, key, table_name) {
  if (any(vapply(table[key], anyNA, logical(1))) || any(duplicated(table[key]))) {
    stop(table_name, " does not have a complete unique key: ", paste(key, collapse = ", "), call. = FALSE)
  }
}

source_data_assert_value_quality <- function(table, table_name) {
  if (any(vapply(table, is.factor, logical(1)))) {
    stop(table_name, " contains a factor column.", call. = FALSE)
  }
  character_columns <- vapply(table, is.character, logical(1))
  for (column in names(table)[character_columns]) {
    values <- table[[column]]
    if (any(!is.na(values) & !nzchar(values))) {
      stop(table_name, " contains an empty string in ", column, ".", call. = FALSE)
    }
  }
  numeric_columns <- vapply(table, is.numeric, logical(1))
  for (column in names(table)[numeric_columns]) {
    values <- table[[column]]
    if (any(is.nan(values)) || any(!is.na(values) & !is.finite(values))) {
      stop(table_name, " contains NaN or an infinite value in ", column, ".", call. = FALSE)
    }
  }
}

source_data_equal_numeric <- function(x, y, tolerance = 1e-12) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance, check.attributes = FALSE))
}

source_data_max_abs_difference <- function(expected, observed) {
  if (length(expected) != length(observed) || !identical(is.na(expected), is.na(observed))) {
    return(Inf)
  }
  keep <- !is.na(expected)
  if (!any(keep)) return(0)
  max(abs(as.numeric(expected[keep]) - as.numeric(observed[keep])))
}

source_data_with_c_numeric_locale <- function(callback) {
  prior_locale <- Sys.getlocale("LC_NUMERIC")
  locale_result <- suppressWarnings(Sys.setlocale("LC_NUMERIC", "C"))
  if (is.na(locale_result) || !identical(Sys.localeconv()[["decimal_point"]], ".")) {
    stop("Could not establish the C numeric locale for deterministic CSV serialization.", call. = FALSE)
  }
  on.exit({
    restore_result <- suppressWarnings(Sys.setlocale("LC_NUMERIC", prior_locale))
    if (is.na(restore_result) || !identical(Sys.getlocale("LC_NUMERIC"), prior_locale)) {
      stop("Could not restore LC_NUMERIC after deterministic CSV serialization.", call. = FALSE)
    }
  }, add = TRUE)
  callback()
}

source_data_double_decimal_tokens <- function(values) {
  if (!is.double(values)) {
    stop("Only double columns may use the deterministic decimal formatter.", call. = FALSE)
  }
  if (any(is.nan(values)) || any(!is.na(values) & !is.finite(values))) {
    stop("Cannot serialize a NaN or infinite double value.", call. = FALSE)
  }
  tokens <- rep(NA_character_, length(values))
  finite <- !is.na(values)
  if (any(finite)) {
    tokens[finite] <- sprintf("%.17g", values[finite])
    tokens[finite] <- sub("E", "e", tokens[finite], fixed = TRUE)
    tokens[finite] <- sub(
      "e([+-])0+([0-9]+)$",
      "e\\1\\2",
      tokens[finite],
      perl = TRUE
    )
    token_pattern <- "^-?(?:0|[1-9][0-9]*)(?:\\.[0-9]+)?(?:e[+-](?:0|[1-9][0-9]*))?$"
    if (any(!grepl(token_pattern, tokens[finite], perl = TRUE))) {
      stop("The deterministic double formatter produced a non-C decimal token.", call. = FALSE)
    }
  }
  tokens
}

source_data_serialization_tables <- function(tables) {
  source_data_with_c_numeric_locale(function() {
    lapply(tables, function(table) {
      serialized <- table
      double_columns <- names(table)[vapply(table, is.double, logical(1))]
      for (column in double_columns) {
        serialized[[column]] <- source_data_double_decimal_tokens(table[[column]])
      }
      serialized
    })
  })
}

source_data_table_key_columns <- function(table_name) {
  switch(
    table_name,
    observations = c("dataset_id", "observation_id"),
    fitted_mean_curves = c("analysis_id", "model_family", "prediction_index"),
    model_comparisons = c("analysis_id", "model_family"),
    model_terms = c("analysis_id", "model_family", "term_order"),
    figure1_illustrative_curves = c("illustrative_family", "point_index"),
    character(0)
  )
}

source_data_row_key <- function(table, table_name, row_number) {
  key_columns <- source_data_table_key_columns(table_name)
  if (length(key_columns) == 0L ||
      !all(key_columns %in% names(table)) ||
      is.na(row_number) || row_number < 1L || row_number > nrow(table)) {
    return("<unavailable>")
  }
  key_values <- vapply(key_columns, function(column) {
    value <- table[[column]][row_number]
    if (is.na(value)) "NA" else as.character(value)
  }, character(1))
  paste(paste0(key_columns, "=", key_values), collapse = ", ")
}

source_data_ordered_double_bytes <- function(value) {
  if (isTRUE(value == 0)) value <- 0
  bytes <- as.integer(writeBin(as.double(value), raw(), size = 8L, endian = "big"))
  if (bitwAnd(bytes[1], 128L) != 0L) {
    incremented <- source_data_increment_bytes(255L - bytes)
    if (is.null(incremented)) {
      stop("Could not construct the ordered representation of a negative double.", call. = FALSE)
    }
    incremented
  } else {
    bytes[1] <- bitwXor(bytes[1], 128L)
    bytes
  }
}

source_data_increment_bytes <- function(bytes) {
  out <- as.integer(bytes)
  for (i in rev(seq_along(out))) {
    if (out[i] < 255L) {
      out[i] <- out[i] + 1L
      return(out)
    }
    out[i] <- 0L
  }
  NULL
}

source_data_adjacent_double_status <- function(expected, observed) {
  if (!isTRUE(.Machine$double.base == 2) ||
      !isTRUE(.Machine$double.digits == 53) ||
      length(expected) != 1L || length(observed) != 1L ||
      is.na(expected) || is.na(observed) ||
      !is.finite(expected) || !is.finite(observed) ||
      isTRUE(expected == observed) ||
      identical(expected, observed, num.eq = FALSE)) {
    return("not applicable")
  }
  expected_bytes <- source_data_ordered_double_bytes(expected)
  observed_bytes <- source_data_ordered_double_bytes(observed)
  first_difference <- which(expected_bytes != observed_bytes)[1]
  if (is.na(first_difference)) return("false")
  expected_is_lower <- expected_bytes[first_difference] < observed_bytes[first_difference]
  lower <- if (expected_is_lower) expected_bytes else observed_bytes
  higher <- if (expected_is_lower) observed_bytes else expected_bytes
  incremented <- source_data_increment_bytes(lower)
  if (!is.null(incremented) && identical(incremented, higher)) "true" else "false"
}

source_data_double_compatibility <- function(
    expected, observed, relative_threshold = 1e-12) {
  if (length(expected) != length(observed)) {
    return(list(
      finite_value_count = length(expected),
      exact_value_count = 0L,
      invalid_rows = 1L,
      failure_reasons = "vector length changed",
      absolute_differences = rep(Inf, length(expected)),
      relative_differences = rep(Inf, length(expected)),
      max_abs_difference = Inf,
      max_relative_difference = Inf
    ))
  }
  missingness_mismatch <- which(is.na(expected) != is.na(observed))
  present <- !is.na(expected) & !is.na(observed)
  nonfinite_rows <- which(present & (!is.finite(expected) | !is.finite(observed)))
  finite <- present & is.finite(expected) & is.finite(observed)
  zero_rows <- which(finite & expected == 0)
  zero_mismatch <- zero_rows[observed[zero_rows] != 0]
  nonzero_rows <- which(finite & expected != 0)
  absolute_differences <- rep(NA_real_, length(expected))
  absolute_differences[finite] <- abs(observed[finite] - expected[finite])
  relative_differences <- rep(NA_real_, length(expected))
  relative_differences[nonzero_rows] <-
    absolute_differences[nonzero_rows] / abs(expected[nonzero_rows])
  relative_failure <- nonzero_rows[
    relative_differences[nonzero_rows] > relative_threshold
  ]
  invalid_rows <- unique(c(
    missingness_mismatch,
    nonfinite_rows,
    zero_mismatch,
    relative_failure
  ))
  failure_reason <- rep(NA_character_, length(expected))
  failure_reason[missingness_mismatch] <- "missingness changed"
  failure_reason[nonfinite_rows] <- "parsed value is not finite"
  failure_reason[zero_mismatch] <- "zero value changed"
  failure_reason[relative_failure] <- paste0(
    "relative difference exceeds ",
    source_data_full_precision_double(relative_threshold)
  )
  max_abs_difference <- if (length(missingness_mismatch) > 0L ||
      length(nonfinite_rows) > 0L) {
    Inf
  } else if (any(finite)) {
    max(absolute_differences[finite])
  } else {
    0
  }
  max_relative_difference <- if (length(missingness_mismatch) > 0L ||
      length(nonfinite_rows) > 0L) {
    Inf
  } else if (length(nonzero_rows) > 0L) {
    max(relative_differences[nonzero_rows])
  } else {
    0
  }
  list(
    finite_value_count = sum(!is.na(expected)),
    exact_value_count = sum(finite & observed == expected),
    invalid_rows = invalid_rows,
    failure_reasons = failure_reason[invalid_rows],
    absolute_differences = absolute_differences,
    relative_differences = relative_differences,
    max_abs_difference = max_abs_difference,
    max_relative_difference = max_relative_difference
  )
}

source_data_full_precision_double <- function(value) {
  if (length(value) != 1L || is.na(value)) return("NA")
  source_data_with_c_numeric_locale(function() sprintf("%.17g", value))
}

source_data_diagnostic_mismatch_row <- function(expected, observed) {
  if (length(expected) != length(observed)) return(NA_integer_)
  missingness_mismatch <- which(is.na(expected) != is.na(observed))
  if (length(missingness_mismatch) > 0L) return(missingness_mismatch[1])
  finite <- !is.na(expected)
  differences <- rep(NA_real_, length(expected))
  differences[finite] <- abs(expected[finite] - observed[finite])
  maximum <- if (any(finite)) max(differences[finite]) else 0
  if (is.finite(maximum) && maximum > 0) {
    return(which(differences == maximum)[1])
  }
  mismatch <- vapply(seq_along(expected), function(i) {
    !identical(expected[[i]], observed[[i]], num.eq = FALSE)
  }, logical(1))
  which(mismatch)[1]
}

source_data_round_trip_mismatch_message <- function(
    phase, table_name, column_name, expected_table, observed_table,
    serialized_table, compatibility, failure_index) {
  row_number <- compatibility$invalid_rows[[failure_index]]
  if (is.na(row_number)) row_number <- 1L
  serialized_token <- if (!is.null(serialized_table) &&
      column_name %in% names(serialized_table) &&
      row_number <= nrow(serialized_table)) {
    serialized_table[[column_name]][row_number]
  } else {
    NA_character_
  }
  if (is.na(serialized_token)) serialized_token <- "NA"
  expected_value <- expected_table[[column_name]][[row_number]]
  observed_value <- observed_table[[column_name]][[row_number]]
  row_difference <- compatibility$absolute_differences[[row_number]]
  row_relative_difference <- compatibility$relative_differences[[row_number]]
  paste0(
    phase, ": ", table_name, ".", column_name,
    " failed typed CSV parser compatibility; row = ", row_number,
    "; key = ", source_data_row_key(expected_table, table_name, row_number),
    "; serialized token = ", serialized_token,
    "; expected = ", source_data_full_precision_double(expected_value),
    "; reread = ", source_data_full_precision_double(observed_value),
    "; row absolute difference = ",
    source_data_full_precision_double(row_difference),
    "; row relative difference = ",
    source_data_full_precision_double(row_relative_difference),
    "; reason = ", compatibility$failure_reasons[[failure_index]], "."
  )
}

source_data_compare_table_sets <- function(
    expected, observed, phase, serialized_tables = NULL) {
  if (!identical(names(expected), names(observed))) {
    stop(phase, ": source-data table names or order changed.", call. = FALSE)
  }
  if (!is.null(serialized_tables) && !identical(names(expected), names(serialized_tables))) {
    stop(phase, ": serialization table names or order changed.", call. = FALSE)
  }

  diagnostics <- list()
  failures <- character(0)
  diagnostic_index <- 0L
  for (table_name in names(expected)) {
    expected_table <- expected[[table_name]]
    observed_table <- observed[[table_name]]
    if (!identical(names(expected_table), names(observed_table)) ||
        nrow(expected_table) != nrow(observed_table)) {
      stop(phase, ": ", table_name, " changed schema, column order, or row count.", call. = FALSE)
    }
    expected_types <- vapply(expected_table, typeof, character(1))
    observed_types <- vapply(observed_table, typeof, character(1))
    if (!identical(expected_types, observed_types)) {
      stop(phase, ": ", table_name, " changed one or more column types.", call. = FALSE)
    }

    for (column in names(expected_table)) {
      expected_values <- expected_table[[column]]
      observed_values <- observed_table[[column]]
      if (is.double(expected_values)) {
        compatibility <- source_data_double_compatibility(expected_values, observed_values)
        if (length(compatibility$invalid_rows) > 0L) {
          failures <- c(failures, vapply(
            seq_along(compatibility$invalid_rows),
            function(failure_index) source_data_round_trip_mismatch_message(
              phase = phase,
              table_name = table_name,
              column_name = column,
              expected_table = expected_table,
              observed_table = observed_table,
              serialized_table = if (is.null(serialized_tables)) NULL else serialized_tables[[table_name]],
              compatibility = compatibility,
              failure_index = failure_index
            ),
            character(1)
          ))
        }
        diagnostic_index <- diagnostic_index + 1L
        diagnostics[[diagnostic_index]] <- tibble::tibble(
          phase = phase,
          table_name = table_name,
          column_name = column,
          finite_value_count = compatibility$finite_value_count,
          exact_value_count = compatibility$exact_value_count,
          max_abs_difference = compatibility$max_abs_difference,
          max_relative_difference = compatibility$max_relative_difference,
          failure_count = length(compatibility$invalid_rows)
        )
      } else if (!identical(expected_values, observed_values)) {
        stop(
          phase, ": ", table_name, ".", column,
          " changed; non-double columns require exact equality.",
          call. = FALSE
        )
      }
    }
  }
  diagnostics <- dplyr::bind_rows(diagnostics)
  if (length(failures) > 0L) {
    stop(
      paste(
        c(
          paste0(
            phase, ": ", length(failures),
            " double value(s) failed the uniform 1e-12 relative parser-compatibility rule."
          ),
          failures
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  diagnostics
}

source_data_report_round_trip_diagnostics <- function(diagnostics) {
  if (is.null(diagnostics) || nrow(diagnostics) == 0L) {
    stop("No double-column round-trip diagnostics were produced.", call. = FALSE)
  }
  if (any(!is.finite(diagnostics$max_abs_difference)) ||
      any(!is.finite(diagnostics$max_relative_difference)) ||
      any(diagnostics$max_relative_difference > 1e-12) ||
      any(diagnostics$failure_count != 0L)) {
    stop("Round-trip diagnostics exceed the uniform 1e-12 relative parser-compatibility rule.", call. = FALSE)
  }
  for (table_name in unique(diagnostics$table_name)) {
    rows <- diagnostics[diagnostics$table_name == table_name, , drop = FALSE]
    message(
      unique(rows$phase), ": ", table_name,
      " finite double values = ", sum(rows$finite_value_count),
      "; exact values = ", sum(rows$exact_value_count),
      "; maximum absolute difference = ",
      source_data_full_precision_double(max(rows$max_abs_difference)),
      "; maximum relative difference = ",
      source_data_full_precision_double(max(rows$max_relative_difference)), "."
    )
  }
  invisible(diagnostics)
}

source_data_assert_reference_compatibility <- function(
    observed, expected, round_trip_diagnostics, phase, table_name, column_name) {
  compatibility <- source_data_double_compatibility(expected, observed)
  if (is.null(round_trip_diagnostics)) {
    compatible <- length(compatibility$invalid_rows) == 0L &&
      compatibility$max_abs_difference == 0
  } else {
    diagnostic <- round_trip_diagnostics[
      round_trip_diagnostics$table_name == table_name &
        round_trip_diagnostics$column_name == column_name,
      ,
      drop = FALSE
    ]
    compatible <- nrow(diagnostic) == 1L &&
      length(compatibility$invalid_rows) == 0L &&
      diagnostic$failure_count == 0L
  }
  if (!compatible) {
    stop(
      phase, ": ", table_name, ".", column_name,
      " differs from its authoritative source beyond the applicable typed CSV rule; ",
      "maximum absolute difference = ",
      source_data_full_precision_double(compatibility$max_abs_difference),
      "; maximum relative difference = ",
      source_data_full_precision_double(compatibility$max_relative_difference), ".",
      call. = FALSE
    )
  }
}

validate_source_data_tables <- function(
    tables, candidate_sets, phase = "in-memory construction",
    round_trip_diagnostics = NULL) {
  expected_columns <- source_data_expected_columns()
  expected_rows <- c(
    observations = 298L,
    fitted_mean_curves = 3000L,
    model_comparisons = 15L,
    model_terms = 47L,
    figure1_illustrative_curves = 603L
  )
  if (!identical(names(tables), names(expected_rows))) {
    stop("The source-data table set is incomplete or out of order.", call. = FALSE)
  }
  for (name in names(tables)) {
    if (!identical(names(tables[[name]]), expected_columns[[name]])) {
      stop(name, " has an unexpected column schema or order.", call. = FALSE)
    }
    if (nrow(tables[[name]]) != expected_rows[[name]]) {
      stop(name, " has an unexpected row count.", call. = FALSE)
    }
    source_data_assert_value_quality(tables[[name]], name)
  }

  observations <- tables$observations
  source_data_assert_unique_key(observations, c("dataset_id", "observation_id"), "observations")
  observation_counts <- table(observations$dataset_id)
  if (!identical(
    names(observation_counts),
    c("allouche_2012", "macarthur_1961")
  ) || !identical(as.integer(observation_counts), c(285L, 13L))) {
    stop("observations must contain 13 MacArthur and 285 Catalonia rows.", call. = FALSE)
  }
  macarthur <- observations[observations$dataset_id == "macarthur_1961", , drop = FALSE]
  catalonia <- observations[observations$dataset_id == "allouche_2012", , drop = FALSE]
  if (!identical(macarthur$observation_id, LETTERS[1:13]) || !identical(macarthur$site_id, LETTERS[1:13])) {
    stop("MacArthur observation identifiers must be sites A-M in order.", call. = FALSE)
  }
  raster_order <- suppressWarnings(as.integer(sub("\\.tif$", "", catalonia$source_fid)))
  if (anyNA(raster_order) || !identical(raster_order, sort(raster_order)) ||
      !identical(catalonia$raster_filename, basename(catalonia$source_fid)) ||
      !identical(catalonia$observation_id, catalonia$utm10)) {
    stop("Catalonia observation lineage or numeric raster ordering is invalid.", call. = FALSE)
  }
  if (!source_data_equal_numeric(
    observations$delta,
    observations$environmental_sample_variance / observations$environmental_mean
  )) {
    stop("Observation delta values must equal sample variance divided by the mean.", call. = FALSE)
  }
  if (any(observations$retained_value_count <= 0L)) {
    stop("Every observation must have a positive retained-value count.", call. = FALSE)
  }

  metadata <- source_data_analysis_metadata()
  families <- names(source_data_candidate_order())
  curves <- tables$fitted_mean_curves
  source_data_assert_unique_key(curves, c("analysis_id", "model_family", "prediction_index"), "fitted_mean_curves")
  curve_counts <- curves %>% dplyr::count(analysis_id, model_family)
  if (nrow(curve_counts) != 15L || any(curve_counts$n != 200L)) {
    stop("Every analysis-family curve must contain 200 predictions.", call. = FALSE)
  }
  if (!identical(unique(curves$analysis_id), metadata$analysis_id) ||
      !identical(unique(curves$model_family), families)) {
    stop("Curve analysis or candidate ordering is invalid.", call. = FALSE)
  }
  interval_rows <- curves[curves$confidence_interval_included, , drop = FALSE]
  noninterval_rows <- curves[!curves$confidence_interval_included, , drop = FALSE]
  primary_ids <- unname(figure1_primary_analysis_ids())
  if (nrow(interval_rows) != 400L ||
      any(!interval_rows$analysis_id %in% primary_ids) ||
      any(interval_rows$model_family != "logarithmic") ||
      any(interval_rows$confidence_interval_type != "pointwise_fitted_mean") ||
      any(interval_rows$confidence_level != 0.95) ||
      any(!is.finite(interval_rows$confidence_lower)) ||
      any(!is.finite(interval_rows$confidence_upper)) ||
      any(interval_rows$confidence_lower > interval_rows$confidence_upper)) {
    stop("Exactly 400 primary logarithmic rows must contain pointwise 95% fitted-mean intervals.", call. = FALSE)
  }
  interval_fields <- c("confidence_interval_type", "confidence_level", "confidence_lower", "confidence_upper")
  if (any(!vapply(noninterval_rows[interval_fields], function(x) all(is.na(x)), logical(1)))) {
    stop("Noninterval curve rows must use NA for every interval field.", call. = FALSE)
  }

  for (i in seq_len(nrow(metadata))) {
    meta <- metadata[i, , drop = FALSE]
    set <- candidate_sets[[meta$analysis_id]]
    rows <- curves[curves$analysis_id == meta$analysis_id, , drop = FALSE]
    fixed <- source_data_fixed_values(set, meta$adjustment_scenario)
    expected_fixed <- list(
      fixed_mean_foliage_height_ft = fixed$foliage,
      fixed_mean_elevation_m = fixed$elevation,
      fixed_mean_elevation_squared_m2 = fixed$elevation_squared
    )
    for (column in names(expected_fixed)) {
      expected <- expected_fixed[[column]]
      if (is.na(expected)) {
        if (!all(is.na(rows[[column]]))) stop("Unexpected fixed prediction value in ", column, ".", call. = FALSE)
      } else {
        source_data_assert_reference_compatibility(
          rows[[column]],
          rep(expected, nrow(rows)),
          round_trip_diagnostics,
          phase,
          "fitted_mean_curves",
          column
        )
      }
    }
  }

  comparisons <- tables$model_comparisons
  source_data_assert_unique_key(comparisons, c("analysis_id", "model_family"), "model_comparisons")
  if (!identical(unique(comparisons$analysis_id), metadata$analysis_id) ||
      !identical(unique(comparisons$model_family), families)) {
    stop("Model-comparison analysis or candidate ordering is invalid.", call. = FALSE)
  }
  for (analysis_id in metadata$analysis_id) {
    rows <- comparisons[comparisons$analysis_id == analysis_id, , drop = FALSE]
    expected_rank <- as.integer(rank(rows$information_criterion_value, ties.method = "min"))
    expected_best <- rows$information_criterion_value == min(rows$information_criterion_value)
    if (!identical(rows$rank, expected_rank) || !identical(rows$best_by_ic, expected_best)) {
      stop("Model-comparison ranks or exact-minimum flags are invalid for ", analysis_id, ".", call. = FALSE)
    }
    if (!identical(figure1_select_plotted_family(
      tibble::tibble(model_family = rows$model_family, delta_IC = rows$delta_ic)
    ), as.character(rows$model_family[which(rows$best_by_ic)[1]]))) {
      stop("Model-comparison tie ordering is inconsistent.", call. = FALSE)
    }
  }

  terms <- tables$model_terms
  source_data_assert_unique_key(terms, c("analysis_id", "model_family", "term_order"), "model_terms")
  expected_term_counts <- c(7L, 10L, 7L, 10L, 13L)
  observed_term_counts <- terms %>%
    dplyr::count(analysis_id) %>%
    dplyr::right_join(
      tibble::tibble(analysis_id = metadata$analysis_id),
      by = "analysis_id"
    ) %>%
    dplyr::arrange(match(analysis_id, metadata$analysis_id)) %>%
    dplyr::pull(n)
  if (!identical(observed_term_counts, expected_term_counts)) {
    stop("Model-term counts must be 7, 10, 7, 10, and 13 by analysis.", call. = FALSE)
  }
  term_groups <- split(terms, interaction(terms$analysis_id, terms$model_family, drop = TRUE))
  if (any(!vapply(term_groups, function(x) identical(x$term_order, seq_len(nrow(x))), logical(1)))) {
    stop("Model terms must preserve model-matrix order within every candidate.", call. = FALSE)
  }

  illustrative <- tables$figure1_illustrative_curves
  source_data_assert_unique_key(illustrative, c("illustrative_family", "point_index"), "figure1_illustrative_curves")
  expected_illustrative <- figure1_illustrative_curve_data()
  if (!identical(
    unique(illustrative$illustrative_family),
    c("linear", "logarithmic", "quadratic")
  )) {
    stop(phase, ": Figure 1 illustrative family order is invalid.", call. = FALSE)
  }
  illustrative_indices <- split(
    illustrative$point_index,
    factor(
      illustrative$illustrative_family,
      levels = c("linear", "logarithmic", "quadratic")
    )
  )
  if (any(!vapply(
    illustrative_indices,
    function(index) identical(index, 1:201),
    logical(1)
  ))) {
    stop(phase, ": Figure 1 illustrative point indices must be 1-201 within each family.", call. = FALSE)
  }
  if (any(!is.finite(illustrative$plotting_x)) || any(!is.finite(illustrative$plotting_y))) {
    stop(phase, ": Figure 1 illustrative coordinates must be finite.", call. = FALSE)
  }
  if (!identical(
    illustrative[c("illustrative_family", "point_index")],
    expected_illustrative[c("illustrative_family", "point_index")]
  )) {
    stop(phase, ": Figure 1 illustrative identifiers do not match the shared helper.", call. = FALSE)
  }
  for (column in c("plotting_x", "plotting_y")) {
    source_data_assert_reference_compatibility(
      illustrative[[column]],
      expected_illustrative[[column]],
      round_trip_diagnostics,
      phase,
      "figure1_illustrative_curves",
      column
    )
  }

  invisible(TRUE)
}

validate_source_data_phase <- function(
    tables, candidate_sets, phase, round_trip_diagnostics = NULL) {
  tryCatch(
    validate_source_data_tables(
      tables,
      candidate_sets,
      phase = phase,
      round_trip_diagnostics = round_trip_diagnostics
    ),
    error = function(error) {
      message_text <- conditionMessage(error)
      if (startsWith(message_text, paste0(phase, ":"))) stop(error)
      stop(phase, ": ", message_text, call. = FALSE)
    }
  )
}

normalize_figure1_readme_field <- function(x) {
  x <- gsub("`", "", as.character(x), fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}

validate_figure1_readme_mapping <- function(readme_path) {
  lines <- readLines(readme_path, warn = FALSE, encoding = "UTF-8")
  header <- "| Display order | Family | Formula | Hypothesis | Prediction | Application |"
  header_index <- which(trimws(lines) == header)
  if (length(header_index) != 1L || length(lines) < header_index + 4L) {
    stop("The source-data README must contain one complete Figure 1a mapping table.", call. = FALSE)
  }
  row_lines <- lines[header_index + 2:4]
  parsed <- lapply(row_lines, function(line) {
    inner <- sub("^\\|", "", sub("\\|$", "", trimws(line)))
    cells <- strsplit(inner, "|", fixed = TRUE)[[1]]
    if (length(cells) != 6L) {
      stop("Each Figure 1a README mapping row must contain exactly six cells.", call. = FALSE)
    }
    normalize_figure1_readme_field(cells)
  })
  parsed <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
  names(parsed) <- c(
    "display_order",
    "model_family",
    "formula",
    "hypothesis",
    "prediction",
    "application"
  )

  expected <- figure1_hypothesis_rows() %>%
    dplyr::transmute(
      display_order = as.character(display_order),
      model_family = as.character(model_family),
      formula = as.character(formula),
      hypothesis = as.character(hypothesis),
      prediction = as.character(prediction),
      application = as.character(application)
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), normalize_figure1_readme_field))

  if (!identical(parsed, as.data.frame(expected, stringsAsFactors = FALSE))) {
    stop("The Figure 1a README mapping does not match figure1_hypothesis_rows().", call. = FALSE)
  }
  invisible(TRUE)
}

source_data_file_has_lf_contract <- function(path) {
  size <- file.info(path)$size
  if (!is.finite(size) || size <= 0L) return(FALSE)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  bytes <- readBin(connection, what = "raw", n = size)
  identical(tail(bytes, 1L), as.raw(10L)) && !any(bytes == as.raw(13L))
}

write_and_reread_source_data <- function(tables, stage_dir, candidate_sets) {
  dir.create(stage_dir, recursive = FALSE, showWarnings = FALSE)
  filenames <- source_data_filenames()
  schemas <- source_data_column_types()
  serialization_tables <- source_data_serialization_tables(tables)
  staged_paths <- file.path(stage_dir, filenames)
  names(staged_paths) <- names(filenames)

  for (name in names(tables)) {
    readr::write_csv(
      serialization_tables[[name]],
      staged_paths[[name]],
      na = "NA",
      quote = "needed",
      eol = "\n"
    )
    if (!source_data_file_has_lf_contract(staged_paths[[name]])) {
      stop(name, " is not UTF-8/LF serialized with a final newline.", call. = FALSE)
    }
    serialized <- paste(readLines(staged_paths[[name]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (grepl('"NA"', serialized, fixed = TRUE)) {
      stop(name, " quotes an NA token instead of using literal unquoted NA.", call. = FALSE)
    }
    if (any(grepl(
      '(^|,)"-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:e[+-][0-9]+)?"(,|$)',
      readLines(staged_paths[[name]], warn = FALSE, encoding = "UTF-8"),
      perl = TRUE
    ))) {
      stop(name, " quotes a deterministic double token.", call. = FALSE)
    }
  }

  reread <- lapply(names(staged_paths), function(name) {
    readr::read_csv(
      staged_paths[[name]],
      col_types = schemas[[name]],
      na = "NA",
      trim_ws = FALSE,
      show_col_types = FALSE,
      progress = FALSE
    )
  })
  names(reread) <- names(staged_paths)
  round_trip_diagnostics <- source_data_compare_table_sets(
    tables,
    reread,
    phase = "staged typed CSV reread",
    serialized_tables = serialization_tables
  )
  validate_source_data_phase(
    reread,
    candidate_sets,
    phase = "staged typed CSV reread",
    round_trip_diagnostics = round_trip_diagnostics
  )
  source_data_report_round_trip_diagnostics(round_trip_diagnostics)
  list(
    paths = staged_paths,
    tables = reread,
    serialization_tables = serialization_tables,
    round_trip_diagnostics = round_trip_diagnostics
  )
}

replace_source_data_transaction <- function(staged_paths, target_paths, run_dir, validate_targets) {
  if (!identical(names(staged_paths), names(target_paths))) {
    stop("Staged and target source-data paths are not aligned.", call. = FALSE)
  }
  prior_present <- file.exists(target_paths)
  if (any(prior_present) && !all(prior_present)) {
    stop("The prior source-data set is partial; refusing to replace it.", call. = FALSE)
  }
  if (!all(file.exists(staged_paths))) {
    stop("The staged source-data set is incomplete.", call. = FALSE)
  }

  backup_dir <- file.path(run_dir, "backup")
  dir.create(backup_dir, recursive = FALSE, showWarnings = FALSE)
  backup_paths <- file.path(backup_dir, basename(target_paths))
  names(backup_paths) <- names(target_paths)
  backup_moved <- setNames(rep(FALSE, length(target_paths)), names(target_paths))
  new_installed <- setNames(rep(FALSE, length(target_paths)), names(target_paths))
  prior_md5 <- if (all(prior_present)) unname(tools::md5sum(target_paths)) else character(0)

  failure <- NULL
  validation_result <- NULL
  tryCatch({
    if (all(prior_present)) {
      for (name in names(target_paths)) {
        if (!file.rename(target_paths[[name]], backup_paths[[name]])) {
          stop("Could not back up prior source-data file: ", target_paths[[name]])
        }
        backup_moved[[name]] <- TRUE
      }
    }
    for (name in names(target_paths)) {
      if (!file.rename(staged_paths[[name]], target_paths[[name]])) {
        stop("Could not install staged source-data file: ", target_paths[[name]])
      }
      new_installed[[name]] <- TRUE
    }
    validation_result <- validate_targets(target_paths)
  }, error = function(error) {
    failure <<- error
  })

  if (!is.null(failure)) {
    removal_failures <- character(0)
    for (name in names(target_paths)[new_installed]) {
      path <- target_paths[[name]]
      if (file.exists(path) && (!file.remove(path) || file.exists(path))) {
        removal_failures <- c(removal_failures, path)
      }
    }
    restore_failures <- character(0)
    for (name in names(target_paths)[backup_moved]) {
      target <- target_paths[[name]]
      if (file.exists(target) || !file.rename(backup_paths[[name]], target)) {
        restore_failures <- c(restore_failures, target)
      }
    }
    restoration_complete <- if (all(prior_present)) {
      all(file.exists(target_paths)) &&
        identical(unname(tools::md5sum(target_paths)), prior_md5)
    } else {
      !any(file.exists(target_paths))
    }
    if (length(removal_failures) > 0L || length(restore_failures) > 0L || !restoration_complete) {
      incomplete_paths <- unique(c(removal_failures, restore_failures))
      if (!restoration_complete && length(incomplete_paths) == 0L) {
        incomplete_paths <- target_paths
      }
      condition <- structure(
        list(
          message = paste0(
            "Source-data replacement failed and rollback was incomplete for: ",
            paste(incomplete_paths, collapse = ", "),
            ". Original error: ", conditionMessage(failure)
          ),
          call = NULL
        ),
        class = c("source_data_rollback_incomplete", "error", "condition")
      )
      stop(condition)
    }
    stop(
      "Source-data replacement failed; the complete prior set was restored. Original error: ",
      conditionMessage(failure),
      call. = FALSE
    )
  }

  validation_result
}

source_data_assert_exporter_run_dir <- function(run_dir, scratch_root, repository_root) {
  normalized_repository_root <- normalizePath(
    repository_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_scratch_root <- file.path(normalized_repository_root, "_local_scratch")
  normalized_scratch_root <- normalizePath(
    scratch_root,
    winslash = "/",
    mustWork = TRUE
  )
  normalized_run_dir <- normalizePath(
    run_dir,
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(normalized_scratch_root, expected_scratch_root) ||
      !identical(dirname(normalized_run_dir), normalized_scratch_root) ||
      !startsWith(basename(normalized_run_dir), "figure_table_source_data_export_")) {
    stop("Refusing to remove a path outside the exporter-owned scratch transaction.", call. = FALSE)
  }
  invisible(normalized_run_dir)
}

source_data_remove_exporter_run_dir <- function(run_dir, scratch_root, repository_root) {
  source_data_assert_exporter_run_dir(run_dir, scratch_root, repository_root)
  unlink(run_dir, recursive = TRUE, force = FALSE)
  if (dir.exists(run_dir) || file.exists(run_dir)) {
    stop("Could not remove the exporter-owned transaction directory: ", run_dir, call. = FALSE)
  }
  invisible(TRUE)
}

source_data_remove_new_empty_scratch <- function(scratch_root, repository_root) {
  normalized_repository_root <- normalizePath(
    repository_root,
    winslash = "/",
    mustWork = TRUE
  )
  expected_scratch_root <- file.path(normalized_repository_root, "_local_scratch")
  normalized_scratch_root <- normalizePath(
    scratch_root,
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(normalized_scratch_root, expected_scratch_root)) {
    stop("Refusing to remove a scratch directory outside this exporter checkout.", call. = FALSE)
  }
  remaining <- list.files(
    normalized_scratch_root,
    all.files = TRUE,
    no.. = TRUE,
    full.names = TRUE
  )
  if (length(remaining) != 0L) {
    stop(
      "The exporter-created _local_scratch directory is not empty after cleanup: ",
      paste(remaining, collapse = ", "),
      call. = FALSE
    )
  }
  unlink(normalized_scratch_root, recursive = TRUE, force = FALSE)
  if (dir.exists(normalized_scratch_root) || file.exists(normalized_scratch_root)) {
    stop("Could not remove the verified empty exporter scratch directory.", call. = FALSE)
  }
  invisible(TRUE)
}

run_figure_table_source_data_exporter <- function() {
  repository_root <- find_source_data_project_root()
  switches <- list(
    RUN_WORKFLOW_ON_SOURCE = FALSE,
    INSTALL_MISSING = FALSE,
    RUN_EXACT_REPLICATION = FALSE,
    RUN_MANUSCRIPT_HIERARCHY = FALSE,
    RUN_MANUSCRIPT_DELIVERABLES = FALSE,
    WRITE_SETUP_METADATA = FALSE,
    WRITE_SOURCE_PREP_OUTPUTS = FALSE,
    RUN_PV_AUDIT_MODE = FALSE,
    WRITE_PV_AUDIT_INTERMEDIATES = FALSE,
    RUN_EXPENSIVE_GLOBAL_RESAMPLING = FALSE,
    PV_AUDIT_OUTPUT_ROOT = NULL,
    PV_AUDIT_RUN_ID = NULL
  )
  list2env(switches, envir = globalenv())
  source(file.path(repository_root, "code", "run_delta_hdr_brief_reanalysis.R"), local = globalenv())
  if (!identical(
    normalizePath(project_root, winslash = "/", mustWork = TRUE),
    repository_root
  )) {
    stop(
      "The workflow project root does not match the exporter checkout: ",
      project_root,
      call. = FALSE
    )
  }

  scratch_existed_before <- dir.exists(scratch_root)
  acquire_workflow_lock()
  lock_acquired <- TRUE
  run_dir <- NULL
  preserve_run_dir <- FALSE
  on.exit({
    if (isTRUE(lock_acquired)) {
      release_workflow_lock()
      if (file.exists(workflow_lock_path)) {
        stop("Could not release the source-data exporter workflow lock.", call. = FALSE)
      }
      lock_acquired <- FALSE
    }
    if (!isTRUE(preserve_run_dir) && !is.null(run_dir) && dir.exists(run_dir)) {
      source_data_remove_exporter_run_dir(run_dir, scratch_root, repository_root)
    }
    if (!scratch_existed_before && !isTRUE(preserve_run_dir) && dir.exists(scratch_root)) {
      source_data_remove_new_empty_scratch(scratch_root, repository_root)
    }
  }, add = TRUE)

  prior_exporter_residue <- list.files(
    scratch_root,
    pattern = "^figure_table_source_data_export_",
    full.names = TRUE,
    all.files = TRUE
  )
  if (length(prior_exporter_residue) > 0L) {
    stop(
      "Exporter transaction evidence already exists beneath _local_scratch; inspect it before rerunning.",
      call. = FALSE
    )
  }
  run_dir <- tempfile("figure_table_source_data_export_", tmpdir = scratch_root)
  if (!dir.create(run_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create the exporter-owned run directory.", call. = FALSE)
  }

  prepared <- run_manuscript_input_preparation(
    root = repo_root,
    write_source_prep_outputs = FALSE
  )
  assign("results_full_R", prepared, envir = globalenv())
  hierarchy <- run_manuscript_analysis_hierarchy(write_outputs = FALSE)
  candidate_sets <- source_data_candidate_sets(hierarchy)

  if (nrow(hierarchy$manifest) != 9L || !all(is.na(hierarchy$manifest$exists))) {
    stop("The no-write hierarchy manifest must retain nine paths and use exists = NA.", call. = FALSE)
  }
  foliage_samples <- load_foliage_samples(repo_root)

  tables <- list(
    observations = build_source_data_observations(prepared, foliage_samples),
    fitted_mean_curves = build_source_data_fitted_mean_curves(candidate_sets),
    model_comparisons = build_source_data_model_comparisons(candidate_sets),
    model_terms = build_source_data_model_terms(candidate_sets),
    figure1_illustrative_curves = figure1_illustrative_curve_data()
  )
  validate_source_data_phase(
    tables,
    candidate_sets,
    phase = "in-memory construction"
  )

  target_dir <- file.path(project_root, "results", "delta_hdr_reanalysis", "source_data")
  readme_path <- file.path(target_dir, "README.md")
  if (!file.exists(readme_path)) {
    stop("The source-data README is missing: ", readme_path, call. = FALSE)
  }
  validate_figure1_readme_mapping(readme_path)

  staged <- write_and_reread_source_data(
    tables,
    stage_dir = file.path(run_dir, "stage"),
    candidate_sets = candidate_sets
  )
  target_paths <- file.path(target_dir, source_data_filenames())
  names(target_paths) <- names(source_data_filenames())
  final_validation <- tryCatch(
    replace_source_data_transaction(
      staged$paths,
      target_paths,
      run_dir,
      validate_targets = function(paths) {
        validated <- lapply(names(paths), function(name) {
          readr::read_csv(
            paths[[name]],
            col_types = source_data_column_types()[[name]],
            na = "NA",
            trim_ws = FALSE,
            show_col_types = FALSE,
            progress = FALSE
          )
        })
        names(validated) <- names(paths)
        installed_diagnostics <- source_data_compare_table_sets(
          tables,
          validated,
          phase = "installed target typed CSV reread",
          serialized_tables = staged$serialization_tables
        )
        source_data_report_round_trip_diagnostics(installed_diagnostics)
        validate_source_data_phase(
          validated,
          candidate_sets,
          phase = "installed target typed CSV reread",
          round_trip_diagnostics = installed_diagnostics
        )
        for (name in names(paths)) {
          if (!source_data_file_has_lf_contract(paths[[name]])) {
            stop("Final source-data serialization failed for ", name, ".", call. = FALSE)
          }
        }
        list(tables = validated, diagnostics = installed_diagnostics)
      }
    ),
    source_data_rollback_incomplete = function(error) {
      preserve_run_dir <<- TRUE
      stop(error)
    }
  )
  final_tables <- final_validation$tables

  message("Figure and table source-data export complete. Files: ", target_dir)
  invisible(list(
    files = target_paths,
    row_counts = vapply(final_tables, nrow, integer(1)),
    staged_round_trip_diagnostics = staged$round_trip_diagnostics,
    installed_round_trip_diagnostics = final_validation$diagnostics,
    scratch_existed_before = scratch_existed_before
  ))
}

if (sys.nframe() == 0L) {
  run_figure_table_source_data_exporter()
}
