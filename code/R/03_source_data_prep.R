# Prepares inputs shared by the main analysis and optional reimplementation.

calculate_elevation_stats <- function(filename, path) {
  raster_path <- safe_file(path, filename)
  if (!file.exists(raster_path)) stop("Required elevation raster is missing: ", raster_path)
  r <- terra::rast(raster_path)
  v <- as.numeric(terra::values(r, mat = FALSE))
  v <- v[is.finite(v) & v > 0]
  n <- length(v)
  if (n == 0) {
    return(tibble(fid = filename, mu = 0, range = 0, sigma2 = 0, CV = 0, delta = 0, n = 0))
  }
  mu <- mean(v); sigma2 <- var(v)
  tibble(fid = filename, mu = mu, range = max(v) - min(v), sigma2 = sigma2,
         CV = sqrt(sigma2) / mu, delta = sigma2 / mu, n = n)
}

bird_richness <- function(birds) {
  birds %>%
    pivot_longer(cols = -UTM10, names_to = "species", values_to = "value", values_drop_na = TRUE) %>%
    mutate(value = suppressWarnings(as.numeric(value))) %>%
    filter(!is.na(value)) %>%
    group_by(UTM10) %>%
    summarise(richness = sum(value > 0), .groups = "drop")
}

pv_audit_csv_schema_row <- function(input_name, path, data) {
  tibble(
    input_name = input_name,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    row_count = nrow(data),
    column_count = ncol(data),
    column_names_or_raster_count = paste(names(data), collapse = "|"),
    checksum_or_bundle_hash = NA_character_
  )
}

pv_audit_allouche_raw_schema <- function(coordinates_path, coordinates, birds_path, birds, elevgrid_path, filenames) {
  bind_rows(
    pv_audit_csv_schema_row("coordinates_csv", coordinates_path, coordinates),
    pv_audit_csv_schema_row("birds_csv", birds_path, birds),
    tibble(
      input_name = "elevation_raster_bundle",
      path = normalizePath(elevgrid_path, winslash = "/", mustWork = FALSE),
      row_count = NA_integer_,
      column_count = NA_integer_,
      column_names_or_raster_count = paste0("raster_count=", length(filenames)),
      checksum_or_bundle_hash = NA_character_
    )
  )
}

pv_audit_macarthur_raw_manifest <- function(root = repo_root) {
  purrr::map_dfr(LETTERS[1:13], function(site) {
    path <- safe_file(root, "05_review", "data", "digitised_points", paste0(site, ".csv"))
    raw <- readr::read_csv(path, show_col_types = FALSE)
    tibble(
      site = site,
      path = normalizePath(path, winslash = "/", mustWork = FALSE),
      row_count = nrow(raw),
      column_names = paste(names(raw), collapse = "|"),
      checksum = NA_character_
    )
  })
}

pv_audit_macarthur_loess_layers <- function(density_df) {
  purrr::map_dfr(split(density_df, density_df$site), function(df) {
    site <- unique(df$site)[1]
    observed <- df %>%
      transmute(
        site = site,
        panel = "observed_points",
        x_density = x,
        y_height = y,
        predicted_density = NA_real_,
        loess_span = NA_real_,
        degree = NA_real_
      )
    us <- seq(min(df$y, na.rm = TRUE), max(df$y, na.rm = TRUE), by = 0.1)
    model <- stats::loess(x ~ y, data = df, span = 0.1, degree = 1)
    predicted <- tibble(
      site = site,
      panel = "loess_prediction",
      x_density = NA_real_,
      y_height = us,
      predicted_density = as.numeric(stats::predict(model, newdata = data.frame(y = us))),
      loess_span = 0.1,
      degree = 1
    )
    bind_rows(observed, predicted)
  })
}

macarthur_density_panel_specs <- function() {
  tibble(
    site = LETTERS[1:13],
    site_label = c(
      "A. Vermont", "B. Penn.", "C. Penn", "D. Panama", "E. Maine", "F. Vermont",
      "G. Vermont", "H. Vermont", "I. Florida", "J. Maryland", "K. Penn.",
      "L. Penn.", "M. Penn."
    ),
    y_min = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    y_max = c(20, 20, 40, 20, 68, 62, 61, 60, 40, 80, 100, 100, 100),
    split_x = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),
    x_left_min = c(0, 0, 0, 0, NA, NA, NA, 0, 0, 0, NA, NA, 0),
    x_left_max = c(0.024, 0.025, 0.063, 0.063, NA, NA, NA, 0.026, 0.085, 0.032, NA, NA, 0.15),
    x_right_min = c(0.04, 0.04, 0.07, 0.07, NA, NA, NA, 0.035, 0.095, 0.04, NA, NA, 0.2),
    x_right_max = c(0.48, 0.734, 0.670, 1.38, NA, NA, NA, 0.3, 0.178, 0.164, NA, NA, 0.4),
    x_left_breaks = list(
      c(0, 0.02), c(0, 0.02), c(0, 0.02, 0.04, 0.06), c(0, 0.02, 0.04, 0.06),
      NULL, NULL, NULL, c(0, 0.02), seq(0, 0.08, by = 0.02), c(0, 0.02),
      NULL, NULL, c(0, 0.1)
    ),
    x_right_breaks = list(
      0.48, 0.734, 0.670, 1.38, NULL, NULL, NULL, 0.3, 0.178, 0.164,
      NULL, NULL, 0.4
    ),
    right_pane_width = c(0.10, 0.10, 0.12, 0.12, NA, NA, NA, 0.12, 0.13, 0.12, NA, NA, 0.12),
    right_axis_text_size = c(4.8, 4.8, 5.0, 5.0, NA, NA, NA, 5.0, 5.0, 5.0, NA, NA, 5.0),
    y_breaks = list(
      c(0, 20), c(0, 20), c(0, 20, 40), c(0, 20), NULL, NULL, NULL,
      c(0, 20, 40, 60), c(0, 20, 40), seq(0, 80, by = 20), NULL, NULL,
      seq(0, 100, by = 20)
    ),
    layout_t = c(1, 2, 3, 5, 6, 6, 6, 6, 1, 3, 3, 1, 1),
    layout_l = c(1, 1, 1, 1, 1, 4, 9, 11, 3, 5, 7, 9, 11),
    layout_b = c(1, 2, 4, 5, 8, 8, 8, 8, 2, 5, 8, 5, 5),
    layout_r = c(2, 2, 4, 4, 3, 6, 10, 12, 7, 6, 8, 10, 12)
  )
}

macarthur_density_prediction_data <- function(loess_layers) {
  loess_layers %>%
    filter(panel == "loess_prediction") %>%
    transmute(
      site,
      y_height = y_height,
      predicted_density = predicted_density
    ) %>%
    filter(is.finite(y_height), is.finite(predicted_density)) %>%
    arrange(site, y_height)
}

macarthur_density_band_polygon_data <- function(plot_data, group_id) {
  plot_data <- plot_data %>%
    filter(is.finite(y_height), is.finite(predicted_density)) %>%
    arrange(y_height)
  if (nrow(plot_data) == 0) {
    return(tibble(x = numeric(), y = numeric(), group = character()))
  }
  tibble(
    x = c(rep(0, nrow(plot_data)), rev(plot_data$predicted_density)),
    y = c(plot_data$y_height, rev(plot_data$y_height)),
    group = group_id
  )
}

macarthur_density_axis_break_marks <- function(p, x_limits, y_limits, side = c("left", "right")) {
  side <- match.arg(side)
  if (is.null(x_limits) || any(!is.finite(x_limits)) || any(!is.finite(y_limits))) {
    return(p)
  }
  dx <- diff(x_limits) * 0.020
  dy <- diff(y_limits) * 0.018
  if (!is.finite(dx) || dx <= 0 || !is.finite(dy) || dy <= 0) return(p)
  if (identical(side, "right")) {
    x0 <- x_limits[2] - dx
    x1 <- x_limits[2]
  } else {
    x0 <- x_limits[1]
    x1 <- x_limits[1] + dx
  }
  p +
    annotate("segment", x = x0, xend = x1, y = y_limits[1] + dy, yend = y_limits[1],
             linewidth = 0.11, colour = "grey50") +
    annotate("segment", x = x0, xend = x1, y = y_limits[2], yend = y_limits[2] - dy,
             linewidth = 0.11, colour = "grey50")
}

macarthur_density_profile_plot <- function(
  plot_data,
  spec,
  x_limits = NULL,
  x_breaks = NULL,
  show_y_axis = TRUE,
  show_title = TRUE,
  break_side = NULL,
  pane_id = "single",
  x_text_size = 6,
  panel_border_colour = "grey55",
  panel_border_linewidth = 0.28
) {
  y_limits <- c(spec$y_min[[1]], spec$y_max[[1]])
  band_data <- macarthur_density_band_polygon_data(plot_data, paste0(spec$site[[1]], "_", pane_id))
  line_data <- plot_data %>%
    filter(is.finite(y_height), is.finite(predicted_density)) %>%
    arrange(y_height)
  pane_margin <- switch(
    pane_id,
    left = margin(0.6, 0.05, 0.6, 0.6),
    right = margin(0.6, 0.6, 0.6, 0.05),
    margin(0.6, 0.6, 0.6, 0.6)
  )
  base_theme <- theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.93, margin = margin(b = 0.5)),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = panel_border_linewidth),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = x_text_size, margin = margin(t = 0.5)),
      axis.text.y = element_text(size = 7),
      axis.ticks.length = grid::unit(0.9, "pt"),
      plot.margin = pane_margin
    )
  if (!isTRUE(show_y_axis)) {
    base_theme <- base_theme +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  p <- ggplot() +
    geom_polygon(
      data = band_data,
      aes(x = x, y = y, group = group),
      inherit.aes = FALSE,
      fill = "grey75",
      alpha = 0.9,
      colour = NA
    ) +
    geom_path(
      data = line_data,
      aes(x = predicted_density, y = y_height),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "#0072B2",
      lineend = "round",
      na.rm = TRUE
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +
    labs(y = NULL, title = if (isTRUE(show_title)) spec$site_label[[1]] else NULL) +
    base_theme
  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks, guide = guide_axis(check.overlap = TRUE))
  }
  y_breaks <- spec$y_breaks[[1]]
  if (!is.null(y_breaks)) {
    p <- p + scale_y_continuous(breaks = y_breaks)
  }
  if (!is.null(break_side)) {
    p <- macarthur_density_axis_break_marks(p, x_limits, y_limits, break_side)
  }
  p
}

macarthur_density_profile_split_plot <- function(plot_data, spec) {
  left_limits <- c(spec$x_left_min[[1]], spec$x_left_max[[1]])
  right_limits <- c(spec$x_right_min[[1]], spec$x_right_max[[1]])
  right_pane_width <- spec$right_pane_width[[1]]
  right_axis_text_size <- spec$right_axis_text_size[[1]]
  if (!is.finite(right_pane_width)) right_pane_width <- 0.12
  if (!is.finite(right_axis_text_size)) right_axis_text_size <- 5
  left <- macarthur_density_profile_plot(
    plot_data,
    spec,
    left_limits,
    x_breaks = spec$x_left_breaks[[1]],
    show_y_axis = TRUE,
    show_title = TRUE,
    break_side = "right",
    pane_id = "left"
  )
  right <- macarthur_density_profile_plot(
    plot_data,
    spec,
    right_limits,
    x_breaks = spec$x_right_breaks[[1]],
    show_y_axis = FALSE,
    show_title = FALSE,
    break_side = "left",
    pane_id = "right",
    x_text_size = right_axis_text_size,
    panel_border_colour = "grey72",
    panel_border_linewidth = 0.18
  ) +
    labs(title = NULL)
  left + right + patchwork::plot_layout(widths = c(1, right_pane_width))
}

macarthur_density_site_plot <- function(plot_data, spec) {
  if (isTRUE(spec$split_x[[1]])) {
    macarthur_density_profile_split_plot(plot_data, spec)
  } else {
    macarthur_density_profile_plot(plot_data, spec, x_limits = NULL, show_y_axis = TRUE)
  }
}

macarthur_density_panel_figure <- function(loess_layers) {
  specs <- macarthur_density_panel_specs()
  density_data <- macarthur_density_prediction_data(loess_layers)
  site_plots <- purrr::map(specs$site, function(site_id) {
    spec <- filter(specs, site == site_id)
    plot_data <- filter(density_data, site == site_id)
    macarthur_density_site_plot(plot_data, spec)
  })
  names(site_plots) <- specs$site
  design <- purrr::pmap(
    specs %>% select(layout_t, layout_l, layout_b, layout_r),
    function(layout_t, layout_l, layout_b, layout_r) {
      patchwork::area(t = layout_t, l = layout_l, b = layout_b, r = layout_r)
    }
  )
  density_grid <- patchwork::wrap_plots(
    site_plots,
    design = do.call(c, design),
    widths = rep(1, 12),
    heights = c(1, 1, 1, 0.22, 1, 1, 1, 0.22)
  ) &
    theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4))
  panel_label <- patchwork::wrap_elements(
    full = grid::textGrob(
      "d",
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.98, "npc"),
      just = c("center", "top"),
      gp = grid::gpar(fontsize = 11, fontface = "bold")
    )
  )
  y_label <- patchwork::wrap_elements(
    full = grid::textGrob(
      "Foliage height (feet)",
      rot = 90,
      gp = grid::gpar(fontsize = 9)
    )
  )
  x_label <- patchwork::wrap_elements(
    full = grid::textGrob(
      "Foliage density",
      gp = grid::gpar(fontsize = 9)
    )
  )
  top_row <- y_label + density_grid +
    patchwork::plot_layout(widths = c(0.020, 1))
  bottom_row <- patchwork::plot_spacer() + x_label +
    patchwork::plot_layout(widths = c(0.020, 1))
  density_body <- top_row / bottom_row +
    patchwork::plot_layout(heights = c(1, 0.040))
  panel <- panel_label + density_body +
    patchwork::plot_layout(widths = c(0.030, 1))
  panel +
    patchwork::plot_annotation(
      theme = theme(
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 0.35),
        plot.margin = margin(2, 2, 2, 2)
      )
    )
}

macarthur_full_figure <- function(scatter_component, density_component) {
  scatter_component / density_component +
    patchwork::plot_layout(heights = c(0.72, 1.60)) +
    patchwork::plot_annotation(
      theme = theme(
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        plot.margin = margin(2, 2, 2, 2)
      )
    )
}

macarthur_scatter_label_positions <- function(metrics, x_var) {
  x_var_name <- x_var
  x_values <- metrics[[x_var]]
  x_span <- diff(range(x_values, na.rm = TRUE))
  y_span <- diff(range(metrics$bird_species_diversity, na.rm = TRUE))
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  label_offsets <- tibble(
    site = LETTERS[1:13],
    nudge_x_frac = c(-0.018, 0.018, -0.020, 0.018, 0.016, 0.014, 0.018, 0.014, 0.014, 0.014, 0.016, -0.016, 0.014),
    nudge_y_frac = c(-0.035, -0.030, 0.032, 0.030, -0.032, 0.030, -0.034, 0.034, -0.032, 0.030, 0.032, -0.030, 0.032)
  )
  variable_offsets <- tibble(
    x_var = "diff_entropy",
    site = c("E", "F", "G", "J", "K", "M"),
    nudge_x_frac_var = c(-0.030, 0.026, 0.030, 0.038, -0.032, 0.024),
    nudge_y_frac_var = c(-0.040, 0.040, -0.040, 0.014, 0.004, 0.042)
  ) %>%
    filter(.data$x_var == x_var_name) %>%
    select(-x_var)
  metrics %>%
    left_join(label_offsets, by = "site") %>%
    left_join(variable_offsets, by = "site") %>%
    mutate(
      nudge_x_frac = coalesce(nudge_x_frac_var, nudge_x_frac, 0.015),
      nudge_y_frac = coalesce(nudge_y_frac_var, nudge_y_frac, 0.030),
      label_x = .data[[x_var]] + coalesce(nudge_x_frac, 0.015) * x_span,
      label_y = bird_species_diversity + coalesce(nudge_y_frac, 0.030) * y_span
    )
}

macarthur_scatter_plot <- function(metrics, x_var, x_label, panel_tag = NULL) {
  label_data <- macarthur_scatter_label_positions(metrics, x_var)
  ggplot(metrics, aes(x = .data[[x_var]], y = bird_species_diversity)) +
    geom_point() +
    geom_text(
      data = label_data,
      aes(x = label_x, y = label_y, label = site),
      size = 3,
      inherit.aes = FALSE
    ) +
    labs(x = x_label, y = "Bird species diversity", tag = panel_tag) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(3, 7, 3, 4)
    )
}

prepare_allouche_data_R <- function(
  root = repo_root,
  min_n = 14000,
  write_outputs = TRUE,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  data_path <- safe_file(root, "03_corrected_HDR", "empirical", "01_carnicer", "data")
  elevgrid_path <- safe_file(data_path, "elevation", "elevgrid")
  bird_path <- safe_file(data_path, "bird")
  filenames <- paste0(0:385, ".tif")
  elevdata <- purrr::map_dfr(filenames, calculate_elevation_stats, path = elevgrid_path)
  coordinates_path <- safe_file(bird_path, "coordinates.csv")
  birds_path <- safe_file(bird_path, "birds.csv")
  coordinates <- readr::read_csv(coordinates_path, show_col_types = FALSE)
  birds <- readr::read_csv(birds_path, na = c("", "NA", "NaN"), show_col_types = FALSE)
  # Pair raster summaries positionally with coordinate rows in numeric filename
  # order (0.tif through 385.tif), then join bird richness by UTM10.
  adf <- bird_richness(birds) %>%
    left_join(bind_cols(coordinates, elevdata), by = c("UTM10" = "Utm10")) %>%
    drop_na() %>%
    mutate(rich_per_area = richness / n) %>%
    filter(n > min_n)
  if (isTRUE(write_outputs)) {
    readr::write_csv(adf, safe_file(intermediate_dir, "allouche_prepared_data_R.csv"))
  }
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    pv_audit_write_contract_csv(
      pv_audit_allouche_raw_schema(coordinates_path, coordinates, birds_path, birds, elevgrid_path, filenames),
      "EC-011",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Catalonia raw input schema summary from R prep"
    )
  }
  pv_audit_write_contract_csv(
    adf %>%
      transmute(
        UTM10,
        richness,
        Utm10 = UTM10,
        X_coord,
        Y_coord,
        fid,
        mu,
        range,
        sigma_or_sigma2 = sigma2,
        CV,
        delta,
        n,
        rich_per_area
      ),
    "EC-012",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Prepared Catalonia adf dataframe"
  )
  adf
}

run_allouche_analysis_R <- function(
  adf,
  write_outputs = TRUE,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  allouche_model_df <- adf
  if (all(c("X_coord", "Y_coord") %in% names(allouche_model_df))) {
    allouche_model_df <- allouche_model_df %>% mutate(x = rescale01(X_coord), y = rescale01(Y_coord))
  }
  models <- list(
    range_rm = lm(richness ~ mu + range + I(range^2), data = allouche_model_df),
    delta_rm = lm(richness ~ mu + delta + I(delta^2), data = allouche_model_df),
    range2_rm = lm(richness ~ mu + I(mu^2) + range + I(range^2), data = allouche_model_df),
    delta2_rm = lm(richness ~ mu + I(mu^2) + delta + I(delta^2), data = allouche_model_df),
    delta_fin = lm(richness ~ mu + I(mu^2) + delta, data = allouche_model_df)
  )
  null_model <- lm(richness ~ 1, data = allouche_model_df)
  p_values <- tibble(
    model = names(models),
    p_value_vs_null = map_dbl(models, ~nested_f_p(null_model, .x))
  ) %>% bind_rows(tibble(
    model = "delta_fin_vs_delta2_rm_nested_test",
    p_value_vs_null = as.numeric(anova(models$delta_fin, models$delta2_rm)$`Pr(>F)`[2])
  ))
  if (isTRUE(write_outputs)) {
    readr::write_csv(p_values, safe_file(tab_dir, "allouche_lm_p_values_R.csv"))
    write_model_terms(models, "allouche_lm_model_terms_R.csv")
    write_model_fit(models, "allouche_lm_model_fit_R.csv")
  }
  pv_audit_write_contract_csv(
    allouche_model_df %>%
      mutate(row_id = dplyr::row_number()) %>%
      transmute(
        row_id,
        richness,
        mu,
        delta,
        range,
        mu_rs = rescale01(mu),
        delta_rs = rescale01(delta),
        range_rs = rescale01(range),
        x = if ("x" %in% names(allouche_model_df)) x else NA_real_,
        y = if ("y" %in% names(allouche_model_df)) y else NA_real_
      ),
    "EC-013",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Allouche model input dataframe"
  )
  pv_audit_write_contract_csv(
    p_values %>%
      mutate(source_order = dplyr::row_number()) %>%
      transmute(
        model,
        test_type = ifelse(model == "delta_fin_vs_delta2_rm_nested_test", "nested_model_f_test", "vs_null_f_test"),
        p_value = p_value_vs_null,
        source_order,
        notes = NA_character_
      ),
    "EC-015",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Allouche labeled F-test p-values"
  )
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    pv_audit_write_contract_csv(
      purrr::imap_dfr(models, function(model, name) {
        broom::tidy(model) %>%
          transmute(
            model = name,
            term,
            estimate,
            std_error = std.error,
            statistic,
            p_value = p.value,
            predictor_scale = "raw_scale"
          )
      }),
      "EC-014",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Allouche model coefficient table"
    )
  }

  theory_df <- adf %>% filter(n > 14000) %>% mutate(richness = as.numeric(richness))
  theory_models <- list(
    richness_mu_mu2 = lm(richness ~ mu + I(mu^2), data = theory_df),
    richness_delta = lm(richness ~ delta, data = theory_df),
    richness_delta_mu_mu2 = lm(richness ~ delta + mu + I(mu^2), data = theory_df)
  )
  if (isTRUE(write_outputs)) {
    write_model_terms(theory_models, "allouche_theory_raw_scale_model_terms_R.csv")
  }
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    invisible(theory_models)
  }

  p_range <- ggplot(allouche_model_df, aes(range, richness)) + geom_point(alpha = 0.65) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE) +
    labs(x = "Elevation range, raw scale", y = "Bird species richness") + theme_bw()
  p_delta <- ggplot(allouche_model_df, aes(delta, richness)) + geom_point(alpha = 0.65) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    labs(x = "delta elevation, raw scale", y = "Bird species richness") + theme_bw()
  if (isTRUE(write_outputs)) {
    save_png(p_range | p_delta, "allouche_range_vs_delta_diagnostics_R.png", width = 10, height = 4.5)
  }
  list(allouche_model_df = allouche_model_df, models = models, p_values = p_values, theory_models = theory_models)
}


read_digitized_points <- function(root, letter) {
  path <- safe_file(root, "05_review", "data", "digitised_points", paste0(letter, ".csv"))
  df <- readr::read_csv(path, show_col_types = FALSE)
  df <- df[, 1:2]
  names(df) <- c("x", "y")
  df %>% mutate(x = as.numeric(x), y = as.numeric(y))
}
add_points <- function(df, x, y) bind_rows(df, tibble(x = as.numeric(x), y = as.numeric(y)))
clean_digitized_data <- function(data) {
  data <- data %>% mutate(x = if_else(x < 0, 0, x), y = if_else(y < 0, 0, y))
  measure_y <- seq(5, 60, by = 5)
  y_add <- measure_y[measure_y > max(data$y, na.rm = TRUE)]
  if (length(y_add) > 0) data <- bind_rows(tibble(x = rep(0, length(y_add)), y = y_add), data)
  data
}
load_macarthur_digitized_data <- function(
  root = repo_root,
  write_outputs = TRUE,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  a <- read_digitized_points(root, "A") %>% add_points(0.48, 0.5) %>% clean_digitized_data()
  b <- read_digitized_points(root, "B") %>% add_points(0.734, 0.5) %>% clean_digitized_data()
  c <- read_digitized_points(root, "C") %>% add_points(0.670, 0.5) %>% clean_digitized_data()
  d <- read_digitized_points(root, "D") %>% add_points(c(1.38, 1.38), c(2.0, 0.5)) %>% clean_digitized_data() %>%
    add_points(rep(1.38, 6), seq(2.0, 0.5, length.out = 6))
  e <- read_digitized_points(root, "E"); e <- e %>% add_points(e$x[nrow(e)], 0.5) %>% clean_digitized_data()
  f <- read_digitized_points(root, "F") %>% clean_digitized_data()
  g <- read_digitized_points(root, "G") %>% clean_digitized_data()
  keep_idx <- c(seq_len(min(91, nrow(g))), if (nrow(g) > 115) 116:nrow(g) else integer(0))
  g <- g[keep_idx, 1:2, drop = FALSE] %>% add_points(c(0.0053, 0.0055, 0.006, 0.0061, 0.0062, 0.029), c(2.55, 2.35, 2.05, 2.0, 1.95, 1.05))
  h <- read_digitized_points(root, "H") %>% add_points(c(0.028, 0.3), c(2.0, 0.5)) %>% clean_digitized_data() %>% add_points(c(0.3, 0.3), c(0.45, 0.4))
  i <- read_digitized_points(root, "I") %>% add_points(0.178, 0.5) %>% clean_digitized_data()
  j <- read_digitized_points(root, "J") %>% add_points(0.164, 0.5) %>% clean_digitized_data()
  k <- read_digitized_points(root, "K") %>% add_points(0.139, 0.5) %>% clean_digitized_data()
  l <- read_digitized_points(root, "L") %>% clean_digitized_data() %>% add_points(0.0389, 2.0)
  m <- read_digitized_points(root, "M") %>% add_points(0.4, 0.5) %>% clean_digitized_data() %>% add_points(c(0.4, 0.4), c(0.48, 0.45))
  out <- list(A=a,B=b,C=c,D=d,E=e,F=f,G=g,H=h,I=i,J=j,K=k,L=l,M=m)
  if (isTRUE(write_outputs)) {
    imap_dfr(out, ~ mutate(.x, site = .y, .before = 1)) %>% write_csv(safe_file(intermediate_dir, "macarthur_digitized_points_cleaned_R.csv"))
  }
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    pv_audit_write_contract_csv(
      pv_audit_macarthur_raw_manifest(root),
      "EC-024",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "MacArthur raw digitized point bundle manifest"
    )
    pv_audit_write_contract_csv(
      imap_dfr(out, function(data, site) {
        data %>%
          mutate(
            site = site,
            row_id = dplyr::row_number(),
            edit_source = "R_port_clean_digitized_data_and_manual_points",
            .before = 1
          ) %>%
          select(site, row_id, x, y, edit_source)
      }),
      "EC-025",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "MacArthur cleaned digitized points"
    )
  }
  out
}
# Read the first 13 columns of the fixed archived foliage table as sites A–M.
load_foliage_samples <- function(root = repo_root) {
  df <- readr::read_csv(safe_file(root, "05_review", "data", "foliage_density.csv"), show_col_types = FALSE)
  map(df[1:13], ~ {x <- as.numeric(.x); x[is.finite(x)]})
}
FHD3samp <- function(y) { p <- c(sum(y <= 2), sum(y > 2 & y <= 25), sum(y > 25)) / length(y); p <- p[p > 0]; -sum(p * log(p)) }
FHD51samp <- function(y) { edges <- seq(0,100,by=2); p <- map_dbl(seq_len(length(edges)-1), ~sum(y > edges[.x] & y <= edges[.x+1]) / length(y)); p <- p[p>0]; -sum(p*log(p)) }
calc_gini_sample <- function(y) sum(abs(outer(y, y, "-"))) / (2 * length(y)^2 * mean(y))
kl_entropy_1d <- function(y, k = 1, jitter_sd = 1e-10) {
  y <- as.numeric(y); y <- y[is.finite(y)]; n <- length(y)
  if (n <= k) return(NA_real_)
  if (any(duplicated(y))) { set.seed(1); y <- y + rnorm(n, sd = jitter_sd) }
  d <- abs(outer(y, y, "-")); diag(d) <- Inf
  eps <- apply(d, 1, function(row) sort(row, partial = k)[k])
  eps <- pmax(eps, .Machine$double.eps)
  digamma(n) - digamma(k) + log(2) + mean(log(eps))
}
calc_mac_metrics <- function(y) {
  y <- as.numeric(y); y <- y[is.finite(y)]
  mu <- mean(y); sigma2 <- var(y); sigma <- sqrt(sigma2)
  tibble(mu = mu, sigma2 = sigma2, delta = sigma2 / mu, FHD3 = FHD3samp(y), FH_entropy_51 = FHD51samp(y),
         gini = calc_gini_sample(y), range = max(y) - min(y), sigma = sigma, CV = sigma / mu, diff_entropy = kl_entropy_1d(y))
}
cor_row <- function(name, x, y) {
  ct <- cor.test(x, y, method = "pearson")
  tibble(test = name, estimate_r = as.numeric(ct$estimate), conf_low = ct$conf.int[1], conf_high = ct$conf.int[2],
         p_value = as.numeric(ct$p.value), statistic_t = as.numeric(ct$statistic), df = as.numeric(ct$parameter), n = length(x))
}
run_macarthur_analysis_R <- function(
  root = repo_root,
  write_outputs = TRUE,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  digitized <- load_macarthur_digitized_data(
    root,
    write_outputs = write_outputs,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_audit_intermediates = write_audit_intermediates
  )
  samples <- load_foliage_samples(root)
  # MacArthur and MacArthur (1961) calculated bird species diversity with the
  # natural logarithm. Pair their published Shannon index values with foliage
  # sites A through M in the same order.
  bsd <- c(0.639,1.266,2.265,1.906,1.712,2.403,1.721,2.739,1.332,2.285,2.277,2.127,2.567)
  metrics <- map_dfr(samples, calc_mac_metrics) %>% mutate(site = LETTERS[1:13], bird_species_diversity = bsd, .before = 1)
  correlations <- bind_rows(
    cor_row("Mean-BSD", metrics$mu, metrics$bird_species_diversity),
    cor_row("delta-BSD", metrics$delta, metrics$bird_species_diversity),
    cor_row("diff.entropy-BSD", metrics$diff_entropy, metrics$bird_species_diversity)
  )
  if (isTRUE(write_outputs)) {
    write_csv(metrics, safe_file(tab_dir, "macarthur_metrics_R.csv"))
    write_csv(correlations, safe_file(tab_dir, "macarthur_correlations_R.csv"))
  }
  pv_audit_write_contract_csv(
    metrics %>% transmute(
      site,
      mu,
      sigma2,
      delta,
      FHD3,
      FH_entropy_51,
      gini,
      range,
      sigma,
      CV,
      diff_entropy,
      BSD = bird_species_diversity
    ),
    "EC-026",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "MacArthur foliage metrics"
  )
  pv_audit_write_contract_csv(
    correlations %>% transmute(
      test,
      estimate_r,
      statistic = statistic_t,
      df,
      p_value,
      conf_low,
      conf_high,
      n
    ),
    "EC-027",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "MacArthur correlations"
  )
  p1 <- macarthur_scatter_plot(metrics, "mu", "Mean foliage height", panel_tag = "a")
  p2 <- macarthur_scatter_plot(metrics, "delta", "delta foliage height", panel_tag = "b")
  p3 <- macarthur_scatter_plot(metrics, "diff_entropy", "Differential entropy foliage height", panel_tag = "c")
  scatter_component <- p1 | p2 | p3
  if (isTRUE(write_outputs)) {
    save_png(scatter_component, "macarthur_scatterplots_R_equivalent.png", width = 12, height = 4)
  }
  density_df <- imap_dfr(digitized, ~ mutate(.x, site = .y))
  loess_layers <- pv_audit_macarthur_loess_layers(density_df)
  pv_audit_write_contract_csv(
    loess_layers,
    "EC-028",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "MacArthur foliage density loess layer data"
  )
  pd <- ggplot(density_df, aes(x, y)) + geom_line() + geom_area(alpha = 0.25) + facet_wrap(~site, scales = "free_x") + labs(x="Foliage density", y="Foliage height") + theme_bw()
  if (isTRUE(write_outputs)) {
    save_png(pd, "macarthur_foliage_density_panels_R_equivalent.png", width = 10, height = 8)
    density_component <- macarthur_density_panel_figure(loess_layers)
    save_png(density_component, "macarthur_foliage_density_full_R_equivalent.png", width = 10, height = 8)
    save_png(macarthur_full_figure(scatter_component, density_component), "macarthur_full_R_equivalent.png", width = 12, height = 12)
  }
  list(metrics = metrics, correlations = correlations, digitized = digitized)
}
