# Implements the optional R reimplementation of selected Pellett & Valbuena
# analyses and its provenance helpers. Legacy identifiers retain "replication."

write_model_terms <- function(models, filename) {
  out <- purrr::imap_dfr(models, function(model, name) broom::tidy(model) %>% mutate(model = name, .before = 1))
  readr::write_csv(out, safe_file(tab_dir, filename))
  out
}

write_model_fit <- function(models, filename) {
  out <- purrr::imap_dfr(models, function(model, name) broom::glance(model) %>% mutate(model = name, .before = 1))
  readr::write_csv(out, safe_file(tab_dir, filename))
  out
}

write_table1_rescaled_model_terms <- function(
  allouche_model_df,
  filename = "stage8lnq_table1_rescaled_model_terms_R.csv"
) {
  required_cols <- c("richness", "mu", "delta", "range")
  missing_cols <- setdiff(required_cols, names(allouche_model_df))
  if (length(missing_cols) > 0) {
    stop("Missing Allouche columns for Table 1 rescaled export: ", paste(missing_cols, collapse = ", "))
  }

  model_df <- allouche_model_df
  if (!"mu_rs" %in% names(model_df)) model_df <- model_df %>% mutate(mu_rs = rescale01(mu))
  if (!"delta_rs" %in% names(model_df)) model_df <- model_df %>% mutate(delta_rs = rescale01(delta))
  if (!"range_rs" %in% names(model_df)) model_df <- model_df %>% mutate(range_rs = rescale01(range))
  model_df <- model_df %>%
    mutate(
      richness = as.numeric(richness),
      mu_rs = as.numeric(mu_rs),
      delta_rs = as.numeric(delta_rs),
      range_rs = as.numeric(range_rs)
    )

  term_map <- c(
    "(Intercept)" = "(Intercept)",
    "mu_rs" = "\u03bcrs",
    "I(mu_rs^2)" = "\u03bcrs ^ 2",
    "delta_rs" = "\u03b4rs",
    "I(delta_rs^2)" = "\u03b4rs ^ 2",
    "range_rs" = "rangers",
    "I(range_rs^2)" = "rangers ^ 2"
  )
  display_map <- c(
    "(Intercept)" = "beta0",
    "mu_rs" = "beta_mu",
    "I(mu_rs^2)" = "beta_mu2",
    "delta_rs" = "beta_delta",
    "I(delta_rs^2)" = "beta_delta2",
    "range_rs" = "beta_range",
    "I(range_rs^2)" = "beta_range2"
  )

  specs <- list(
    range_linear = list(
      table1_id = "#1",
      model_label = "#1 richness ~ mu_rs + rangers + rangers^2",
      formula = richness ~ mu_rs + range_rs + I(range_rs^2),
      rescale_basis = "mu_rs=rescale01(mu); rangers=rescale01(range)"
    ),
    range_quadratic = list(
      table1_id = "#2",
      model_label = "#2 richness ~ mu_rs + mu_rs^2 + rangers + rangers^2",
      formula = richness ~ mu_rs + I(mu_rs^2) + range_rs + I(range_rs^2),
      rescale_basis = "mu_rs=rescale01(mu); rangers=rescale01(range)"
    ),
    delta_linear = list(
      table1_id = "#3",
      model_label = "#3 richness ~ mu_rs + delta_rs + delta_rs^2",
      formula = richness ~ mu_rs + delta_rs + I(delta_rs^2),
      rescale_basis = "mu_rs=rescale01(mu); delta_rs=rescale01(delta)"
    ),
    delta_quadratic = list(
      table1_id = "#4",
      model_label = "#4 richness ~ mu_rs + mu_rs^2 + delta_rs + delta_rs^2",
      formula = richness ~ mu_rs + I(mu_rs^2) + delta_rs + I(delta_rs^2),
      rescale_basis = "mu_rs=rescale01(mu); delta_rs=rescale01(delta)"
    ),
    delta_final = list(
      table1_id = "#5",
      model_label = "#5 richness ~ mu_rs + mu_rs^2 + delta_rs",
      formula = richness ~ mu_rs + I(mu_rs^2) + delta_rs,
      rescale_basis = "mu_rs=rescale01(mu); delta_rs=rescale01(delta)"
    )
  )

  out <- purrr::imap_dfr(specs, function(spec, model_id) {
    model <- stats::lm(spec$formula, data = model_df)
    fit <- broom::glance(model)
    broom::tidy(model) %>%
      mutate(
        raw_term = term,
        term = dplyr::coalesce(unname(term_map[raw_term]), raw_term),
        display_term = dplyr::coalesce(unname(display_map[raw_term]), raw_term)
      ) %>%
      transmute(
        model_id = model_id,
        model_label = spec$model_label,
        term,
        display_term,
        estimate,
        std_error = std.error,
        statistic,
        p_value = p.value,
        r_squared = fit$r.squared,
        adjusted_r_squared = fit$adj.r.squared,
        predictor_scale = "rescaled_table1",
        rescale_basis = spec$rescale_basis,
        table1_cell_id = paste0(spec$table1_id, "_", dplyr::row_number(), "_", gsub("[^A-Za-z0-9]+", "_", raw_term)),
        notes = "Table 1-specific rescaled lm export; raw-scale audit outputs and p-values are not altered."
      )
  })
  readr::write_csv(out, safe_file(tab_dir, filename))
  out
}

nested_f_p <- function(null_model, full_model) {
  a <- stats::anova(null_model, full_model)
  as.numeric(a$`Pr(>F)`[2])
}

save_png <- function(plot, filename, width = 8, height = 5) {
  ggplot2::ggsave(safe_file(fig_dir, filename), plot, width = width, height = height, device = "png", dpi = 300)
}

pv_visual_theme <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      strip.background = element_rect(fill = "grey95", colour = "grey55"),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1)
    )
}

pv_visual_schematic_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(size = base_size - 1),
      axis.text = element_text(size = base_size - 2, colour = "grey35"),
      axis.line = element_line(linewidth = 0.25, colour = "grey35"),
      axis.ticks = element_line(linewidth = 0.2, colour = "grey35"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.key.height = grid::unit(0.10, "inches"),
      legend.key.width = grid::unit(0.20, "inches"),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

pv_visual_panel_label <- function(plot, label, size = 7.5) {
  plot +
    labs(tag = label) +
    theme(
      plot.tag = element_text(face = "bold", size = size),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(2, 2, 2, 2)
    )
}

pv_visual_empty_panel <- function(label, subtitle = NULL) {
  subtitle <- if (is.null(subtitle)) "" else subtitle
  ggplot(tibble(x = 0.5, y = 0.5), aes(x, y)) +
    annotate("text", x = 0.5, y = 0.58, label = label, fontface = "bold", size = 3) +
    annotate("text", x = 0.5, y = 0.42, label = subtitle, size = 2.4) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_size = 7)
}

build_figure1_arrow_panel <- function() {
  heterogeneity_colours <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")
  ggplot(tibble(x = 0, y = 0), aes(x, y)) +
    annotate(
      "segment",
      x = 0.16, xend = 0.84, y = 0.34, yend = 0.34,
      arrow = grid::arrow(length = grid::unit(0.060, "inches"), type = "closed"),
      linewidth = 0.36,
      colour = "black"
    ) +
    annotate(
      "segment",
      x = 0.16, xend = 0.84, y = 0.64, yend = 0.64,
      linewidth = 0.26,
      colour = "grey20"
    ) +
    annotate(
      "segment",
      x = c(0.42, 0.42), xend = c(0.16, 0.16), y = c(0.69, 0.59), yend = c(0.69, 0.59),
      arrow = grid::arrow(length = grid::unit(0.052, "inches"), type = "closed"),
      linewidth = 0.34,
      colour = heterogeneity_colours[3:4]
    ) +
    annotate(
      "segment",
      x = c(0.58, 0.58), xend = c(0.84, 0.84), y = c(0.69, 0.59), yend = c(0.69, 0.59),
      arrow = grid::arrow(length = grid::unit(0.052, "inches"), type = "closed"),
      linewidth = 0.34,
      colour = heterogeneity_colours[1:2]
    ) +
    annotate("text", x = 0.50, y = 0.82, label = "Increase in heterogeneity", size = 2.1, colour = "grey20") +
    annotate("text", x = 0.50, y = 0.12, label = "Increase in mean", size = 2.1) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_size = 6)
}

build_figure1_density_schematic <- function(pdf_df) {
  pdf_df %>%
    mutate(mu_lab = factor(mu_lab, levels = paste0("mu = ", sort(unique(mu))))) %>%
    ggplot(aes(x, density)) +
    geom_area(fill = "grey88", colour = NA) +
    geom_line(linewidth = 0.32, colour = "grey15") +
    facet_wrap(~mu_lab, nrow = 1) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(xlim = c(-2, 12), ylim = c(-0.04, 0.65), expand = FALSE) +
    pv_visual_schematic_theme(base_size = 6.2) +
    theme(
      strip.text = element_text(size = 5.8, face = "plain"),
      axis.text.x = element_text(size = 4.8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}

build_figure1_mbh_panel <- function(mbh_df) {
  measure_colours <- c("Coef. var." = "#0072B2", "Gini" = "#E69F00", "Range" = "#009E73", "Variance" = "#CC79A7")
  ggplot(mbh_df, aes(mu, heterogeneity, colour = measure)) +
    geom_line(linewidth = 0.52, lineend = "round") +
    scale_colour_manual(values = measure_colours, breaks = names(measure_colours), drop = FALSE) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    labs(x = "Mean", y = "Heterogeneity", colour = "Heterogeneity measures") +
    pv_visual_schematic_theme(base_size = 6.2) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", override.aes = list(linewidth = 0.65))) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = 5.2),
      legend.text = element_text(size = 4.8),
      legend.key.width = grid::unit(0.22, "inches"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 1, 0)
    )
}

build_figure1_true_relationship_panel <- function(p_c_df) {
  p_c_df %>%
    mutate(relationship = factor(relationship, levels = c("True MDR", "True HDR"))) %>%
    ggplot(aes(x, y)) +
    geom_line(linewidth = 0.42, colour = "grey15") +
    facet_wrap(~relationship, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = "Diversity") +
    pv_visual_schematic_theme(base_size = 6.2) +
    theme(
      strip.text = element_text(size = 5.8),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
}

build_figure1_observed_hdr_panel <- function(p_d_df) {
  measure_colours <- c("Coef. var." = "#0072B2", "Gini" = "#E69F00", "Range" = "#009E73", "Variance" = "#CC79A7")
  ggplot(p_d_df, aes(x, y, colour = measure)) +
    geom_line(linewidth = 0.52, lineend = "round") +
    scale_colour_manual(values = measure_colours, breaks = names(measure_colours), drop = FALSE) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    labs(x = "Observed heterogeneity", y = "Diversity", colour = "Heterogeneity measures") +
    pv_visual_schematic_theme(base_size = 6.2) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", override.aes = list(linewidth = 0.65))) +
    theme(legend.position = "none")
}

assemble_figure1_visual_replication <- function(pdf_df, mbh_df, p_c_df, p_d_df) {
  panel_a <- build_figure1_arrow_panel() / build_figure1_density_schematic(pdf_df) +
    patchwork::plot_layout(heights = c(0.22, 1))
  panel_b <- build_figure1_mbh_panel(mbh_df)
  panel_c <- build_figure1_true_relationship_panel(p_c_df)
  panel_d <- build_figure1_observed_hdr_panel(p_d_df)
  (
    (pv_visual_panel_label(panel_a, "a") | pv_visual_panel_label(panel_b, "b")) /
      (pv_visual_panel_label(panel_c, "c") | pv_visual_panel_label(panel_d, "d"))
  ) +
    patchwork::plot_layout(widths = c(1.05, 1), heights = c(1.12, 0.72)) &
    theme(plot.margin = margin(2, 2, 2, 2))
}

pv_visual_delta_colours <- function(n = 3L) {
  rep_len(c("#0072B2", "#E69F00", "#009E73", "#CC79A7"), n)
}

pv_visual_source_data_path <- function(root, ...) {
  path <- safe_file(root, ...)
  if (!file.exists(path)) {
    return(NA_character_)
  }
  path
}

pv_visual_raster_data <- function(path, max_cells = 250000L, drop_top_fraction = 0) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch({
    r <- terra::rast(path)
    if (terra::nlyr(r) > 1L) {
      r <- r[[1]]
    }
    if (is.finite(drop_top_fraction) && drop_top_fraction > 0 && drop_top_fraction < 1) {
      y_cut <- terra::ymax(r) - (terra::ymax(r) - terra::ymin(r)) * drop_top_fraction
      r <- terra::crop(r, terra::ext(terra::xmin(r), terra::xmax(r), terra::ymin(r), y_cut))
    }
    fact <- max(1L, as.integer(ceiling(sqrt(terra::ncell(r) / max_cells))))
    if (fact > 1L) {
      r <- terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
    }
    df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
    value_col <- setdiff(names(df), c("x", "y"))[1]
    if (is.na(value_col)) {
      return(NULL)
    }
    names(df)[names(df) == value_col] <- "value"
    df <- df[, c("x", "y", "value")]
    df[is.finite(df$value), , drop = FALSE]
  }, error = function(e) NULL)
}

pv_visual_raster_context_panel <- function(
  root,
  path_parts,
  title,
  legend_title,
  palette = "elevation",
  value_limits = NULL,
  drop_top_fraction = 0,
  max_cells = 250000L,
  clip_low = FALSE,
  clip_high = FALSE,
  use_full_tile_extent = FALSE
) {
  path <- do.call(safe_file, c(list(root), as.list(path_parts)))
  raster_df <- pv_visual_raster_data(path, max_cells = max_cells, drop_top_fraction = drop_top_fraction)
  if (is.null(raster_df) || nrow(raster_df) == 0L) {
    return(pv_visual_empty_panel(title, "source raster unavailable"))
  }
  if (!is.null(value_limits) && length(value_limits) == 2L) {
    if (isTRUE(clip_low)) {
      raster_df$value[raster_df$value < value_limits[1]] <- NA_real_
    }
    if (isTRUE(clip_high)) {
      raster_df$value[raster_df$value > value_limits[2]] <- NA_real_
    }
  }
  x_step <- stats::median(diff(sort(unique(raster_df$x))), na.rm = TRUE)
  y_step <- stats::median(diff(sort(unique(raster_df$y))), na.rm = TRUE)
  if (!is.finite(x_step) || x_step <= 0) {
    x_step <- NULL
  }
  if (!is.finite(y_step) || y_step <= 0) {
    y_step <- NULL
  }
  coord_layer <- if (isTRUE(use_full_tile_extent)) {
    x_limits <- range(raster_df$x, finite = TRUE)
    y_limits <- range(raster_df$y, finite = TRUE)
    if (!is.null(x_step)) {
      x_limits <- x_limits + c(-0.5, 0.5) * x_step
    }
    if (!is.null(y_step)) {
      y_limits <- y_limits + c(-0.5, 0.5) * y_step
    }
    coord_equal(xlim = x_limits, ylim = y_limits, expand = FALSE)
  } else {
    coord_equal(expand = FALSE)
  }
  p <- ggplot(raster_df, aes(x, y, fill = value)) +
    geom_tile(width = x_step, height = y_step) +
    coord_layer +
    labs(title = title, x = NULL, y = NULL, fill = legend_title) +
    pv_visual_theme(base_size = 6) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 5.8),
      legend.text = element_text(size = 5.2),
      plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5),
      legend.key.height = grid::unit(0.24, "inches"),
      legend.key.width = grid::unit(0.08, "inches"),
      panel.background = element_rect(fill = "grey96", colour = "grey55", linewidth = 0.18),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25)
    )
  if (identical(palette, "figure2_crop")) {
    p +
      scale_fill_gradientn(
        colours = c("#bfe86b", "#8bd35f", "#52b84f", "#238b45", "#006d2c", "#00441b"),
        limits = value_limits,
        na.value = "grey96"
      ) +
      guides(fill = guide_colourbar(barheight = grid::unit(0.98, "inches"), barwidth = grid::unit(0.10, "inches")))
  } else if (identical(palette, "crop")) {
    p +
      scale_fill_gradientn(
        colours = c("#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#005a32"),
        limits = value_limits,
        na.value = "grey96"
      ) +
      guides(fill = guide_colourbar(barheight = grid::unit(0.98, "inches"), barwidth = grid::unit(0.10, "inches")))
  } else if (identical(palette, "figure2_elevation")) {
    p +
      scale_fill_gradientn(
        colours = c("#6fa35b", "#b6c96b", "#d9c27a", "#c58b45", "#8a5a2b", "#5f4a3a", "#3b342d"),
        limits = value_limits,
        na.value = "grey96"
      ) +
      guides(fill = guide_colourbar(barheight = grid::unit(0.98, "inches"), barwidth = grid::unit(0.10, "inches")))
  } else {
    p +
      scale_fill_gradientn(
        colours = c("#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494", "#8c510a", "#4d2f0c"),
        limits = value_limits,
        na.value = "grey96"
      ) +
      guides(fill = guide_colourbar(barheight = grid::unit(0.98, "inches"), barwidth = grid::unit(0.10, "inches")))
  }
}

figure2_density_panel <- function(density_layers, distribution, title) {
  if (identical(distribution, "gamma")) {
    mean_levels <- c(1.5, 4.0, 7.0, 10.0)
    delta_levels <- c(0.5, 1.5, 3.5)
  } else {
    mean_levels <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    delta_levels <- c(0.03, 0.15, 0.33)
  }
  mean_labels <- paste0("mean = ", mean_levels)
  delta_labels <- paste0("delta = ", delta_levels)
  density_layers %>%
    filter(distribution == !!distribution) %>%
    mutate(
      mean_label = factor(paste0("mean = ", signif(mu, 3)), levels = mean_labels),
      delta_label = factor(paste0("delta = ", signif(delta_or_alpha, 3)), levels = delta_labels)
    ) %>%
    ggplot(aes(x, density, colour = delta_label, group = delta_label)) +
    geom_line(linewidth = 0.35) +
    facet_grid(delta_label ~ mean_label, scales = "free_y", drop = FALSE) +
    scale_colour_manual(values = pv_visual_delta_colours(length(delta_labels)), drop = FALSE) +
    labs(title = title, x = "x", y = "Density", colour = "delta") +
    pv_visual_theme(base_size = 6) +
    theme(legend.position = "none")
}

figure2_measure_specs <- function(kind) {
  if (identical(kind, "gamma")) {
    tibble::tribble(
      ~measure, ~label, ~y_min, ~y_max,
      "Entropy", "Entropy", -0.5, 3.5,
      "Gini", "Gini", -0.05, 1.05,
      "Range", "Range", NA_real_, NA_real_,
      "Variance", "Var.", NA_real_, NA_real_,
      "Std. dev.", "Std. dev.", NA_real_, NA_real_,
      "Coef. var.", "Coef. var.", -0.5, 4.5,
      "delta_L", "delta_L", -0.5, 5.0
    )
  } else {
    tibble::tribble(
      ~measure, ~label, ~y_min, ~y_max,
      "Entropy", "Entropy", -5.0, 0.1,
      "Gini", "Gini", -0.05, 1.05,
      "Range", "Range", -0.05, 1.05,
      "Variance", "Var.", NA_real_, NA_real_,
      "Std. dev.", "Std. dev.", NA_real_, NA_real_,
      "Coef. var.", "Coef. var.", -0.05, 3.5,
      "delta_2", "delta_2", -0.05, 0.5
    )
  }
}

figure2_single_measure_panel <- function(data, measure_name, measure_label, delta_label, y_min, y_max) {
  panel <- data %>%
    filter(measure == !!measure_name) %>%
    ggplot(aes(mu, value, colour = factor(delta), group = delta)) +
    geom_line(linewidth = 0.34) +
    scale_colour_manual(values = pv_visual_delta_colours(length(unique(data$delta))), drop = FALSE) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    labs(x = "Mean", y = measure_label, colour = delta_label) +
    pv_visual_theme(base_size = 5.35) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 4.9),
      axis.text = element_text(size = 4.45),
      plot.margin = margin(0.55, 0.75, 0.55, 0.75)
    )
  if (is.finite(y_min) && is.finite(y_max)) {
    panel <- panel + coord_cartesian(ylim = c(y_min, y_max))
  }
  panel
}

pv_visual_formula_legend_panel <- function(delta_label, formula_label, delta_values) {
  values_numeric <- sort(unique(as.numeric(delta_values)))
  values_numeric <- values_numeric[is.finite(values_numeric)]
  if (!length(values_numeric)) {
    values_numeric <- sort(unique(delta_values))
  }
  colours <- pv_visual_delta_colours(length(values_numeric))
  values <- format(values_numeric, trim = TRUE, scientific = FALSE)
  swatches <- tibble::tibble(
    x0 = rep(0.14, length(values)),
    x1 = rep(0.33, length(values)),
    y = seq(0.62, 0.30, length.out = length(values)),
    label_x = rep(0.39, length(values)),
    label = values,
    colour = colours
  )
  ggplot(tibble(x = 0, y = 0), aes(x, y)) +
    geom_segment(
      data = swatches,
      aes(x = x0, xend = x1, y = y, yend = y, colour = colour),
      inherit.aes = FALSE,
      linewidth = 0.70
    ) +
    geom_text(
      data = swatches,
      aes(x = label_x, y = y, label = label),
      inherit.aes = FALSE,
      size = 1.75,
      hjust = 0
    ) +
    annotate("text", x = 0.06, y = 0.86, label = formula_label, parse = TRUE, size = 1.92, hjust = 0) +
    scale_colour_identity() +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_size = 5.6) +
    theme(plot.margin = margin(0.15, 0.20, 0.15, 0.20))
}

figure2_measure_axis <- function(data, title, delta_label, panel_labels = NULL) {
  kind <- if (identical(delta_label, "delta_L")) "gamma" else "beta"
  specs <- figure2_measure_specs(kind)
  plots <- purrr::pmap(specs, function(measure, label, y_min, y_max) {
    figure2_single_measure_panel(data, measure, label, delta_label, y_min, y_max)
  })
  if (!is.null(panel_labels)) {
    plots <- purrr::map2(plots, panel_labels[seq_along(plots)], pv_visual_panel_label)
  }
  formula_label <- if (identical(kind, "gamma")) {
    "delta[L] == frac(sigma^2, mu - L)"
  } else {
    "delta[2] == frac(sigma^2, (mu - L) * (U - mu))"
  }
  legend <- pv_visual_formula_legend_panel(delta_label, formula_label, sort(unique(data$delta)))
  top <- patchwork::wrap_plots(plots[1:6], ncol = 3)
  bottom <- patchwork::wrap_plots(c(plots[7], list(legend), list(patchwork::plot_spacer())), ncol = 3,
                                  widths = c(0.86, 0.88, 1.26))
  (top / bottom) +
    patchwork::plot_layout(heights = c(2.05, 1.00)) +
    patchwork::plot_annotation(title = title) &
    theme(plot.title = element_text(size = 5.9, face = "bold", hjust = 0.5),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))
}

figure2_empirical_context_panel <- function(data, title, x_col, y_col) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.35, size = 0.35, colour = "grey25") +
    labs(title = title, x = x_col, y = y_col) +
    pv_visual_theme(base_size = 6) +
    theme(legend.position = "none")
}

figure2_empirical_panel_specs <- function(kind) {
  if (identical(kind, "crop")) {
    tibble::tribble(
      ~measure, ~label, ~bins, ~threshold, ~colour_min, ~colour_max, ~x_min, ~x_max, ~y_min, ~y_max, ~x_breaks, ~y_breaks, ~count_ticks,
      "Variance", "Variance (e+3)", list(100L), 1, 1, 200, 0, 100, 0, 2.4, list(c(0, 50, 100)), list(c(0, 1, 2)), list(c(1, 10, 100)),
      "Std. dev.", "Std. Dev.", list(100L), 1, 1, 200, 0, 100, 0, 50, list(c(0, 50, 100)), list(c(0, 20, 40)), list(c(1, 10, 100)),
      "Coef. var.", "Coef. Var.", list(120L), 1, 1, 200, 0, 100, -0.5, 9, list(c(0, 50, 100)), list(c(0, 4, 8)), list(c(1, 10, 100)),
      "delta_2", "delta_2", list(100L), 1, 1, 200, 0, 100, -0.05, 1.05, list(c(0, 50, 100)), list(c(0, 0.5, 1)), list(c(1, 10, 100))
    )
  } else {
    tibble::tribble(
      ~measure, ~label, ~bins, ~threshold, ~colour_min, ~colour_max, ~x_min, ~x_max, ~y_min, ~y_max, ~x_breaks, ~y_breaks, ~count_ticks,
      "Variance", "Variance (e+4)", list(c(200L, 300L)), 4, 4, 280, 0, 4000, 0, 9, list(c(0, 1500, 3000)), list(c(0, 4, 8)), list(c(4, 15, 150)),
      "Std. dev.", "Std. Dev. (e+2)", list(c(150L, 150L)), 4, 4, 280, 0, 4000, 0, 3.5, list(c(0, 1500, 3000)), list(c(0, 1, 2, 3)), list(c(4, 15, 150)),
      "Coef. var.", "Coef. Var.", list(c(150L, 150L)), 4, 4, 280, 0, 4000, 0, 0.7, list(c(0, 1500, 3000)), list(c(0, 0.3, 0.6)), list(c(4, 15, 150)),
      "delta_L", "delta_L", list(c(200L, 300L)), 4, 4, 280, 0, 4000, 0, 30, list(c(0, 1500, 3000)), list(c(0, 10, 20, 30)), list(c(4, 15, 150))
    )
  }
}

figure2_empirical_layout_spec <- function(kind) {
  if (identical(kind, "crop")) {
    list(
      grid_heights = c(1.0, 0.72),
      guide_widths = c(1.0, 0.31, 0.04),
      context_widths = c(1.34, 1.06),
      context_heights = c(0.96, 0.90),
      title_size = 5.8,
      map_legend_height = 0.98,
      map_legend_width = 0.10
    )
  } else {
    list(
      grid_heights = c(1.0, 0.72),
      guide_widths = c(1.0, 0.31, 0.04),
      context_widths = c(1.34, 1.06),
      context_heights = c(0.96, 0.90),
      title_size = 5.8,
      map_legend_height = 0.98,
      map_legend_width = 0.10
    )
  }
}

figure2_empirical_axis_label_spec <- function(kind, measure_name, x_label = NULL) {
  spec <- figure2_empirical_panel_specs(kind) %>% filter(measure == !!measure_name)
  y_label <- if (nrow(spec)) spec$label[[1]] else measure_name
  default_x <- if (identical(kind, "crop")) "Mean cover (%)" else "Mean elev. (m)"
  list(
    x_label = if (is.null(x_label) || !nzchar(x_label)) default_x else x_label,
    y_label = y_label
  )
}

figure2_empirical_kind <- function(data) {
  if ("delta_2" %in% unique(data$measure)) "crop" else "elevation"
}

pv_visual_compact_count_guide <- function(kind) {
  spec <- figure2_empirical_panel_specs(kind)
  value_range <- range(c(spec$colour_min, spec$colour_max), finite = TRUE)
  value_range <- as.numeric(value_range)
  if (length(value_range) < 2L || !all(is.finite(value_range)) || any(value_range <= 0)) {
    value_range <- c(1, 100)
  }
  ticks <- pv_visual_axis_breaks(spec$count_ticks[[1]])
  ticks <- ticks[is.finite(ticks) & ticks > 0]
  if (!length(ticks)) {
    ticks <- pretty(value_range, n = 3)
    ticks <- ticks[is.finite(ticks) & ticks > 0]
  }
  if (!length(ticks)) {
    ticks <- value_range
  }
  ticks <- sort(unique(as.numeric(ticks)))
  vals <- exp(seq(log(value_range[1]), log(value_range[2]), length.out = 80))
  guide_df <- tibble::tibble(x = 1, y = seq_along(vals), value = vals)
  ggplot(guide_df, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(trans = "log10", limits = value_range, breaks = ticks, na.value = "grey96") +
    annotate("text", x = 0.46, y = length(vals) + 4.6, label = "N obs.", size = 1.60, hjust = 0) +
    annotate("text", x = 1.48, y = approx(log(vals), seq_along(vals), xout = log(ticks), rule = 2)$y,
             label = ticks, size = 1.38, hjust = 0) +
    coord_cartesian(xlim = c(0.38, 2.05), ylim = c(0, length(vals) + 6.5), expand = FALSE) +
    theme_void(base_size = 5.0) +
    theme(legend.position = "none", plot.margin = margin(0.25, 0.25, 0.25, 0.25))
}


pv_visual_numeric_values <- function(x) {
  # Recursively extract numeric values from list or tuple specifications used by
  # Julia and Makie without changing empirical source data. This is used only
  # to give ggplot2 hexbin display parameters the required type.
  if (is.null(x)) {
    return(numeric())
  }
  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }
  if (is.list(x) || is.pairlist(x)) {
    vals <- unlist(lapply(x, pv_visual_numeric_values), recursive = FALSE, use.names = FALSE)
    return(suppressWarnings(as.numeric(vals)))
  }
  suppressWarnings(tryCatch(as.numeric(x), error = function(e) numeric()))
}

pv_visual_normalize_hex_bins <- function(bins, default = 35L) {
  # This guard affects only visual output from ggplot2::geom_hex/stat_binhex.
  # Bin specifications from Julia and Makie may arrive as nested lists or tuples
  # containing values for two axes. ggplot2 needs a finite positive numeric scalar here.
  vals <- pv_visual_numeric_values(bins)
  vals <- vals[is.finite(vals) & vals > 0]
  if (!length(vals)) {
    return(as.integer(default))
  }
  as.integer(round(vals[[1]]))
}

pv_visual_axis_breaks <- function(x) {
  # This guard affects only breaks in ggplot2 continuous scales.
  # Source axis specifications from Julia may be lists or tuples; ggplot2 scale breaks
  # must be an atomic numeric vector or NULL.
  vals <- pv_visual_numeric_values(x)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(NULL)
  }
  sort(unique(as.numeric(vals)))
}

pv_visual_axis_limits <- function(x) {
  # This guard affects only limits in ggplot2 continuous scales.
  # Return a numeric range of two values or NULL. This preserves values from the
  # source while preventing lists from being passed as scale parameters.
  vals <- pv_visual_numeric_values(x)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2L) {
    return(NULL)
  }
  as.numeric(c(min(vals), max(vals)))
}

figure2_empirical_measure_panel <- function(data, measure_name, x_label, kind = NULL, show_legend = FALSE) {
  if (is.null(kind)) {
    kind <- figure2_empirical_kind(data)
  }
  spec <- figure2_empirical_panel_specs(kind) %>% filter(measure == !!measure_name)
  if (nrow(spec) == 0L) {
    spec <- tibble::tibble(
      label = measure_name, bins = list(34L), threshold = 1,
      colour_min = 1, colour_max = NA_real_, x_min = NA_real_, x_max = NA_real_,
      y_min = NA_real_, y_max = NA_real_, x_breaks = list(NULL), y_breaks = list(NULL),
      count_ticks = list(NULL)
    )
  }
  threshold <- spec$threshold[[1]]
  colour_limits <- c(spec$colour_min[[1]], spec$colour_max[[1]])
  axis_labels <- figure2_empirical_axis_label_spec(kind, measure_name, x_label)
  panel <- data %>%
    filter(measure == !!measure_name) %>%
    ggplot(aes(mean, value)) +
    geom_hex(
      bins = pv_visual_normalize_hex_bins(spec$bins[[1]]),
      linewidth = 0,
      aes(fill = after_stat(ifelse(count >= threshold, count, NA_real_)))
    ) +
    scale_fill_viridis_c(
      trans = "log10",
      limits = if (all(is.finite(colour_limits))) colour_limits else NULL,
      breaks = pv_visual_axis_breaks(spec$count_ticks[[1]]),
      na.value = "white"
    ) +
    scale_x_continuous(breaks = pv_visual_axis_breaks(spec$x_breaks[[1]])) +
    scale_y_continuous(breaks = pv_visual_axis_breaks(spec$y_breaks[[1]])) +
    coord_cartesian(
      xlim = c(spec$x_min[[1]], spec$x_max[[1]]),
      ylim = c(spec$y_min[[1]], spec$y_max[[1]]),
      expand = FALSE
    ) +
    labs(x = axis_labels$x_label, y = axis_labels$y_label, fill = "N obs.") +
    pv_visual_theme(base_size = 5.15) +
    theme(
      legend.position = if (isTRUE(show_legend)) "bottom" else "none",
      legend.title = element_text(size = 4.7),
      legend.text = element_text(size = 4.25),
      axis.title = element_text(size = 4.65),
      axis.text = element_text(size = 4.2),
      plot.margin = margin(0.45, 0.45, 0.45, 0.45)
    )
  if (isTRUE(show_legend)) {
    panel <- panel +
      guides(fill = guide_colourbar(
        barheight = grid::unit(0.08, "inches"),
        barwidth = grid::unit(0.42, "inches"),
        title.position = "top"
      ))
  }
  panel
}

figure2_empirical_measure_grid <- function(data, title, x_label = "Mean", panel_labels = NULL, kind = NULL) {
  if (is.null(kind)) {
    kind <- figure2_empirical_kind(data)
  }
  layout <- figure2_empirical_layout_spec(kind)
  measure_order <- c("Variance", "Std. dev.", "Coef. var.", "delta_L", "delta_2")
  measures <- intersect(measure_order, unique(data$measure))
  plots <- purrr::map(measures, ~ figure2_empirical_measure_panel(data, .x, x_label, kind = kind))
  if (!is.null(panel_labels)) {
    plots <- purrr::map2(plots, panel_labels[seq_along(plots)], pv_visual_panel_label)
  }
  if (length(plots) >= 4L) {
    panel_grid <- patchwork::wrap_plots(plots[1:4], ncol = 2)
    grid_with_guide <- patchwork::wrap_plots(
      patchwork::wrap_elements(full = panel_grid),
      pv_visual_compact_count_guide(kind),
      ncol = 2,
      widths = c(1, 0.14)
    )
  } else {
    grid_with_guide <- patchwork::wrap_plots(plots, ncol = min(2, length(plots)))
  }
  grid_with_guide +
    patchwork::plot_annotation(title = title) &
    theme(plot.title = element_text(size = layout$title_size, face = "bold", hjust = 0.5),
          legend.position = "none",
          plot.margin = margin(0.4, 0.4, 0.4, 0.4))
}

figure2_empirical_context_block <- function(data, title, x_label, context_panel, panel_labels, context_label, kind = NULL) {
  if (is.null(kind)) {
    kind <- figure2_empirical_kind(data)
  }
  layout <- figure2_empirical_layout_spec(kind)
  measure_grid <- figure2_empirical_measure_grid(data, title, x_label, panel_labels = panel_labels, kind = kind)
  (
    patchwork::wrap_elements(full = measure_grid) /
      pv_visual_panel_label(context_panel, context_label)
  ) +
    patchwork::plot_layout(heights = layout$context_heights)
}

figure2_srtm_context_panel <- function(root = repo_root) {
  layout <- figure2_empirical_layout_spec("elevation")
  pv_visual_raster_context_panel(
    root,
    c("02_MBH_comb", "data", "elevation", "lowres_full_world.tif"),
    "SRTM elevation context",
    "Elevation",
    palette = "figure2_elevation",
    value_limits = c(0, 4000),
    max_cells = 260000L,
    clip_low = TRUE,
    use_full_tile_extent = TRUE
  ) +
    labs(title = NULL, fill = "Elevation (m)") +
    guides(fill = guide_colourbar(
      barheight = grid::unit(layout$map_legend_height, "inches"),
      barwidth = grid::unit(layout$map_legend_width, "inches")
    )) +
    theme(
      plot.title = element_blank(),
      legend.title = element_text(size = 4.6),
      legend.text = element_text(size = 4.15),
      legend.key.height = grid::unit(layout$map_legend_height, "inches"),
      legend.key.width = grid::unit(layout$map_legend_width, "inches"),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25)
    )
}

figure2_crop_context_panel <- function(root = repo_root) {
  layout <- figure2_empirical_layout_spec("crop")
  pv_visual_raster_context_panel(
    root,
    c("02_MBH_comb", "data", "crop_cover", "lowres_crop.tif"),
    "Crop-cover context",
    "Crop cover (%)",
    palette = "figure2_crop",
    value_limits = c(0, 100),
    max_cells = 260000L,
    clip_high = TRUE,
    use_full_tile_extent = TRUE
  ) +
    labs(title = NULL, fill = "Crop cover (%)") +
    guides(fill = guide_colourbar(
      barheight = grid::unit(layout$map_legend_height, "inches"),
      barwidth = grid::unit(layout$map_legend_width, "inches")
    )) +
    theme(
      plot.title = element_blank(),
      legend.title = element_text(size = 4.6),
      legend.text = element_text(size = 4.15),
      legend.key.height = grid::unit(layout$map_legend_height, "inches"),
      legend.key.width = grid::unit(layout$map_legend_width, "inches"),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25)
    )
}

assemble_figure2_visual_replication <- function(gamma_df, beta_df, density_layers, elev_long, crop_long, elev, crop, root = repo_root) {
  gamma_block <- pv_visual_panel_label(figure2_density_panel(density_layers, "gamma", "Lower-bounded variables"), "a") /
    figure2_measure_axis(gamma_df, "Gamma theoretical measures", "delta_L", panel_labels = letters[2:8]) +
    patchwork::plot_layout(heights = c(0.8, 1.55))
  beta_block <- pv_visual_panel_label(figure2_density_panel(density_layers, "beta", "Double-bounded variables"), "i") /
    figure2_measure_axis(beta_df, "Beta theoretical measures", "delta_2", panel_labels = letters[10:16]) +
    patchwork::plot_layout(heights = c(0.8, 1.55))
  elev_block <- figure2_empirical_context_block(
    elev_long,
    "Land elevation heterogeneity from the SRTM",
    "Mean elev. (m)",
    figure2_srtm_context_panel(root),
    panel_labels = letters[17:20],
    context_label = "u",
    kind = "elevation"
  )
  crop_block <- figure2_empirical_context_block(
    crop_long,
    "Crop cover heterogeneity from ESA PROBA-V",
    "Mean crop cover (%)",
    figure2_crop_context_panel(root),
    panel_labels = letters[22:25],
    context_label = "z",
    kind = "crop"
  )
  ((gamma_block | beta_block) / (elev_block | crop_block)) +
    patchwork::plot_layout(heights = c(1.20, 1.24), widths = c(1, 1)) &
    theme(plot.margin = margin(1.2, 1.2, 1.2, 1.2))
}

pv_even_integer_indices <- function(n_rows, max_n) {
  n_rows <- as.integer(n_rows)
  max_n <- as.integer(max_n)
  if (is.na(n_rows) || n_rows <= 0L || is.na(max_n) || max_n <= 0L) {
    return(integer(0))
  }
  if (n_rows <= max_n) {
    return(seq_len(n_rows))
  }
  unique(as.integer(round(seq(1L, n_rows, length.out = max_n))))
}

figure3_sample_grid <- function(data, max_n = 20000) {
  rows <- data %>%
    filter(if_all(everything(), ~ !is.nan(.x)))
  rows %>%
    slice(pv_even_integer_indices(nrow(rows), max_n)) %>%
    as_tibble()
}

figure3_selected_levels <- function(x, n = 7L) {
  vals <- sort(unique(x[is.finite(x)]))
  if (length(vals) <= n) {
    return(vals)
  }
  vals[unique(as.integer(round(seq(1, length(vals), length.out = n))))]
}

figure3_source_levels <- function() {
  list(
    mu_range = c(64, 2422),
    H_range = c(0.64, 144),
    mu_lines = seq(64, 2422, length.out = 7L),
    H_lines = seq(0.64, 144, length.out = 7L)
  )
}

figure3_source_line_targets <- function(col) {
  source_levels <- figure3_source_levels()
  if (identical(col, "mu")) {
    return(source_levels$mu_lines)
  }
  if (identical(col, "H")) {
    return(source_levels$H_lines)
  }
  NULL
}

figure3_source_display_spec <- function() {
  source_levels <- figure3_source_levels()
  list(
    line_family = list(
      H_colour = "#0072B2",
      mu_colour = "#E69F00",
      other_colour = "grey35",
      alpha = seq(0.22, 1, length.out = 7L),
      H_linetype = "solid",
      mu_linetype = "22",
      outline_width = 0.78,
      underlay_width = 0.56,
      inner_width = 0.34,
      caption = "solid: heterogeneity; dashed: mean"
    ),
    theory = list(
      base_size = 5.45,
      axis_text_size = 4.2,
      title_size = 5.55,
      model_width = 0.5,
      row_heights = c(0.74, 1, 1, 1, 1),
      panel_margin = c(0.55, 0.55, 0.55, 0.55)
    ),
    empirical = list(
      base_size = 5.35,
      point_outer_size = 0.95,
      point_inner_size = 0.34,
      point_alpha = 0.82,
      fit_outer_width = 0.78,
      fit_inner_width = 0.46,
      fit_inner_colour = "grey82",
      row_heights = c(0.74, 1, 1, 1, 1),
      map_width = 0.52,
      panel_margin = c(0.55, 0.55, 0.55, 0.55)
    ),
    context_map = list(
      title = "Catalonia",
      legend_title = "Elevation (km)",
      value_limits = c(0, 2500),
      max_cells = 220000L,
      legend_height = 0.5,
      legend_width = 0.055
    ),
    layout = list(
      theory_width = 1.5,
      empirical_width = 1,
      figure_margin = c(1.2, 1.2, 1.2, 1.2)
    ),
    source_levels = source_levels
  )
}

pv_visual_nearest_levels <- function(values, targets) {
  vals <- sort(unique(values[is.finite(values)]))
  if (!length(vals)) {
    return(vals)
  }
  if (is.null(targets) || !length(targets)) {
    return(figure3_selected_levels(vals, n = 7L))
  }
  vals[unique(vapply(targets, function(target) {
    which.min(abs(vals - target))
  }, integer(1)))]
}

figure3_visual_theory_grid_julia_scaled <- function(n = 500L) {
  # This visual grid matches the density and measure display units used by Julia
  # without changing audit exports.
  theory_grid_R(n = n) %>%
    mutate(
      range_mbh = range_mbh / 100,
      cv_mbh = cv_mbh * 10,
      var_mbh = var_mbh / 10000,
      delta_mbh = delta_mbh / 10
    )
}

figure3_source_axis_limits <- function(col) {
  source_levels <- figure3_source_levels()
  if (identical(col, "mu")) {
    return(source_levels$mu_range)
  }
  if (identical(col, "H")) {
    return(source_levels$H_range)
  }
  NULL
}

figure3_line_family_data <- function(
  data,
  x_col,
  y_col,
  family_col,
  linetype = "solid",
  spec = figure3_source_display_spec(),
  role = NULL
) {
  family_levels <- pv_visual_nearest_levels(data[[family_col]], figure3_source_line_targets(family_col))
  if (is.null(role)) {
    role <- family_col
  }
  base_colour <- if (identical(role, "H")) {
    spec$line_family$H_colour
  } else if (identical(role, "mu")) {
    spec$line_family$mu_colour
  } else {
    spec$line_family$other_colour
  }
  alpha_values <- spec$line_family$alpha[seq_along(family_levels)]
  level_colours <- vapply(
    alpha_values,
    function(alpha_value) grDevices::adjustcolor(base_colour, alpha.f = alpha_value),
    character(1)
  )
  names(level_colours) <- as.character(family_levels)
  data %>%
    filter(.data[[family_col]] %in% family_levels) %>%
    arrange(.data[[family_col]], .data[[x_col]], .data[[y_col]]) %>%
    mutate(
      .line_group = factor(paste(role, .data[[family_col]], sep = "_")),
      .line_linetype = linetype,
      .line_role = role,
      .line_colour = unname(level_colours[as.character(.data[[family_col]])])
    )
}

figure3_dual_line_family_layers <- function(
  data,
  x_col,
  y_col,
  family_col,
  secondary_family_col = NULL,
  spec = figure3_source_display_spec()
) {
  family_role <- if (identical(family_col, "mu")) "mu" else if (identical(family_col, "H")) "H" else "other"
  family_linetype <- if (identical(family_role, "mu")) spec$line_family$mu_linetype else spec$line_family$H_linetype
  line_data <- figure3_line_family_data(
    data,
    x_col,
    y_col,
    family_col,
    linetype = family_linetype,
    spec = spec,
    role = family_role
  )
  layers <- list(
    geom_line(
      data = line_data,
      aes(x = .data[[x_col]], y = .data[[y_col]], group = .line_group),
      linewidth = spec$line_family$outline_width,
      colour = "black",
      linetype = family_linetype,
      lineend = "round"
    ),
    geom_line(
      data = line_data,
      aes(x = .data[[x_col]], y = .data[[y_col]], group = .line_group),
      linewidth = spec$line_family$underlay_width,
      colour = "white",
      linetype = family_linetype,
      lineend = "round"
    ),
    geom_line(
      data = line_data,
      aes(x = .data[[x_col]], y = .data[[y_col]], group = .line_group, colour = .line_colour),
      linewidth = spec$line_family$inner_width,
      linetype = family_linetype,
      lineend = "round"
    )
  )
  if (!is.null(secondary_family_col)) {
    secondary_role <- if (identical(secondary_family_col, "mu")) "mu" else if (identical(secondary_family_col, "H")) "H" else "other"
    secondary_linetype <- if (identical(secondary_role, "mu")) spec$line_family$mu_linetype else spec$line_family$H_linetype
    secondary_data <- figure3_line_family_data(
      data,
      x_col,
      y_col,
      secondary_family_col,
      linetype = secondary_linetype,
      spec = spec,
      role = secondary_role
    )
    layers <- c(
      layers,
      list(
        geom_line(
          data = secondary_data,
          aes(x = .data[[x_col]], y = .data[[y_col]], group = .line_group),
          linewidth = spec$line_family$outline_width,
          colour = "black",
          linetype = secondary_linetype,
          lineend = "round"
        ),
        geom_line(
          data = secondary_data,
          aes(x = .data[[x_col]], y = .data[[y_col]], group = .line_group),
          linewidth = spec$line_family$underlay_width,
          colour = "white",
          linetype = secondary_linetype,
          lineend = "round"
        ),
        geom_line(
          data = secondary_data,
          aes(x = .data[[x_col]], y = .data[[y_col]], group = .line_group, colour = .line_colour),
          linewidth = spec$line_family$inner_width,
          linetype = secondary_linetype,
          lineend = "round"
        )
      )
    )
  }
  layers
}

figure3_compact_model_inset <- function(model, measure_label, spec = figure3_source_display_spec()) {
  co <- stats::coef(model)
  model_x <- model$model$x
  x_range <- range(model_x, finite = TRUE)
  x_seq <- seq(x_range[1], x_range[2], length.out = 120L)
  curve_df <- tryCatch(
    tibble::tibble(x = x_seq, y = as.numeric(stats::predict(model, newdata = data.frame(x = x_seq)))),
    error = function(e) tibble::tibble(x = numeric(), y = numeric())
  )
  label <- sprintf(
    "%s\n%.1f %+.2g x %+.2g x^2",
    measure_label,
    unname(co[1]),
    unname(co[2]),
    unname(co[3])
  )
  if (!nrow(curve_df)) {
    return(pv_visual_empty_panel("Quadratic fit", measure_label))
  }
  ggplot(curve_df, aes(x, y)) +
    geom_line(linewidth = 0.44, colour = "grey15") +
    annotate(
      "label",
      x = mean(x_range),
      y = max(curve_df$y, na.rm = TRUE),
      label = label,
      size = 1.18,
      hjust = 0.5,
      vjust = 1,
      label.size = 0,
      fill = "grey92",
      label.padding = grid::unit(0.06, "lines")
    ) +
    labs(x = NULL, y = NULL) +
    pv_visual_theme(base_size = spec$theory$base_size) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 3.75),
      axis.ticks = element_line(linewidth = 0.1),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.18, colour = "grey35"),
      plot.margin = margin(0.45, 0.45, 0.45, 0.45)
    )
}

figure3_empirical_style_spec <- function() {
  figure3_source_display_spec()$empirical
}

figure3_context_map_style_spec <- function() {
  figure3_source_display_spec()$context_map
}

figure3_theory_line_family_panel <- function(
  data,
  x_col,
  y_col,
  family_col,
  x_label,
  y_label,
  colour_label,
  title,
  secondary_family_col = NULL,
  show_legend = TRUE
) {
  spec <- figure3_source_display_spec()
  panel <- ggplot() +
    figure3_dual_line_family_layers(data, x_col, y_col, family_col, secondary_family_col, spec = spec) +
    scale_colour_identity() +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      caption = if (isTRUE(show_legend)) spec$line_family$caption else NULL
    ) +
    pv_visual_theme(base_size = spec$theory$base_size) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = spec$theory$base_size - 0.25),
      axis.text = element_text(size = spec$theory$axis_text_size),
      plot.title = element_text(size = spec$theory$title_size, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 3.7, hjust = 0.5, colour = "grey35"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.08, colour = "grey88"),
      plot.margin = do.call(margin, as.list(spec$theory$panel_margin))
    )
  x_limits <- figure3_source_axis_limits(x_col)
  y_limits <- figure3_source_axis_limits(y_col)
  if (!is.null(x_limits) || !is.null(y_limits)) {
    panel <- panel + coord_cartesian(xlim = x_limits, ylim = y_limits)
  }
  panel
}

figure3_equation_inset_panel <- function(model, measure_label) {
  figure3_compact_model_inset(model, measure_label)
}

figure3_model_label_panel <- function(model, measure_label) {
  figure3_equation_inset_panel(model, measure_label)
}

figure3_theory_observed_hdr_line_family_panel <- function(data, measure_col, measure_label) {
  spec <- figure3_source_display_spec()
  figure3_theory_line_family_panel(
    data,
    measure_col,
    "richness",
    "mu",
    measure_label,
    "Richness",
    "Mean",
    "Observed HDR",
    secondary_family_col = "H"
  ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = spec$theory$title_size, face = "bold", hjust = 0.5),
      plot.margin = do.call(margin, as.list(spec$theory$panel_margin))
    )
}

figure3_theory_mdr_hdr_row <- function(tg) {
  p_mdr <- figure3_theory_line_family_panel(
    tg,
    "mu",
    "richness",
    "H",
    "Mean",
    "Richness",
    "Heterogeneity",
    "Mean-Diversity (MDR)",
    secondary_family_col = "mu"
  )
  p_hdr <- figure3_theory_line_family_panel(
    tg,
    "H",
    "richness",
    "mu",
    "Heterogeneity",
    "Richness",
    "Mean",
    "Heterogeneity-Diversity (HDR)",
    secondary_family_col = "H"
  )
  pv_visual_panel_label(p_mdr, "a") | pv_visual_panel_label(p_hdr, "b")
}

figure3_theory_mbh_hdr_row <- function(tg, measure_col, measure_label, panel_labels = NULL) {
  model <- stats::lm(richness ~ x + I(x^2), data = data.frame(x = tg[[measure_col]], richness = tg$richness))
  p_mhr <- figure3_theory_line_family_panel(
    tg,
    "mu",
    measure_col,
    "H",
    "Mean",
    measure_label,
    "Heterogeneity",
    "Mean-heterogeneity",
    secondary_family_col = "mu"
  )
  p_hdr <- figure3_theory_observed_hdr_line_family_panel(tg, measure_col, measure_label)
  p_model <- figure3_model_label_panel(model, measure_label)
  row <- p_mhr | p_hdr | p_model
  if (!is.null(panel_labels) && length(panel_labels) == 3L) {
    row <- pv_visual_panel_label(p_mhr, panel_labels[1]) |
      pv_visual_panel_label(p_hdr, panel_labels[2]) |
      pv_visual_panel_label(p_model, panel_labels[3])
  }
  row + patchwork::plot_layout(widths = c(1, 1, figure3_source_display_spec()$theory$model_width))
}

figure3_empirical_mdr_row <- function(allouche_model_df) {
  spec <- figure3_empirical_style_spec()
  ggplot(allouche_model_df, aes(mu, richness)) +
    geom_point(size = spec$point_outer_size, alpha = spec$point_alpha, colour = "black") +
    geom_point(size = spec$point_inner_size, alpha = 0.9, colour = "white") +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, linewidth = spec$fit_outer_width, colour = "black") +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, linewidth = spec$fit_inner_width, colour = spec$fit_inner_colour) +
    labs(title = "Mean-Diversity (MDR)", x = "Mean elevation", y = "Richness") +
    pv_visual_theme(base_size = spec$base_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = spec$base_size + 0.2, face = "bold", hjust = 0.5),
      plot.margin = do.call(margin, as.list(spec$panel_margin))
    )
}

figure3_catalonia_context_panel <- function(root = repo_root) {
  spec <- figure3_context_map_style_spec()
  pv_visual_raster_context_panel(
    root,
    c("03_corrected_HDR", "empirical", "01_carnicer", "data", "elevation", "catalunya_lowres.tif"),
    spec$title,
    spec$legend_title,
    palette = "elevation",
    value_limits = spec$value_limits,
    max_cells = spec$max_cells
  ) +
    guides(fill = guide_colourbar(
      barheight = grid::unit(spec$legend_height, "inches"),
      barwidth = grid::unit(spec$legend_width, "inches")
    )) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 4.7),
      legend.text = element_text(size = 4.1),
      plot.title = element_text(size = 5.2, face = "bold", hjust = 0.5),
      plot.margin = margin(0.4, 0.4, 0.4, 0.4)
    )
}

figure3_empirical_map_panel <- function(allouche_model_df, root = repo_root) {
  spec <- figure3_context_map_style_spec()
  raster_panel <- figure3_catalonia_context_panel(root)
  if (!inherits(raster_panel, "ggplot") && all(c("X_coord", "Y_coord") %in% names(allouche_model_df))) {
    return(
      ggplot(allouche_model_df, aes(X_coord, Y_coord, colour = richness)) +
        geom_point(size = 0.55, alpha = 0.85) +
        coord_equal() +
        scale_colour_viridis_c() +
        labs(title = spec$title, x = NULL, y = NULL, colour = "Richness") +
        pv_visual_theme(base_size = 5.4) +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          plot.margin = margin(0.4, 0.4, 0.4, 0.4)
        )
    )
  }
  raster_panel
}

figure3_empirical_hdr_row <- function(emp_layers, measure_name, measure_label, panel_labels = NULL) {
  spec <- figure3_empirical_style_spec()
  layer <- emp_layers %>%
    filter(measure == !!measure_name) %>%
    arrange(x)
  fit_layer <- layer %>%
    filter(is.finite(fit_y)) %>%
    arrange(x)
  p_mhr <- ggplot(layer, aes(raw_mu, x)) +
    geom_point(size = spec$point_outer_size, alpha = spec$point_alpha, colour = "black") +
    geom_point(size = spec$point_inner_size, alpha = 0.9, colour = "white") +
    labs(title = "Mean-heterogeneity", x = "Mean elevation", y = measure_label) +
    pv_visual_theme(base_size = spec$base_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = spec$base_size + 0.2, face = "bold", hjust = 0.5),
      plot.margin = do.call(margin, as.list(spec$panel_margin))
    )
  p_hdr <- ggplot(layer, aes(x, richness)) +
    geom_point(size = spec$point_outer_size, alpha = spec$point_alpha, colour = "black") +
    geom_point(size = spec$point_inner_size, alpha = 0.9, colour = "white") +
    geom_line(data = fit_layer, aes(x, fit_y), linewidth = spec$fit_outer_width, colour = "black") +
    geom_line(data = fit_layer, aes(x, fit_y), linewidth = spec$fit_inner_width, colour = spec$fit_inner_colour) +
    labs(title = "Empirical HDR", x = measure_label, y = "Richness") +
    pv_visual_theme(base_size = spec$base_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = spec$base_size + 0.2, face = "bold", hjust = 0.5),
      plot.margin = do.call(margin, as.list(spec$panel_margin))
    )
  row <- p_mhr | p_hdr
  if (!is.null(panel_labels) && length(panel_labels) == 2L) {
    row <- pv_visual_panel_label(p_mhr, panel_labels[1]) | pv_visual_panel_label(p_hdr, panel_labels[2])
  }
  row
}

assemble_figure3_visual_replication <- function(tg, allouche_model_df, empirical_layers, root = repo_root) {
  spec <- figure3_source_display_spec()
  tg_visual <- figure3_visual_theory_grid_julia_scaled(n = 500L)
  theory <- (
    figure3_theory_mdr_hdr_row(tg_visual) /
      figure3_theory_mbh_hdr_row(tg_visual, "range_mbh", "Range (e+2)", panel_labels = c("c", "d", "e")) /
      figure3_theory_mbh_hdr_row(tg_visual, "cv_mbh", "Coef. Var. (e-1)", panel_labels = c("f", "g", "h")) /
      figure3_theory_mbh_hdr_row(tg_visual, "var_mbh", "Variance (e+4)", panel_labels = c("i", "j", "k")) /
      figure3_theory_mbh_hdr_row(tg_visual, "delta_mbh", "delta (e+1)", panel_labels = c("l", "m", "n"))
  ) +
    patchwork::plot_layout(heights = spec$theory$row_heights) +
    patchwork::plot_annotation(title = "Theoretical")
  empirical_top <- pv_visual_panel_label(figure3_empirical_mdr_row(allouche_model_df), "o") |
    pv_visual_panel_label(figure3_empirical_map_panel(allouche_model_df, root), "p")
  empirical_top <- empirical_top + patchwork::plot_layout(widths = c(1, spec$empirical$map_width))
  empirical <- (
    empirical_top /
      figure3_empirical_hdr_row(empirical_layers, "range", "Range (e+2)", panel_labels = c("q", "r")) /
      figure3_empirical_hdr_row(empirical_layers, "coefficient_of_variation", "Coef. Var. (e-1)", panel_labels = c("s", "t")) /
      figure3_empirical_hdr_row(empirical_layers, "variance", "Variance (e+4)", panel_labels = c("u", "v")) /
      figure3_empirical_hdr_row(empirical_layers, "delta", "delta (e+1)", panel_labels = c("w", "x"))
  ) +
    patchwork::plot_layout(heights = spec$empirical$row_heights) +
    patchwork::plot_annotation(title = "Empirical")
  (theory | empirical) +
    patchwork::plot_layout(widths = c(spec$layout$theory_width, spec$layout$empirical_width)) &
    theme(plot.margin = do.call(margin, as.list(spec$layout$figure_margin)))
}

pv_audit_write_contract_metadata_sidecar <- function(
  data,
  export_contract_id,
  filename,
  metadata_category,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  notes = NA_character_
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }
  row <- pv_audit_contract_row(export_contract_id)
  data <- tibble::as_tibble(data)
  data <- tibble::add_column(
    data,
    export_contract_id = export_contract_id,
    target_id = row$target_id,
    .before = 1
  )
  path <- file.path(audit_paths$root, "metadata", safe_audit_filename(filename))
  ensure_parent_dir(path)
  readr::write_csv(data, path)
  # Compute after the sidecar write so the manifest records the actual file bytes.
  metadata_sidecar_sha256 <- if (file.exists(path)) {
    as.character(tools::sha256sum(path)[[1]])
  } else {
    NA_character_
  }
  write_audit_manifest_row(
    path,
    object_name = paste0(export_contract_id, "_", metadata_category),
    artifact_type = "metadata_sidecar",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = notes
  )
  append_audit_csv_row(
    tibble::tibble(
      audit_run_id = audit_paths$run_id,
      export_contract_id = export_contract_id,
      target_id = row$target_id,
      file = normalizePath(path, winslash = "/", mustWork = FALSE),
      n_rows = nrow(data),
      n_cols = ncol(data),
      status = "metadata_sidecar_written",
      recorded_time = audit_now_utc_string(),
      notes = notes,
      metadata_category = metadata_category,
      metadata_sidecar_file = normalizePath(path, winslash = "/", mustWork = FALSE),
      metadata_sidecar_sha256 = metadata_sidecar_sha256,
      metadata_sidecar_status = "written"
    ),
    audit_paths$audit_contract_export_manifest
  )
  invisible(path)
}

pv_audit_sha256_or_na <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) {
    return(NA_character_)
  }
  as.character(tools::sha256sum(path)[[1]])
}

pv_audit_write_contract_provenance_csv <- function(
  data,
  export_contract_id,
  filename,
  subdir,
  object_name,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  notes = NA_character_
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }
  aligned <- pv_audit_validate_export_schema(data, export_contract_id)
  aligned <- pv_audit_sort_export(aligned, export_contract_id)
  path <- write_audit_csv(
    aligned,
    filename,
    subdir = subdir,
    object_name = object_name,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = notes
  )
  pv_audit_record_contract_export(
    export_contract_id,
    path,
    nrow(aligned),
    ncol(aligned),
    status = "provenance_harmonized_written",
    notes = notes,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates
  )
  invisible(path)
}

# 02 R REIMPLEMENTATION OF SELECTED PELLETT & VALBUENA ANALYSES ####
# ============================================================ #
# 02A distribution helpers reimplemented from selected Julia procedures
# CPGamma(mu, delta) = Gamma(shape = mu / delta, scale = delta)
gamma_shape <- function(mu, delta) mu / delta
gamma_scale <- function(mu, delta) delta

gamma_var_pv <- function(mu, delta) mu * delta
gamma_sd_pv <- function(mu, delta) sqrt(mu * delta)
gamma_cv_pv <- function(mu, delta) sqrt(delta / mu)
gamma_delta_pv <- function(mu, delta) delta
gamma_entropy_pv <- function(mu, delta) {
  a <- gamma_shape(mu, delta)
  theta <- gamma_scale(mu, delta)
  a + log(theta) + lgamma(a) + (1 - a) * digamma(a)
}
gamma_gini_pv <- function(mu, delta) {
  a <- gamma_shape(mu, delta)
  exp(lgamma(a + 0.5) - 0.5 * log(pi) - lgamma(a + 1))
}
gamma_range_pv <- function(mu, delta, probs = c(0.025, 0.975)) {
  qgamma(probs[2], shape = gamma_shape(mu, delta), scale = gamma_scale(mu, delta)) -
    qgamma(probs[1], shape = gamma_shape(mu, delta), scale = gamma_scale(mu, delta))
}

# CPBeta(mu, delta) = Beta(mu * (1/delta - 1), (1 - mu) * (1/delta - 1))
beta_phi <- function(delta) 1 / delta - 1
beta_p <- function(mu, delta) mu * beta_phi(delta)
beta_q <- function(mu, delta) (1 - mu) * beta_phi(delta)
beta_var_pv <- function(mu, delta) mu * (1 - mu) * delta
beta_sd_pv <- function(mu, delta) sqrt(beta_var_pv(mu, delta))
beta_cv_pv <- function(mu, delta) beta_sd_pv(mu, delta) / mu
beta_delta_pv <- function(mu, delta) delta
beta_entropy_pv <- function(mu, delta) {
  p <- beta_p(mu, delta); q <- beta_q(mu, delta)
  lbeta(p, q) - (p - 1) * digamma(p) - (q - 1) * digamma(q) + (p + q - 2) * digamma(p + q)
}
beta_gini_pv <- function(mu, delta) {
  p <- beta_p(mu, delta); q <- beta_q(mu, delta)
  2 * exp(lbeta(2 * p, 2 * q) - log(p) - 2 * lbeta(p, q))
}
beta_range_pv <- function(mu, delta, probs = c(0.025, 0.975)) {
  qbeta(probs[2], shape1 = beta_p(mu, delta), shape2 = beta_q(mu, delta)) -
    qbeta(probs[1], shape1 = beta_p(mu, delta), shape2 = beta_q(mu, delta))
}

# Helper measures for the Pellett & Valbuena concept figure use 0.01 and 0.99 quantiles for range.
mbh_range_gamma <- function(mu, delta) gamma_range_pv(mu, delta, probs = c(0.01, 0.99))
mbh_cv_gamma <- gamma_cv_pv
mbh_var_gamma <- gamma_var_pv
mbh_gini_gamma <- gamma_gini_pv
mbh_delta_gamma <- gamma_delta_pv

# 02B read the archived empirical files containing heterogeneity indices
read_elevation_heterogeneity <- function(root = repo_root) {
  p <- safe_file(root, "02_MBH_comb", "data", "elevation", "heterogeneity_inds.csv")
  readr::read_csv(p, show_col_types = FALSE) %>%
    rename(elev_var = var, iod = iod) %>%
    mutate(across(everything(), as.numeric))
}

read_crop_heterogeneity <- function(root = repo_root) {
  p <- safe_file(root, "02_MBH_comb", "data", "crop_cover", "heterogeneity_inds.csv")
  readr::read_csv(p, show_col_types = FALSE) %>%
    rename(crop_var = var) %>%
    mutate(across(everything(), as.numeric))
}

# 02C partial R ports of the optional procedures for global sampling in Julia
balanced_sample_indices <- function(x, n, bins, seed = 1) {
  set.seed(seed)
  edges <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = bins)
  idx <- seq_along(x)
  n_per_bin <- ceiling(n / bins + 10)
  pieces <- vector("list", length(edges) - 1)
  for (i in seq_len(length(edges) - 1)) {
    in_bin <- idx[x > edges[i] & x < edges[i + 1]]
    if (length(in_bin) == 0) {
      pieces[[i]] <- integer(0)
    } else {
      pieces[[i]] <- sample(in_bin, size = min(n_per_bin, length(in_bin)), replace = FALSE)
    }
  }
  unlist(pieces)
}

calc_elev_window_indices <- function(vals) {
  v <- as.numeric(vals)
  v <- v[is.finite(v) & v > 0]
  if (length(v) < 300) return(tibble(mean = NA_real_, var = NA_real_, sd = NA_real_, cv = NA_real_, iod = NA_real_))
  mu <- mean(v); sig2 <- var(v); sig <- sqrt(sig2)
  tibble(mean = mu, var = sig2, sd = sig, cv = sig / mu, iod = sig2 / mu)
}

calc_crop_window_indices <- function(vals) {
  v <- as.numeric(vals)
  v <- v[is.finite(v) & v < 101]
  if (length(v) < 300) return(tibble(mean = NA_real_, var = NA_real_, sd = NA_real_, cv = NA_real_, delta = NA_real_))
  mu <- mean(v); sig2 <- var(v); sig <- sqrt(sig2)
  tibble(mean = mu, var = sig2, sd = sig, cv = sig / mu, delta = sig2 / (mu * (100 - mu)))
}

# These optional functions are partial R ports of selected Julia procedures for
# global sampling and require the original global rasters. The optional Pellett &
# Valbuena workflow instead reads the archived files containing heterogeneity indices.
resample_global_elevation_indices <- function(root = repo_root, n = GLOBAL_RESAMPLING_N, seed = 1) {
  tif <- safe_file(root, "02_MBH_comb", "data", "elevation", "srtm", "srtm_world.vrt")
  r <- terra::rast(tif)
  set.seed(seed)
  xs <- sample(seq_len(terra::ncol(r)), n, replace = TRUE)
  ys <- sample(seq_len(terra::nrow(r)), n, replace = TRUE)
  centers <- terra::extract(r, data.frame(x = xs, y = ys), cells = FALSE)
  # The extraction above is a partial algorithmic port, not a complete global sample.
  stop("This partial R port does not complete global elevation sampling. The optional workflow uses the archived heterogeneity_inds.csv file instead.")
}

resample_global_crop_indices <- function(root = repo_root, n = GLOBAL_RESAMPLING_N, seed = 1) {
  stop("This partial R port does not complete global crop sampling. The optional workflow uses the archived heterogeneity_inds.csv file instead.")
}

pv_audit_figure2_density_layers <- function() {
  gamma_x <- seq(0.001, 18.999, by = 0.01)
  beta_x <- seq(0.001, 0.999, by = 0.001)
  gamma_layers <- expand_grid(
    distribution = "gamma",
    mu = c(1.5, 4.0, 7.0, 10.0),
    delta_or_alpha = c(0.5, 1.5, 3.5),
    x = gamma_x
  ) %>%
    mutate(
      panel_id = paste0("gamma_mu_", mu, "_delta_", delta_or_alpha),
      density = dgamma(x, shape = gamma_shape(mu, delta_or_alpha), scale = gamma_scale(mu, delta_or_alpha))
    )
  beta_layers <- expand_grid(
    distribution = "beta",
    mu = c(0.1, 0.3, 0.5, 0.7, 0.9),
    delta_or_alpha = c(0.03, 0.15, 0.33),
    x = beta_x
  ) %>%
    mutate(
      panel_id = paste0("beta_mu_", mu, "_delta_", delta_or_alpha),
      density = dbeta(x, shape1 = beta_p(mu, delta_or_alpha), shape2 = beta_q(mu, delta_or_alpha))
    )
  bind_rows(gamma_layers, beta_layers) %>%
    transmute(distribution, panel_id, mu, delta_or_alpha, x, density)
}

pv_audit_figure3_theoretical_fits <- function(models, grid_n = 250) {
  purrr::imap_dfr(models, function(model, measure) {
    terms <- broom::tidy(model) %>%
      transmute(
        measure = measure,
        term,
        estimate,
        std_error = std.error,
        pred_x = NA_real_,
        pred_y = NA_real_,
        grid_n = grid_n
      )
    pred_x <- seq(
      min(model$model[[measure]], na.rm = TRUE),
      max(model$model[[measure]], na.rm = TRUE),
      length.out = 200
    )
    pred_data <- data.frame(pred_x)
    names(pred_data) <- measure
    predictions <- tibble(
      measure = measure,
      term = "fitted_line",
      estimate = NA_real_,
      std_error = NA_real_,
      pred_x = pred_x,
      pred_y = as.numeric(stats::predict(model, newdata = pred_data)),
      grid_n = grid_n
    )
    bind_rows(terms, predictions)
  })
}

pv_audit_julia_polyfit_model <- function(data, predictor_col, response_col = "richness") {
  model_data <- data.frame(
    richness = data[[response_col]],
    mbh = data[[predictor_col]]
  )
  stats::lm(richness ~ mbh + I(mbh^2), data = model_data)
}

pv_audit_julia_polyfit_term_label <- function(term) {
  labels <- c(
    "(Intercept)" = "intercept",
    "mbh" = "mbh",
    "I(mbh^2)" = "mbh_squared"
  )
  out <- unname(labels[term])
  out[is.na(out)] <- term[is.na(out)]
  out
}

pv_audit_figure3_theoretical_fits_julia_style <- function(grid, grid_n = nrow(grid), pred_x_n = 100) {
  measures <- c("range_mbh", "cv_mbh", "var_mbh", "delta_mbh")
  purrr::map_dfr(measures, function(measure) {
    model <- pv_audit_julia_polyfit_model(grid, measure)
    terms <- broom::tidy(model) %>%
      transmute(
        measure = measure,
        term = pv_audit_julia_polyfit_term_label(term),
        estimate,
        std_error = std.error,
        pred_x = NA_real_,
        pred_y = NA_real_,
        grid_n = grid_n
      )
    pred_x <- seq(
      min(grid[[measure]], na.rm = TRUE),
      max(grid[[measure]], na.rm = TRUE),
      length.out = pred_x_n
    )
    predictions <- tibble::tibble(
      measure = measure,
      term = "fitted_line",
      estimate = NA_real_,
      std_error = NA_real_,
      pred_x = pred_x,
      pred_y = as.numeric(stats::predict(model, newdata = data.frame(mbh = pred_x))),
      grid_n = grid_n
    )
    bind_rows(terms, predictions)
  })
}

pv_audit_figure3_empirical_layers <- function(allouche_model_df) {
  raw <- allouche_model_df %>%
    mutate(row_id = dplyr::row_number()) %>%
    transmute(
      row_id,
      richness,
      raw_mu = mu,
      raw_range = range,
      CV,
      sigma_or_sigma2 = sigma2,
      delta,
      mean_elev = mu,
      range = range / 100,
      CV_plot = CV * 10,
      variance = sigma2 / 10000,
      delta_plot = delta / 10
    )
  long <- bind_rows(
    raw %>% transmute(row_id, measure = "mean_elev", x = mean_elev, richness, raw_mu, raw_range, CV, sigma_or_sigma2, delta),
    raw %>% transmute(row_id, measure = "range", x = range, richness, raw_mu, raw_range, CV, sigma_or_sigma2, delta),
    raw %>% transmute(row_id, measure = "CV", x = CV_plot, richness, raw_mu, raw_range, CV, sigma_or_sigma2, delta),
    raw %>% transmute(row_id, measure = "variance", x = variance, richness, raw_mu, raw_range, CV, sigma_or_sigma2, delta),
    raw %>% transmute(row_id, measure = "delta", x = delta_plot, richness, raw_mu, raw_range, CV, sigma_or_sigma2, delta)
  )
  long %>%
    group_by(measure) %>%
    mutate(fit_y = as.numeric(stats::predict(lm(richness ~ x + I(x^2), data = cur_data()), newdata = cur_data()))) %>%
    ungroup()
}

pv_audit_empirical_poisson_formula <- function(degree) {
  if (degree == 1L) {
    return(richness ~ x)
  }
  richness ~ x + I(x^2)
}

pv_audit_figure3_empirical_layers_julia_style <- function(allouche_model_df) {
  sigma_col <- if ("sigma2" %in% names(allouche_model_df)) {
    "sigma2"
  } else if ("sigma" %in% names(allouche_model_df)) {
    "sigma"
  } else {
    stop("Allouche model dataframe lacks sigma2/sigma column for EC-018 harmonized output.", call. = FALSE)
  }
  raw_input <- allouche_model_df
  raw_input$sigma_or_sigma2_exact <- raw_input[[sigma_col]]
  raw <- raw_input %>%
    mutate(row_id = dplyr::row_number()) %>%
    transmute(
      row_id,
      richness,
      raw_mu = mu,
      raw_range = range,
      CV,
      sigma_or_sigma2 = sigma_or_sigma2_exact,
      delta,
      range_scaled = range / 100,
      cv_scaled = CV * 10,
      variance_scaled = sigma_or_sigma2_exact / 10000,
      delta_scaled = delta / 10
    )
  measure_specs <- tibble::tribble(
    ~measure, ~x_column, ~degree,
    "range", "range_scaled", 2L,
    "coefficient_of_variation", "cv_scaled", 2L,
    "variance", "variance_scaled", 2L,
    "delta", "delta_scaled", 1L
  )
  purrr::pmap_dfr(measure_specs, function(measure, x_column, degree) {
    layer <- raw %>%
      transmute(
        row_id,
        measure = measure,
        x = .data[[x_column]],
        richness,
        raw_mu,
        raw_range,
        CV,
        sigma_or_sigma2,
        delta
      )
    model <- stats::glm(
      pv_audit_empirical_poisson_formula(degree),
      data = layer,
      family = stats::poisson(link = "log")
    )
    layer %>%
      mutate(fit_y = as.numeric(stats::predict(model, newdata = layer, type = "response")))
  })
}

pv_audit_review2_summary <- function(data, simulation, seed_status, distribution = NULL, notes = NA_character_) {
  quantiles <- paste(
    paste0(
      names(stats::quantile(data$mu, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),
      "_mu=",
      signif(as.numeric(stats::quantile(data$mu, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)), 8)
    ),
    collapse = ";"
  )
  base <- tibble(
    simulation = simulation,
    seed_status = seed_status,
    n = nrow(data),
    mean_mu = mean(data$mu, na.rm = TRUE),
    mean_delta = mean(data$delta, na.rm = TRUE),
    mean_var = mean(data$var, na.rm = TRUE),
    quantiles = quantiles,
    notes = notes
  )
  if (!is.null(distribution)) {
    base <- mutate(base, distribution = distribution, .after = n)
  }
  base
}

pv_audit_hex_loess_layers <- function(data, simulation, panels, bins = 120, span = 0.75) {
  purrr::map_dfr(panels, function(panel) {
    x <- data$mu
    y <- data[[panel]]
    hb <- hexbin::hexbin(x = x, y = y, xbins = bins)
    centers <- hexbin::hcell2xy(hb)
    hex_rows <- tibble(
      simulation = simulation,
      panel = panel,
      bin_x = centers$x,
      bin_y = centers$y,
      count = as.numeric(hb@count),
      loess_x = NA_real_,
      loess_y = NA_real_,
      bins = bins,
      span = span
    )
    loess_x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)
    model <- stats::loess(y ~ x, span = span)
    loess_rows <- tibble(
      simulation = simulation,
      panel = panel,
      bin_x = NA_real_,
      bin_y = NA_real_,
      count = NA_real_,
      loess_x = loess_x,
      loess_y = as.numeric(stats::predict(model, newdata = data.frame(x = loess_x))),
      bins = bins,
      span = span
    )
    bind_rows(hex_rows, loess_rows)
  })
}

pv_audit_ec034_phase_marker <- function(
  audit_paths,
  phase,
  status = "marker",
  simulation = NA_character_,
  panel = NA_character_,
  input_n = NA_integer_,
  loess_x_n = NA_integer_,
  bins = NA_integer_,
  span = NA_real_,
  elapsed_seconds = NA_real_,
  details = NA_character_
) {
  tryCatch({
    if (is.null(audit_paths) || is.null(audit_paths$root) || !nzchar(audit_paths$root)) {
      return(invisible(FALSE))
    }
    if (!dir.exists(audit_paths$root)) {
      return(invisible(FALSE))
    }
    marker_path <- file.path(audit_paths$root, "metadata", "ec034_phase_markers.csv")
    recorded_time <- if (exists("audit_now_utc_string", mode = "function")) {
      audit_now_utc_string()
    } else {
      format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    }
    append_audit_csv_row(
      tibble::tibble(
        recorded_time = recorded_time,
        phase = phase,
        status = status,
        simulation = simulation,
        panel = panel,
        input_n = input_n,
        loess_x_n = loess_x_n,
        bins = bins,
        span = span,
        elapsed_seconds = elapsed_seconds,
        details = details
      ),
      marker_path,
      schema = c(
        "recorded_time", "phase", "status", "simulation", "panel",
        "input_n", "loess_x_n", "bins", "span", "elapsed_seconds", "details"
      )
    )
    invisible(TRUE)
  }, error = function(e) invisible(FALSE))
}

pv_audit_julia_style_loess_layers <- function(
  data,
  simulation,
  panel_map,
  bins = 300,
  span = 0.4,
  loess_x = seq(0.001, 0.999, by = 0.001),
  audit_paths = NULL,
  phase_prefix = NULL,
  loess_control = NULL
) {
  if (is.null(loess_control)) {
    loess_control <- stats::loess.control(
      surface = "interpolate",
      statistics = "none",
      trace.hat = "approximate"
    )
  }
  purrr::imap_dfr(panel_map, function(panel_column, panel_label) {
    phase_simulation <- if (!is.null(phase_prefix) && length(phase_prefix) == 1L && !is.na(phase_prefix) && nzchar(phase_prefix)) {
      phase_prefix
    } else {
      simulation
    }
    phase_panel <- paste0(phase_simulation, "_", panel_label)
    input_n <- nrow(data)
    loess_x_n <- length(loess_x)
    pv_audit_ec034_phase_marker(
      audit_paths,
      paste0("before_", phase_panel, "_loess"),
      simulation = simulation,
      panel = panel_label,
      input_n = input_n,
      loess_x_n = loess_x_n,
      bins = bins,
      span = span,
      details = paste0("panel_column=", panel_column, ";before data extraction and loess fit")
    )
    x <- data$mu
    y <- data[[panel_column]]
    loess_start <- proc.time()[["elapsed"]]
    model <- stats::loess(y ~ x, span = span, control = loess_control)
    pv_audit_ec034_phase_marker(
      audit_paths,
      paste0("after_", phase_panel, "_loess"),
      simulation = simulation,
      panel = panel_label,
      input_n = input_n,
      loess_x_n = loess_x_n,
      bins = bins,
      span = span,
      elapsed_seconds = proc.time()[["elapsed"]] - loess_start,
      details = paste0("panel_column=", panel_column, ";loess fit complete")
    )
    pv_audit_ec034_phase_marker(
      audit_paths,
      paste0("before_", phase_panel, "_predict"),
      simulation = simulation,
      panel = panel_label,
      input_n = input_n,
      loess_x_n = loess_x_n,
      bins = bins,
      span = span,
      details = paste0("panel_column=", panel_column, ";before loess prediction")
    )
    predict_start <- proc.time()[["elapsed"]]
    loess_y <- as.numeric(stats::predict(model, newdata = data.frame(x = loess_x)))
    pv_audit_ec034_phase_marker(
      audit_paths,
      paste0("after_", phase_panel, "_predict"),
      simulation = simulation,
      panel = panel_label,
      input_n = input_n,
      loess_x_n = loess_x_n,
      bins = bins,
      span = span,
      elapsed_seconds = proc.time()[["elapsed"]] - predict_start,
      details = paste0("panel_column=", panel_column, ";loess prediction complete")
    )
    layer <- tibble(
      simulation = simulation,
      panel = panel_label,
      bin_x = NA_real_,
      bin_y = NA_real_,
      count = NA_real_,
      loess_x = loess_x,
      loess_y = loess_y,
      bins = bins,
      span = span
    )
    pv_audit_ec034_phase_marker(
      audit_paths,
      paste0("after_", phase_panel, "_layer_assembly"),
      simulation = simulation,
      panel = panel_label,
      input_n = input_n,
      loess_x_n = loess_x_n,
      bins = bins,
      span = span,
      details = paste0("panel_column=", panel_column, ";layer rows=", nrow(layer))
    )
    layer
  })
}

pv_audit_review2_ec034_julia_layers <- function(
  pqdf,
  mudf,
  bins = 300,
  span = 0.4,
  loess_x = seq(0.001, 0.999, by = 0.001),
  audit_paths = NULL,
  loess_control = NULL
) {
  bind_rows(
    pv_audit_julia_style_loess_layers(
      pqdf,
      simulation = "pqsolution",
      panel_map = c(variance = "var", delta = "delta"),
      bins = bins,
      span = span,
      loess_x = loess_x,
      audit_paths = audit_paths,
      phase_prefix = "pqsolution",
      loess_control = loess_control
    ),
    pv_audit_julia_style_loess_layers(
      mudf,
      simulation = "muvarsolution",
      panel_map = c(variance = "var", delta = "delta"),
      bins = bins,
      span = span,
      loess_x = loess_x,
      audit_paths = audit_paths,
      phase_prefix = "muvarsolution",
      loess_control = loess_control
    )
  )
}

pv_audit_write_review2_ec034_julia_layers <- function(
  pqdf,
  mudf,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  bins = 300,
  span = 0.4,
  loess_x = seq(0.001, 0.999, by = 0.001),
  loess_control = NULL
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  pv_audit_ec034_phase_marker(
    audit_paths,
    "before_ec034_layer_construction",
    input_n = nrow(pqdf) + nrow(mudf),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    details = "before pv_audit_review2_ec034_julia_layers"
  )
  layer_start <- proc.time()[["elapsed"]]
  ec034_layers <- pv_audit_review2_ec034_julia_layers(
    pqdf,
    mudf,
    bins = bins,
    span = span,
    loess_x = loess_x,
    audit_paths = audit_paths,
    loess_control = loess_control
  )
  pv_audit_ec034_phase_marker(
    audit_paths,
    "after_ec034_layer_construction",
    input_n = nrow(pqdf) + nrow(mudf),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    elapsed_seconds = proc.time()[["elapsed"]] - layer_start,
    details = paste0("ec034_layer_rows=", nrow(ec034_layers))
  )
  pv_audit_ec034_phase_marker(
    audit_paths,
    "before_ec034_contract_csv_write",
    input_n = nrow(ec034_layers),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    details = "before pv_audit_write_contract_csv"
  )
  contract_path <- pv_audit_write_contract_csv(
    ec034_layers,
    "EC-034",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = "Review2 Julia-style loess-only layer data"
  )
  pv_audit_ec034_phase_marker(
    audit_paths,
    "after_ec034_contract_csv_write",
    input_n = nrow(ec034_layers),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    details = paste0("contract_path=", contract_path)
  )
  # Metadata records the harmonization policy; layer values are computed from R data.
  pv_audit_ec034_phase_marker(
    audit_paths,
    "before_ec034_metadata_sidecar_write",
    input_n = nrow(ec034_layers),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    details = "before pv_audit_write_contract_metadata_sidecar"
  )
  sidecar_path <- pv_audit_write_contract_metadata_sidecar(
    tibble(
      method_policy = "julia_style_review2_loess_layers",
      output_schema = "simulation,panel,bin_x,bin_y,count,loess_x,loess_y,bins,span",
      key_fields = "simulation,panel,loess_x",
      loess_only_policy = TRUE,
      intentionally_empty_fields = "bin_x,bin_y,count",
      bins = bins,
      span = span,
      loess_x_min = min(loess_x, na.rm = TRUE),
      loess_x_max = max(loess_x, na.rm = TRUE),
      loess_x_step = unique(round(diff(loess_x), 12))[1],
      loess_x_n = length(loess_x),
      expected_total_rows = 4L * length(loess_x),
      expected_rows_per_simulation_panel = length(loess_x),
      simulation_label_mapping = "pqdf->pqsolution;mudf->muvarsolution",
      panel_label_mapping = "var->variance;delta->delta",
      value_generation_policy = "computed_from_R_generated_review2_data_not_julia_csv"
    ),
    "EC-034",
    "stage8lnq_ec034_julia_style_loess_metadata_r.csv",
    "julia_style_loess_metadata",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = "Stage 8L-NQ EC-034 Julia-style loess metadata sidecar"
  )
  pv_audit_ec034_phase_marker(
    audit_paths,
    "after_ec034_metadata_sidecar_write",
    input_n = nrow(ec034_layers),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    details = paste0("sidecar_path=", sidecar_path)
  )
  pv_audit_ec034_phase_marker(
    audit_paths,
    "after_ec034_final_artifact_checks",
    input_n = nrow(ec034_layers),
    loess_x_n = length(loess_x),
    bins = bins,
    span = span,
    details = paste0(
      "contract_exists=", file.exists(contract_path),
      ";sidecar_exists=", file.exists(sidecar_path)
    )
  )
  invisible(ec034_layers)
}

# 02D R reimplementation of selected Figure 1 procedures from concept_col1.jl
run_figure1_concept_R <- function(
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  x <- seq(0.01, 13, by = 0.01)
  pdf_df <- expand_grid(mu = c(0.5, 1.5, 3.0), x = x) %>%
    mutate(density = dgamma(x, shape = gamma_shape(mu, 1), scale = 1), mu_lab = paste0("mu = ", mu))

  p_a <- ggplot(pdf_df, aes(x, density)) +
    geom_area(fill = "grey85") + geom_line() + facet_wrap(~mu_lab, nrow = 1) +
    labs(x = "x", y = "Density", title = "a. Lower-bounded distributions") + theme_bw()

  mur <- seq(1, 41, by = 0.01); H <- 1.5
  mbh_df <- tibble(
    mu = mur,
    `Coef. var.` = rescale01(mbh_cv_gamma(mur, H)),
    Range = rescale01(mbh_range_gamma(mur, H)),
    Variance = rescale01(mbh_var_gamma(mur, H)) / 1.3,
    Gini = rescale01(mbh_gini_gamma(mur, H + 5))
  ) %>% pivot_longer(-mu, names_to = "measure", values_to = "heterogeneity")
  p_b <- ggplot(mbh_df, aes(mu, heterogeneity, linetype = measure)) +
    geom_line() + labs(x = "Mean", y = "Heterogeneity", title = "b. Mean-biased measures") + theme_bw()

  gen_mdr <- function(mu, intercept = 50) mu * (intercept / 1.6 - mu) * 0.1 + intercept
  gen_hdr <- function(h, intercept = 50) h
  p_c_df <- bind_rows(
    tibble(x = mur, y = rescale01(gen_mdr(mur)), relationship = "True MDR"),
    tibble(x = mur, y = rescale01(gen_hdr(mur)), relationship = "True HDR")
  )
  p_c <- p_c_df %>% ggplot(aes(x, y)) + geom_line() + facet_wrap(~relationship, scales = "free_x") +
    labs(x = "Mean or heterogeneity", y = "Diversity", title = "c. True relationships") + theme_bw()

  coefvar <- seq(0, 5.4, by = 0.1); rangev <- seq(0, 15.9, by = 0.5)
  variance <- seq(0, 15.4, by = 0.5); gini <- seq(0, 1, by = 0.01)
  p_d_df <- bind_rows(
    tibble(x = rescale01(coefvar), y = 47 + 57.8 * coefvar - 9.65 * coefvar^2, measure = "Coef. var."),
    tibble(x = rescale01(rangev), y = 82 + 5 * rangev - 0.32 * rangev^2, measure = "Range"),
    tibble(x = rescale01(variance), y = 96 + 2.2 * variance - 0.3 * variance^2, measure = "Variance"),
    tibble(x = rescale01(gini), y = 23 + 196.6 * gini - 93.65 * gini^2, measure = "Gini")
  )
  p_d <- ggplot(p_d_df, aes(x, y, linetype = measure)) + geom_line() +
    labs(x = "Observed heterogeneity", y = "Diversity", title = "d. Observed mean-biased HDRs") + theme_bw()

  pv_audit_write_contract_csv(
    pdf_df %>% transmute(panel_id = paste0("mu_", mu), mu, delta = 1, x, density, layer_type = "gamma_density"),
    "EC-001",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 1 panel a gamma densities"
  )
  pv_audit_write_contract_csv(
    mbh_df %>% mutate(H = H, scale_note = "rescale01; variance divided by 1.3; gini uses H+5") %>%
      transmute(measure, mu, heterogeneity, H, scale_note),
    "EC-002",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 1 MBH lines"
  )
  pv_audit_write_contract_csv(
    bind_rows(
      p_c_df %>%
        mutate(
          panel_id = "panel_c",
          measure = relationship,
          formula_label = ifelse(relationship == "True MDR", "gen_mdr(mu)", "gen_hdr(H)")
        ) %>%
        transmute(panel_id, relationship, measure, x, y, formula_label),
      p_d_df %>%
        mutate(
          panel_id = "panel_d",
          relationship = "Observed mean-biased HDR",
          formula_label = paste0("observed_hdr_", measure)
        ) %>%
        transmute(panel_id, relationship, measure, x, y, formula_label)
    ),
    "EC-003",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 1 MDR/HDR and observed HDR layers"
  )

  diagnostic_fig <- (p_a / (p_b | p_c | p_d)) + plot_layout(heights = c(1, 1.2))
  fig <- assemble_figure1_visual_replication(pdf_df, mbh_df, p_c_df, p_d_df)
  save_png(diagnostic_fig, "Figure1_R_diagnostic_summary.png", width = 11, height = 7)
  save_png(fig, "Figure1_R_visual_replication.png", width = 6.4, height = 4.6)
  save_png(fig, "Figure1_R_equivalent.png", width = 6.4, height = 4.6)
  invisible(fig)
}

# 02E R reimplementation of selected Figure 2 procedures from MBH_comb_figure.jl
make_gamma_measure_data <- function() {
  mu <- seq(0.01, 15, by = 0.01)
  deltas <- c(0.5, 1.5, 3.5)
  expand_grid(mu = mu, delta = deltas) %>%
    mutate(
      Variance = gamma_var_pv(mu, delta),
      `Std. dev.` = gamma_sd_pv(mu, delta),
      `Coef. var.` = gamma_cv_pv(mu, delta),
      `delta_L` = delta,
      Entropy = gamma_entropy_pv(mu, delta),
      Gini = gamma_gini_pv(mu, delta),
      Range = gamma_range_pv(mu, delta)
    ) %>% pivot_longer(c(Variance, `Std. dev.`, `Coef. var.`, delta_L, Entropy, Gini, Range),
                       names_to = "measure", values_to = "value")
}

make_beta_measure_data <- function() {
  mu <- seq(0.001, 0.999, by = 0.001)
  deltas <- c(0.03, 0.15, 0.33)
  expand_grid(mu = mu, delta = deltas) %>%
    mutate(
      Variance = beta_var_pv(mu, delta),
      `Std. dev.` = beta_sd_pv(mu, delta),
      `Coef. var.` = beta_cv_pv(mu, delta),
      `delta_2` = delta,
      Entropy = beta_entropy_pv(mu, delta),
      Gini = beta_gini_pv(mu, delta),
      Range = beta_range_pv(mu, delta)
    ) %>% pivot_longer(c(Variance, `Std. dev.`, `Coef. var.`, delta_2, Entropy, Gini, Range),
                       names_to = "measure", values_to = "value")
}

run_figure2_MBH_comb_R <- function(
  root = repo_root,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  gamma_df <- make_gamma_measure_data()
  beta_df <- make_beta_measure_data()
  readr::write_csv(gamma_df, safe_file(intermediate_dir, "figure2_gamma_measure_curves_R.csv"))
  readr::write_csv(beta_df, safe_file(intermediate_dir, "figure2_beta_measure_curves_R.csv"))
  pv_audit_write_contract_csv(
    gamma_df %>% mutate(distribution = "gamma", .before = 1) %>% transmute(distribution, measure, mu, delta, value),
    "EC-004",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 2 gamma measure curves"
  )
  pv_audit_write_contract_csv(
    beta_df %>% mutate(distribution = "beta", .before = 1) %>% transmute(distribution, measure, mu, delta, value),
    "EC-005",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 2 beta measure curves"
  )
  pv_audit_write_contract_csv(
    pv_audit_figure2_density_layers(),
    "EC-006",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 2 gamma and beta density panels"
  )

  p_gamma <- ggplot(gamma_df, aes(mu, value, group = delta, linetype = factor(delta))) +
    geom_line() + facet_wrap(~measure, scales = "free_y", ncol = 3) +
    labs(title = "Lower-bounded gamma variables", x = "Mean", y = "Measure", linetype = "delta_L") + theme_bw()
  p_beta <- ggplot(beta_df, aes(mu, value, group = delta, linetype = factor(delta))) +
    geom_line() + facet_wrap(~measure, scales = "free_y", ncol = 3) +
    labs(title = "Double-bounded beta variables", x = "Mean", y = "Measure", linetype = "delta_2") + theme_bw()

  elev <- read_elevation_heterogeneity(root)
  crop <- read_crop_heterogeneity(root)
  elev_long <- elev %>%
    transmute(mean, Variance = elev_var / 10000, `Std. dev.` = sd / 100, `Coef. var.` = cv, delta_L = iod) %>%
    pivot_longer(-mean, names_to = "measure", values_to = "value")
  crop_long <- crop %>%
    transmute(mean, Variance = crop_var / 1000, `Std. dev.` = sd, `Coef. var.` = cv, delta_2 = delta) %>%
    pivot_longer(-mean, names_to = "measure", values_to = "value")
  pv_audit_write_contract_csv(
    elev %>%
      mutate(source_row_id = dplyr::row_number()) %>%
      {
        bind_rows(
          transmute(., source_row_id, mean, measure = "Variance", value = elev_var / 10000, raw_var = elev_var, sd, cv, iod, scale_note = "elev_var/10000"),
          transmute(., source_row_id, mean, measure = "Std. dev.", value = sd / 100, raw_var = elev_var, sd, cv, iod, scale_note = "sd/100"),
          transmute(., source_row_id, mean, measure = "Coef. var.", value = cv, raw_var = elev_var, sd, cv, iod, scale_note = "cv raw"),
          transmute(., source_row_id, mean, measure = "delta_L", value = iod, raw_var = elev_var, sd, cv, iod, scale_note = "iod raw")
        )
      },
    "EC-007",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 2 elevation empirical layer"
  )
  pv_audit_write_contract_csv(
    crop %>%
      mutate(source_row_id = dplyr::row_number()) %>%
      {
        bind_rows(
          transmute(., source_row_id, mean, measure = "Variance", value = crop_var / 1000, raw_var = crop_var, sd, cv, delta, scale_note = "crop_var/1000"),
          transmute(., source_row_id, mean, measure = "Std. dev.", value = sd, raw_var = crop_var, sd, cv, delta, scale_note = "sd raw"),
          transmute(., source_row_id, mean, measure = "Coef. var.", value = cv, raw_var = crop_var, sd, cv, delta, scale_note = "cv raw"),
          transmute(., source_row_id, mean, measure = "delta_2", value = delta, raw_var = crop_var, sd, cv, delta, scale_note = "delta raw")
        )
      },
    "EC-008",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 2 crop empirical layer"
  )

  p_elev <- ggplot(elev_long, aes(mean, value)) + geom_hex(bins = 90) +
    scale_fill_viridis_c(trans = "log10") + facet_wrap(~measure, scales = "free_y", ncol = 2) +
    labs(title = "Empirical SRTM land elevation", x = "Mean elevation", y = "Measure") + theme_bw()
  p_crop <- ggplot(crop_long, aes(mean, value)) + geom_hex(bins = 80) +
    scale_fill_viridis_c(trans = "log10") + facet_wrap(~measure, scales = "free_y", ncol = 2) +
    labs(title = "Empirical crop cover", x = "Mean crop cover (%)", y = "Measure") + theme_bw()

  p_elev_context <- figure2_srtm_context_panel(root)
  p_crop_context <- figure2_crop_context_panel(root)
  diagnostic_fig <- (p_gamma | p_beta) / (p_elev | p_crop)
	  fig_visual <- assemble_figure2_visual_replication(
	    gamma_df,
	    beta_df,
	    pv_audit_figure2_density_layers(),
	    elev_long,
	    crop_long,
	    elev,
	    crop,
	    root = root
	  )
  save_png(p_gamma, "Figure2_gamma_theory_R_equivalent.png", width = 11, height = 7)
  save_png(p_beta, "Figure2_beta_theory_R_equivalent.png", width = 11, height = 7)
  save_png(p_elev, "Figure2_elevation_empirical_R_equivalent.png", width = 8, height = 7)
  save_png(p_crop, "Figure2_crop_empirical_R_equivalent.png", width = 8, height = 7)
  save_png(p_elev_context, "Figure2_elevation_context_map_R_equivalent.png", width = 5.2, height = 2.4)
  save_png(p_crop_context, "Figure2_crop_context_map_R_equivalent.png", width = 5.2, height = 2.4)
  save_png(diagnostic_fig, "Figure2_R_diagnostic_summary.png", width = 16, height = 14)
  save_png(fig_visual, "Figure2_R_visual_replication.png", width = 11.3, height = 10.0)
  save_png(fig_visual, "Figure2_R_equivalent.png", width = 11.3, height = 10.0)
  invisible(list(
    gamma = gamma_df,
    beta = beta_df,
    elevation = elev,
    crop = crop,
    elevation_context_map = p_elev_context,
    crop_context_map = p_crop_context,
    visual_replication = fig_visual,
    diagnostic_summary = diagnostic_fig
  ))
}

# 02F 02_MBH_comb/empirical_hypo_test.jl
paired_residual_test_row <- function(name, h0, h1) {
  tt <- t.test(h0, h1, paired = TRUE)
  tibble(
    test = name,
    mean_h0_res2 = mean(h0, na.rm = TRUE),
    mean_h1_res2 = mean(h1, na.rm = TRUE),
    mean_difference_h0_minus_h1 = mean(h0 - h1, na.rm = TRUE),
    t_statistic = as.numeric(tt$statistic),
    df = as.numeric(tt$parameter),
    p_value = as.numeric(tt$p.value),
    conf_low = tt$conf.int[1],
    conf_high = tt$conf.int[2]
  )
}

run_mean_biased_hypothesis_tests_R <- function(
  root = repo_root,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  elev <- read_elevation_heterogeneity(root)
  crop <- read_crop_heterogeneity(root)
  out <- list()

  elev_h0_var <- (elev$elev_var - mean(elev$elev_var))^2
  elev_h1_var <- (elev$elev_var - (mean(elev$iod) * elev$mean))^2
  out[["elevation_variance"]] <- paired_residual_test_row("elevation_variance", elev_h0_var, elev_h1_var)

  elev_h0_sd <- (elev$sd - mean(elev$sd))^2
  elev_h1_sd <- (elev$sd - sqrt(mean(elev$iod) * elev$mean))^2
  out[["elevation_sd"]] <- paired_residual_test_row("elevation_sd", elev_h0_sd, elev_h1_sd)

  elev_h0_cv <- (elev$cv - mean(elev$cv))^2
  elev_h1_cv <- (elev$cv - sqrt(median(elev$iod) * elev$mean) / elev$mean)^2
  out[["elevation_cv"]] <- paired_residual_test_row("elevation_cv", elev_h0_cv, elev_h1_cv)

  crop_h0_var <- (crop$crop_var - mean(crop$crop_var))^2
  crop_h1_var <- (crop$crop_var - (mean(crop$delta) * (crop$mean * (100 - crop$mean))))^2
  out[["crop_variance"]] <- paired_residual_test_row("crop_variance", crop_h0_var, crop_h1_var)

  crop_h0_sd <- (crop$sd - mean(crop$sd))^2
  crop_h1_sd <- (crop$sd - sqrt(mean(crop$delta) * (crop$mean * (100 - crop$mean))))^2
  out[["crop_sd"]] <- paired_residual_test_row("crop_sd", crop_h0_sd, crop_h1_sd)

  crop_h0_cv <- (crop$cv - mean(crop$cv))^2
  crop_h1_cv <- (crop$cv - sqrt(median(crop$delta) * (crop$mean * (100 - crop$mean))) / crop$mean)^2
  out[["crop_cv"]] <- paired_residual_test_row("crop_cv", crop_h0_cv, crop_h1_cv)

  ans <- bind_rows(out)
  residual_vectors <- bind_rows(
    tibble(test = "elevation_variance", row_id = seq_along(elev_h0_var), h0_res2 = elev_h0_var, h1_res2 = elev_h1_var, input_mean = elev$mean, input_measure = elev$elev_var, source = "elevation"),
    tibble(test = "elevation_sd", row_id = seq_along(elev_h0_sd), h0_res2 = elev_h0_sd, h1_res2 = elev_h1_sd, input_mean = elev$mean, input_measure = elev$sd, source = "elevation"),
    tibble(test = "elevation_cv", row_id = seq_along(elev_h0_cv), h0_res2 = elev_h0_cv, h1_res2 = elev_h1_cv, input_mean = elev$mean, input_measure = elev$cv, source = "elevation"),
    tibble(test = "crop_variance", row_id = seq_along(crop_h0_var), h0_res2 = crop_h0_var, h1_res2 = crop_h1_var, input_mean = crop$mean, input_measure = crop$crop_var, source = "crop_cover"),
    tibble(test = "crop_sd", row_id = seq_along(crop_h0_sd), h0_res2 = crop_h0_sd, h1_res2 = crop_h1_sd, input_mean = crop$mean, input_measure = crop$sd, source = "crop_cover"),
    tibble(test = "crop_cv", row_id = seq_along(crop_h0_cv), h0_res2 = crop_h0_cv, h1_res2 = crop_h1_cv, input_mean = crop$mean, input_measure = crop$cv, source = "crop_cover")
  )
  pv_audit_write_contract_csv(
    residual_vectors,
    "EC-009",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Mean-biased residual-squared vectors"
  )
  pv_audit_write_contract_csv(
    ans %>% transmute(test, t_statistic, df, p_value, mean_h0_res2, mean_h1_res2, method = "paired_t_test"),
    "EC-010",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Mean-biased test table"
  )
  readr::write_csv(ans, safe_file(tab_dir, "mean_biased_hypothesis_tests_R.csv"))
  ans
}

# 02G 03_corrected_HDR/empirical/prepare_data.jl + analysis.jl

theory_grid_R <- function(n = 500) {
  intercept <- 66.3947 / 2
  mu_range <- seq(64, 2422, length.out = n)
  H_range <- seq(0.64, 144, length.out = n)
  grid <- expand_grid(mu = mu_range, H = H_range)
  MDR <- function(mu) intercept + mu * 0.0303389 + mu^2 * -1.18841e-5
  HDR <- function(H) H * 0.107633
  grid %>% mutate(richness = intercept + MDR(mu) + HDR(H),
                  range_mbh = mbh_range_gamma(mu, H),
                  cv_mbh = mbh_cv_gamma(mu, H),
                  var_mbh = mbh_var_gamma(mu, H),
                  delta_mbh = mbh_delta_gamma(mu, H))
}

run_figure3_corrected_HDR_R <- function(
  allouche_results,
  root = repo_root,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  # Theoretical summaries of model fits for mean-biased and mean-independent heterogeneity.
  # Use a moderate grid for speed; the equations are the same as the Julia script.
  tg <- theory_grid_R(n = 250)
  theory_models <- list(
    range_mbh = lm(richness ~ range_mbh + I(range_mbh^2), data = tg),
    cv_mbh = lm(richness ~ cv_mbh + I(cv_mbh^2), data = tg),
    var_mbh = lm(richness ~ var_mbh + I(var_mbh^2), data = tg),
    delta_mbh = lm(richness ~ delta_mbh + I(delta_mbh^2), data = tg)
  )
  theory_terms <- write_model_terms(theory_models, "figure3_theoretical_polyfits_R.csv")
  pv_audit_write_contract_csv(
    tg %>% mutate(grid_n = 250) %>% transmute(mu, H, richness, range_mbh, cv_mbh, var_mbh, delta_mbh, grid_n),
    "EC-016",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 3 theoretical full grid"
  )
  # Stage 8R-NQ instrumentation: EC-016 metadata sidecar for units and scales only.
  pv_audit_write_contract_metadata_sidecar(
    tibble::tribble(
      ~column_name, ~column_role, ~unit_or_scale, ~scale_basis, ~formula_source, ~grid_n, ~metadata_note,
      "mu", "input_axis", "original mu grid units", "seq(64, 2422, length.out = grid_n)", "mu_range", 250, "No primary EC-016 values are changed.",
      "H", "input_axis", "original H grid units", "seq(0.64, 144, length.out = grid_n)", "H_range", 250, "No primary EC-016 values are changed.",
      "richness", "response", "predicted richness", "intercept + MDR(mu) + HDR(H)", "richness expression in theory_grid_R", 250, "No primary EC-016 values are changed.",
      "range_mbh", "computed_measure", "mbh_range_gamma output", "0.01/0.99 gamma quantile range helper output", "mbh_range_gamma(mu, H)", 250, "No primary EC-016 values are changed.",
      "cv_mbh", "computed_measure", "mbh_cv_gamma output", "gamma coefficient-of-variation helper output", "mbh_cv_gamma(mu, H)", 250, "No primary EC-016 values are changed.",
      "var_mbh", "computed_measure", "mbh_var_gamma output", "gamma variance helper output", "mbh_var_gamma(mu, H)", 250, "No primary EC-016 values are changed.",
      "delta_mbh", "computed_measure", "mbh_delta_gamma output", "gamma delta helper output", "mbh_delta_gamma(mu, H)", 250, "No primary EC-016 values are changed.",
      "grid_n", "metadata", "row-generation grid size", "n passed to theory_grid_R", "grid_n = 250", 250, "No primary EC-016 values are changed."
    ),
    "EC-016",
    "stage8rnq_patch_EC016_unit_metadata_r.csv",
    "unit_annotation",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-016 unit/scale metadata sidecar"
  )
  pv_audit_write_contract_csv(
    pv_audit_figure3_theoretical_fits(theory_models, grid_n = 250),
    "EC-017",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 3 theoretical model terms and fitted lines"
  )
  # Stage 8R-NQ instrumentation: EC-017 metadata sidecar for model statistics only.
  pv_audit_write_contract_metadata_sidecar(
    purrr::imap_dfr(theory_models, function(model, measure) {
      pred_range <- range(model$model[[measure]], na.rm = TRUE)
      tibble::tibble(
        measure = measure,
        model_formula = paste(deparse(stats::formula(model)), collapse = " "),
        fit_engine = "stats::lm",
        predictor_column = measure,
        term_convention = "linear term plus I(measure^2)",
        statistic_source = "broom::tidy for model terms; fitted_line rows have NA statistic fields",
        std_error_source = "broom::tidy std.error for model terms",
        p_value_source = "not included in EC-017 primary output",
        prediction_source = "stats::predict over seq(min predictor, max predictor, length.out = 200)",
        pred_x_min = pred_range[1],
        pred_x_max = pred_range[2],
        pred_x_n = 200,
        grid_n = 250
      )
    }),
    "EC-017",
    "stage8rnq_patch_EC017_model_statistic_metadata_r.csv",
    "model_statistic_metadata",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-017 model-statistic metadata sidecar"
  )

  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    # This harmonized export, used only for provenance, mirrors the grid density
    # and fitted line convention of the P&V Julia audit for EC-017.
    # The comparison export uses Julia's MBH predictor scaling; the main
    # theory_grid_R helper and manuscript outputs remain on R scales.
    tg_exact <- theory_grid_R(n = 500) %>%
      mutate(
        range_mbh = range_mbh / 100,
        cv_mbh = cv_mbh * 10,
        var_mbh = var_mbh / 10000,
        delta_mbh = delta_mbh / 10
      )
    ec017_exact <- pv_audit_figure3_theoretical_fits_julia_style(
      tg_exact,
      grid_n = nrow(tg_exact),
      pred_x_n = 100
    )
    ec017_exact_path <- pv_audit_write_contract_provenance_csv(
      ec017_exact,
      "EC-017",
      "stage8lnq_exact_EC017_theoretical_fits_grid250000_r.csv",
      subdir = "tables",
      object_name = "EC-017_harmonized_theoretical_fits_grid250000",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Stage 8L-NQ exact-replication harmonized EC-017 output using Julia-style grid_n=250000 and fitted grid n=100"
    )
    pv_audit_write_contract_metadata_sidecar(
      tibble::tibble(
        output_name = "stage8lnq_exact_EC017_theoretical_fits_grid250000_r.csv",
        output_sha256 = pv_audit_sha256_or_na(ec017_exact_path),
        output_row_count = nrow(ec017_exact),
        output_column_count = ncol(ec017_exact),
        model_family = "Gaussian linear model via stats::lm",
        model_basis = "intercept + mbh + mbh^2",
        coefficient_terms = "intercept;mbh;mbh_squared",
        grid_axis_n = 500L,
        grid_n = nrow(tg_exact),
        pred_x_n = 100L,
        mbh_scale_convention = "range_mbh=mbh_range_gamma/100; cv_mbh=mbh_cv_gamma*10; var_mbh=mbh_var_gamma/10000; delta_mbh=mbh_delta_gamma/10",
        key_fields = "coefficients: measure,term; fitted values: measure,pred_x",
        source_function = "pv_audit_figure3_theoretical_fits_julia_style",
        provenance_only = TRUE,
        manuscript_facing_output_changed = FALSE,
        metadata_note = "grid_n is total rows, matching Julia length(mu)=500*500; MBH predictors are scaled to Julia exact-replication units; no Julia output values are hard-coded."
      ),
      "EC-017",
      "stage8lnq_exact_EC017_theoretical_fits_grid250000_metadata_r.csv",
      "exact_replication_harmonization_metadata",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Stage 8L-NQ exact-replication EC-017 harmonized output metadata"
    )
  }

  theory_long <- tg %>%
    sample_n(min(25000, nrow(tg))) %>%
    pivot_longer(c(range_mbh, cv_mbh, var_mbh, delta_mbh), names_to = "measure", values_to = "heterogeneity")
  p_theory <- ggplot(theory_long, aes(heterogeneity, richness)) +
    geom_point(alpha = 0.06, size = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
    facet_wrap(~measure, scales = "free_x", ncol = 2) +
    labs(title = "Theoretical mean-biased and mean-independent HDRs", x = "Heterogeneity measure", y = "Richness") + theme_bw()

  allouche_model_df <- allouche_results$allouche_model_df
  emp_long <- allouche_model_df %>%
    transmute(richness, mean_elev = mu, range = range / 100, CV = CV * 10, variance = sigma2 / 10000, delta = delta / 10) %>%
    pivot_longer(c(mean_elev, range, CV, variance, delta), names_to = "measure", values_to = "x")
  pv_audit_write_contract_csv(
    pv_audit_figure3_empirical_layers(allouche_model_df),
    "EC-018",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Figure 3 empirical plot-layer dataframe"
  )
  ec018_exact <- NULL
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    ec018_exact <- pv_audit_figure3_empirical_layers_julia_style(allouche_model_df)
    ec018_exact_path <- pv_audit_write_contract_provenance_csv(
      ec018_exact,
      "EC-018",
      "stage8lnq_exact_EC018_empirical_layers_poisson_glm_r.csv",
      subdir = "intermediate",
      object_name = "EC-018_harmonized_empirical_layers_poisson_glm",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Stage 8L-NQ exact-replication harmonized EC-018 output using Julia-style Poisson GLM fitted values"
    )
    pv_audit_write_contract_metadata_sidecar(
      tibble::tribble(
        ~output_name, ~output_sha256, ~output_row_count, ~output_column_count, ~model_family, ~link_function, ~prediction_scale, ~measure, ~degree, ~x_scale, ~key_fields, ~comparison_key_fields, ~row_id_role, ~source_function, ~provenance_only, ~manuscript_facing_output_changed, ~metadata_note,
        "stage8lnq_exact_EC018_empirical_layers_poisson_glm_r.csv", pv_audit_sha256_or_na(ec018_exact_path), nrow(ec018_exact), ncol(ec018_exact), "Poisson GLM via stats::glm", "log", "response", "range", 2L, "raw_range / 100", "measure,x,richness rounded to 8 decimals", "measure,x,richness rounded to 8 decimals", "row_id ignored for cross-language pairing; diagnostic/order-only", "pv_audit_figure3_empirical_layers_julia_style", TRUE, FALSE, "Mirrors Julia measure range=data.range/100; no offset or weights.",
        "stage8lnq_exact_EC018_empirical_layers_poisson_glm_r.csv", pv_audit_sha256_or_na(ec018_exact_path), nrow(ec018_exact), ncol(ec018_exact), "Poisson GLM via stats::glm", "log", "response", "coefficient_of_variation", 2L, "CV * 10", "measure,x,richness rounded to 8 decimals", "measure,x,richness rounded to 8 decimals", "row_id ignored for cross-language pairing; diagnostic/order-only", "pv_audit_figure3_empirical_layers_julia_style", TRUE, FALSE, "Mirrors Julia measure coefficient_of_variation=data.CV*10; no offset or weights.",
        "stage8lnq_exact_EC018_empirical_layers_poisson_glm_r.csv", pv_audit_sha256_or_na(ec018_exact_path), nrow(ec018_exact), ncol(ec018_exact), "Poisson GLM via stats::glm", "log", "response", "variance", 2L, "sigma2 / 10000", "measure,x,richness rounded to 8 decimals", "measure,x,richness rounded to 8 decimals", "row_id ignored for cross-language pairing; diagnostic/order-only", "pv_audit_figure3_empirical_layers_julia_style", TRUE, FALSE, "R source column sigma2 is recorded as the sigma_or_sigma2 alias for Julia data.sigma in the audit layer.",
        "stage8lnq_exact_EC018_empirical_layers_poisson_glm_r.csv", pv_audit_sha256_or_na(ec018_exact_path), nrow(ec018_exact), ncol(ec018_exact), "Poisson GLM via stats::glm", "log", "response", "delta", 1L, "delta / 10", "measure,x,richness rounded to 8 decimals", "measure,x,richness rounded to 8 decimals", "row_id ignored for cross-language pairing; diagnostic/order-only", "pv_audit_figure3_empirical_layers_julia_style", TRUE, FALSE, "Mirrors Julia measure delta=data.delta/10; no offset or weights."
      ),
      "EC-018",
      "stage8lnq_exact_EC018_empirical_layers_poisson_glm_metadata_r.csv",
      "exact_replication_harmonization_metadata",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Stage 8L-NQ exact-replication EC-018 harmonized output metadata"
    )
  }
  p_emp <- ggplot(emp_long, aes(x, richness)) +
    geom_point(alpha = 0.65, size = 1) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
    facet_wrap(~measure, scales = "free_x", ncol = 2) +
    labs(title = "Empirical Catalonia bird richness", x = "Predictor", y = "Richness") + theme_bw()

  diagnostic_fig <- p_theory | p_emp
  empirical_visual_layers <- if (is.null(ec018_exact)) {
    pv_audit_figure3_empirical_layers_julia_style(allouche_model_df)
  } else {
    ec018_exact
  }
	  fig_visual <- assemble_figure3_visual_replication(tg, allouche_model_df, empirical_visual_layers, root = root)
  save_png(p_theory, "Figure3_theoretical_R_equivalent.png", width = 9, height = 7)
  save_png(p_emp, "Figure3_empirical_R_equivalent.png", width = 8, height = 7)
  save_png(diagnostic_fig, "Figure3_R_diagnostic_summary.png", width = 16, height = 7)
  save_png(fig_visual, "Figure3_R_visual_replication.png", width = 11.3, height = 10.0)
  save_png(fig_visual, "Figure3_R_equivalent.png", width = 11.3, height = 10.0)
  list(theory_grid = tg, theory_models = theory_models, empirical_long = emp_long, visual_replication = fig_visual, diagnostic_summary = diagnostic_fig)
}

# 02I R reimplementation of selected SOM S1 procedures from fixedCV.jl
run_som_s1_fixedCV_R <- function(
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  mus <- c(1.5, 4, 7, 11, 15)
  cvs <- c(0.5, 1.5, 2.5)
  x <- seq(0.01, 14.9, by = 0.01)
  # Julia labels call these CV, but FAGamma(mu, alpha) sets shape = alpha, scale = mu/alpha.
  df <- expand_grid(mu = mus, alpha = cvs, x = x) %>%
    mutate(density = dgamma(x, shape = alpha, scale = mu / alpha),
           density_display = if_else(alpha == 0.5, pmin(density, 0.8), density),
           mu_lab = factor(sprintf("mu = %.1f", mu), levels = sprintf("mu = %.1f", mus)),
           cv_lab = factor(sprintf("CV = %.1f", alpha), levels = sprintf("CV = %.1f", cvs)))
  cv_colours <- c("CV = 0.5" = "#0072B2", "CV = 1.5" = "#E69F00", "CV = 2.5" = "#009E73")
  p <- ggplot(df, aes(x, density_display)) +
    geom_area(fill = "grey88", colour = NA) +
    geom_line(aes(colour = cv_lab), linewidth = 0.48, show.legend = FALSE) +
    scale_colour_manual(values = cv_colours) +
    scale_x_continuous(breaks = c(0, 5, 10, 15), expand = expansion(mult = c(0, 0.015))) +
    scale_y_continuous(breaks = function(x) pretty(x, n = 3), expand = expansion(mult = c(0, 0.035))) +
    facet_grid(cv_lab ~ mu_lab, scales = "free_y") +
    coord_cartesian(xlim = c(0, 14.9), expand = FALSE) +
    labs(x = "x", y = "Probability Density") +
    pv_visual_schematic_theme(base_size = 6.8) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "grey30", linewidth = 0.22),
      strip.text.x = element_text(face = "plain", size = 5.8, margin = margin(0.5, 1, 0.5, 1)),
      strip.text.y.right = element_text(face = "plain", size = 5.8, angle = 90, margin = margin(1, 0.5, 1, 0.5)),
      panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.28),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.13),
      panel.spacing = grid::unit(0.022, "inches"),
      axis.title = element_text(size = 6.2),
      axis.text = element_text(size = 5.2, colour = "grey25"),
      plot.margin = margin(0.6, 0.6, 0.6, 0.6)
    )
  pv_audit_write_contract_csv(
    df %>% transmute(mu, alpha_label = alpha, x, density, row_label = cv_lab),
    "EC-020",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "SOM S1 gamma density layers"
  )
  save_png(p, "SOM_s1_fixedCV_R_equivalent.png", width = 4.85, height = 2.72)
  invisible(df)
}

# 02J R reimplementation of selected SOM S2–S5 procedures from s2.jl
pv_audit_som_theory_grid_julia_style <- function(MDR_fun, HDR_fun, n = 500) {
  intercept <- 66.3947 / 2
  mu_range <- seq(64, 2422, length.out = n)
  H_range <- seq(0.64, 144, length.out = n)
  expand_grid(mu = mu_range, H = H_range) %>%
    mutate(
      MDR = MDR_fun(mu, intercept),
      HDR = HDR_fun(H, intercept),
      richness = intercept + MDR + HDR,
      range_mbh = mbh_range_gamma(mu, H) / 100,
      cv_mbh = mbh_cv_gamma(mu, H) * 10,
      var_mbh = mbh_var_gamma(mu, H) / 10000
    )
}

pv_audit_som_polyfit_terms_julia_style <- function(grid, scenario, grid_n = nrow(grid)) {
  model_columns <- c(range_mbh = "range_mbh", cv_mbh = "cv_mbh", var_mbh = "var_mbh")
  purrr::imap_dfr(model_columns, function(column, model_name) {
    model <- pv_audit_julia_polyfit_model(grid, column)
    broom::tidy(model) %>%
      transmute(
        scenario = scenario,
        model = model_name,
        term = pv_audit_julia_polyfit_term_label(term),
        estimate,
        std_error = std.error,
        statistic,
        p_value = p.value,
        grid_n = grid_n
      )
  })
}

pv_som_measure_specs <- function() {
  list(
    range = list(column = "range_mbh_visual", label = "Range (e+2)", title = "Range"),
    cv = list(column = "cv_mbh_visual", label = "Coef. Var. (e-1)", title = "Coef. Var."),
    var = list(column = "var_mbh_visual", label = "Variance (e+4)", title = "Variance")
  )
}

pv_som_visual_grid <- function(grid) {
  grid %>%
    mutate(
      range_mbh_visual = range_mbh / 100,
      cv_mbh_visual = cv_mbh * 10,
      var_mbh_visual = var_mbh / 10000
    )
}

pv_som_line_family_colours <- function(n) {
  list(
    heterogeneity = grDevices::colorRampPalette(c("#9CC5E8", "#0072B2"))(n),
    mean = grDevices::colorRampPalette(c("#F6C267", "#E69F00"))(n)
  )
}

pv_som_line_family_plot_theme <- function(base_size = 5.7) {
  pv_visual_schematic_theme(base_size = base_size) +
    theme(
      panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.36),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.13),
      axis.line = element_line(linewidth = 0.23, colour = "grey20"),
      axis.ticks = element_line(linewidth = 0.20, colour = "grey20"),
      axis.title = element_text(size = base_size - 0.55),
      axis.text = element_text(size = base_size - 1.35, colour = "grey20"),
      plot.title = element_text(size = base_size - 0.15, face = "plain", hjust = 0.5, margin = margin(0, 0, 0.8, 0)),
      plot.margin = margin(0.45, 0.45, 0.45, 0.45)
    )
}

pv_som_panel_letter_layer <- function(label, size = 2.25) {
  if (is.null(label) || !nzchar(label)) {
    return(NULL)
  }
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = label,
    hjust = -0.35,
    vjust = 1.18,
    fontface = "bold",
    size = size,
    colour = "grey10"
  )
}

pv_som_true_relationship_panel <- function(MDR_fun, HDR_fun, x_role = c("mean", "heterogeneity"), n = 180, nlines = 7) {
  x_role <- match.arg(x_role)
  intercept <- 66.3947 / 2
  mu_range <- seq(64, 2422, length.out = n)
  H_range <- seq(0.64, 144, length.out = n)
  mu_lines <- seq(64, 2422, length.out = nlines)
  H_lines <- seq(0.64, 144, length.out = nlines)
  cols <- pv_som_line_family_colours(nlines)

  fixed_H <- purrr::map2_dfr(H_lines, seq_along(H_lines), function(H_value, idx) {
    if (identical(x_role, "mean")) {
      tibble(
        x = mu_range,
        y = intercept + MDR_fun(mu_range, intercept) + HDR_fun(H_value, intercept),
        family = paste0("H ", idx),
        family_role = "Heterogeneity",
        colour = cols$heterogeneity[[idx]],
        linetype = "solid"
      )
    } else {
      tibble(
        x = rep(H_value, length(mu_range)),
        y = intercept + MDR_fun(mu_range, intercept) + HDR_fun(H_value, intercept),
        family = paste0("H ", idx),
        family_role = "Heterogeneity",
        colour = cols$heterogeneity[[idx]],
        linetype = "solid"
      )
    }
  })

  fixed_mu <- purrr::map2_dfr(mu_lines, seq_along(mu_lines), function(mu_value, idx) {
    if (identical(x_role, "mean")) {
      tibble(
        x = rep(mu_value, length(H_range)),
        y = intercept + MDR_fun(mu_value, intercept) + HDR_fun(H_range, intercept),
        family = paste0("Mean ", idx),
        family_role = "Mean",
        colour = cols$mean[[idx]],
        linetype = "22"
      )
    } else {
      tibble(
        x = H_range,
        y = intercept + MDR_fun(mu_value, intercept) + HDR_fun(H_range, intercept),
        family = paste0("Mean ", idx),
        family_role = "Mean",
        colour = cols$mean[[idx]],
        linetype = "22"
      )
    }
  })

  bind_rows(fixed_H, fixed_mu)
}

pv_som_line_family_plot <- function(data, x_label, y_label, title, tag = NULL) {
  tag <- if (is.null(tag)) "" else tag
  ggplot(data, aes(x, y, group = interaction(family_role, family), colour = colour, linetype = linetype)) +
    geom_line(aes(linetype = linetype), linewidth = 0.72, colour = "grey70", alpha = 0.45, show.legend = FALSE) +
    geom_line(linewidth = 0.42, show.legend = FALSE) +
    pv_som_panel_letter_layer(tag) +
    scale_colour_identity() +
    scale_linetype_identity() +
    labs(x = x_label, y = y_label, title = title) +
    pv_som_line_family_plot_theme(base_size = 5.55)
}

pv_som_mbh_line_data <- function(grid, measure_column, panel_role = c("mean_heterogeneity", "observed_hdr"), nlines = 7) {
  panel_role <- match.arg(panel_role)
  mu_lines <- seq(min(grid$mu, na.rm = TRUE), max(grid$mu, na.rm = TRUE), length.out = nlines)
  H_lines <- seq(min(grid$H, na.rm = TRUE), max(grid$H, na.rm = TRUE), length.out = nlines)
  cols <- pv_som_line_family_colours(nlines)

  fixed_H <- purrr::map2_dfr(H_lines, seq_along(H_lines), function(H_value, idx) {
    rows <- grid %>%
      filter(abs(H - H_value) == min(abs(H - H_value), na.rm = TRUE)) %>%
      arrange(mu)
    (if (identical(panel_role, "mean_heterogeneity")) {
      tibble(x = rows$mu, y = rows[[measure_column]])
    } else {
      tibble(x = rows[[measure_column]], y = rows$richness)
    }) %>%
      mutate(family = paste0("H ", idx), family_role = "Heterogeneity", colour = cols$heterogeneity[[idx]], linetype = "solid")
  })

  fixed_mu <- purrr::map2_dfr(mu_lines, seq_along(mu_lines), function(mu_value, idx) {
    rows <- grid %>%
      filter(abs(mu - mu_value) == min(abs(mu - mu_value), na.rm = TRUE)) %>%
      arrange(H)
    (if (identical(panel_role, "mean_heterogeneity")) {
      tibble(x = rows$mu, y = rows[[measure_column]])
    } else {
      tibble(x = rows[[measure_column]], y = rows$richness)
    }) %>%
      mutate(family = paste0("Mean ", idx), family_role = "Mean", colour = cols$mean[[idx]], linetype = "22")
  })

  bind_rows(fixed_H, fixed_mu)
}

pv_som_polyfit_panel <- function(grid, measure_column, measure_label, tag = NULL, title = NULL) {
  tag <- if (is.null(tag)) "" else tag
  fit_df <- grid %>%
    transmute(x = .data[[measure_column]], richness) %>%
    filter(is.finite(x), is.finite(richness))
  if (!nrow(fit_df)) {
    return(pv_visual_empty_panel("Polynomial fit", measure_label) + pv_som_panel_letter_layer(tag))
  }
  model <- stats::lm(richness ~ x + I(x^2), data = fit_df)
  x_range <- range(fit_df$x, finite = TRUE)
  curve_df <- tibble(x = seq(x_range[1], x_range[2], length.out = 160))
  curve_df$richness <- as.numeric(stats::predict(model, newdata = curve_df))
  co <- stats::coef(model)
  label <- sprintf("%.1f %+.2g x %+.2g x^2", co[[1]], co[[2]], co[[3]])
  y_range <- range(curve_df$richness, finite = TRUE)
  y_span <- diff(y_range)

  ggplot(curve_df, aes(x, richness)) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = y_range[2] - 0.18 * y_span,
      ymax = Inf,
      fill = "grey94",
      colour = NA
    ) +
    geom_line(linewidth = 0.52, colour = "grey18") +
    pv_som_panel_letter_layer(tag) +
    annotate(
      "label",
      x = x_range[1] + 0.05 * diff(x_range),
      y = y_range[2] - 0.08 * y_span,
      label = label,
      size = 1.08,
      linewidth = 0.10,
      fill = "grey96",
      hjust = 0,
      vjust = 0.5,
      label.padding = grid::unit(0.05, "lines")
    ) +
    labs(x = measure_label, y = NULL, title = title) +
    pv_som_line_family_plot_theme(base_size = 5.55) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

pv_som_line_family_guide <- function() {
  cols <- pv_som_line_family_colours(7)
  bars <- bind_rows(
    tibble(x = 0.32, y = seq_along(cols$heterogeneity), colour = cols$heterogeneity, guide = "Heterogeneity"),
    tibble(x = 0.70, y = seq_along(cols$mean), colour = cols$mean, guide = "Mean (e+3)")
  )
  ggplot(bars, aes(x, y)) +
    geom_tile(aes(fill = colour), width = 0.10, height = 0.90, colour = "grey25", linewidth = 0.06) +
    annotate("segment", x = 0.62, xend = 0.74, y = 1, yend = 1, linetype = "22", linewidth = 0.35, colour = "#E69F00") +
    annotate("segment", x = 0.62, xend = 0.74, y = 7, yend = 7, linetype = "22", linewidth = 0.35, colour = "#E69F00") +
    annotate("text", x = 0.32, y = 7.47, label = "Heterogeneity", angle = 90, size = 1.35, fontface = "bold") +
    annotate("text", x = 0.68, y = 7.47, label = "Mean (e+3)", angle = 90, size = 1.35, fontface = "bold") +
    annotate("text", x = 0.19, y = c(1, 3, 5, 7), label = c("1", "37", "73", "144"), size = 1.25, hjust = 1) +
    annotate("text", x = 0.82, y = c(1, 3.5, 6), label = c("0.1", "1.3", "1.9"), size = 1.25, hjust = 0) +
    scale_fill_identity() +
    xlim(0, 1) +
    ylim(0.55, 7.85) +
    theme_void(base_size = 4.8) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

pv_som_full_revised_figure <- function(name, grid, MDR_fun, HDR_fun) {
  vgrid <- pv_som_visual_grid(grid)
  specs <- pv_som_measure_specs()

  top_row <- (
    pv_som_line_family_plot(
      pv_som_true_relationship_panel(MDR_fun, HDR_fun, x_role = "mean"),
      "Mean",
      "Richness",
      "Mean-Diversity (MDR)",
      "a"
    ) |
      pv_som_line_family_plot(
        pv_som_true_relationship_panel(MDR_fun, HDR_fun, x_role = "heterogeneity"),
        "Heterogeneity",
        "Richness",
        "Heterogeneity-Diversity (HDR)",
        "b"
      ) |
      pv_som_line_family_guide()
  ) + patchwork::plot_layout(widths = c(1, 1, 0.22))

  row_tags <- list(
    range = c("c", "d", "e"),
    cv = c("f", "g", "h"),
    var = c("i", "j", "k")
  )
  measure_rows <- purrr::imap(specs, function(spec, measure_name) {
    tags <- row_tags[[measure_name]]
    left_title <- if (identical(measure_name, "range")) "Mean-heterogeneity" else NULL
    middle_title <- if (identical(measure_name, "range")) "Observed mean-biased HDRs" else NULL
    (
      pv_som_line_family_plot(
        pv_som_mbh_line_data(vgrid, spec$column, "mean_heterogeneity"),
        "Mean",
        spec$label,
        left_title,
        tags[[1]]
      ) |
        pv_som_line_family_plot(
          pv_som_mbh_line_data(vgrid, spec$column, "observed_hdr"),
          spec$label,
          "Richness",
          middle_title,
          tags[[2]]
        ) |
        pv_som_polyfit_panel(vgrid, spec$column, spec$label, tags[[3]])
    ) + patchwork::plot_layout(widths = c(1, 1, 0.72))
  })

  (top_row / measure_rows$range / measure_rows$cv / measure_rows$var) +
    patchwork::plot_layout(heights = c(0.58, 0.95, 0.95, 0.95)) +
    patchwork::plot_annotation(
      theme = theme(
        plot.background = element_rect(fill = NA, colour = "grey25", linewidth = 0.35),
        plot.margin = margin(1, 1, 1, 1)
      )
    ) &
    theme(plot.margin = margin(0.55, 0.55, 0.55, 0.55))
}

make_som_theory_scenario <- function(
  name,
  MDR_fun,
  HDR_fun,
  n = 250,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  intercept <- 66.3947 / 2
  mu_range <- seq(64, 2422, length.out = n)
  H_range <- seq(0.64, 144, length.out = n)
  grid <- expand_grid(mu = mu_range, H = H_range) %>%
    mutate(MDR = MDR_fun(mu, intercept), HDR = HDR_fun(H, intercept), richness = intercept + MDR + HDR,
           range_mbh = mbh_range_gamma(mu, H), cv_mbh = mbh_cv_gamma(mu, H), var_mbh = mbh_var_gamma(mu, H))
  long <- grid %>%
    sample_n(min(20000, nrow(grid))) %>%
    pivot_longer(c(range_mbh, cv_mbh, var_mbh), names_to = "measure", values_to = "heterogeneity")
  p <- ggplot(long, aes(heterogeneity, richness)) + geom_point(alpha = 0.05, size = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
    facet_wrap(~measure, scales = "free_x") +
    labs(title = paste0("SOM ", name, " R equivalent"), x = "Mean-biased heterogeneity", y = "Richness") + theme_bw()
  save_png(p, paste0("SOM_", name, "_R_equivalent.png"), width = 8, height = 5)
  p_full <- pv_som_full_revised_figure(name, grid, MDR_fun, HDR_fun)
  save_png(p_full, paste0("SOM_", name, "_full_R_equivalent.png"), width = 5.0, height = 6.25)
  models <- list(
    range = lm(richness ~ range_mbh + I(range_mbh^2), data = grid),
    cv = lm(richness ~ cv_mbh + I(cv_mbh^2), data = grid),
    var = lm(richness ~ var_mbh + I(var_mbh^2), data = grid)
  )
  terms <- imap_dfr(models, ~ broom::tidy(.x) %>% mutate(scenario = name, model = .y, .before = 1))
  readr::write_csv(terms, safe_file(tab_dir, paste0("SOM_", name, "_polyfits_R.csv")))
  list(data = grid, plot = p, full_plot = p_full, models = models, terms = terms, grid_n = n)
}

run_som_s2_s5_R <- function(
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  scenarios <- list(
    s2r = list(MDR = function(mu, intercept) intercept + mu * 0.03 + mu^2 * -1.2e-5,
               HDR = function(H, intercept) H * 0.21 + H^2 * -1.2e-3),
    s3r = list(MDR = function(mu, intercept) intercept + mu * 0.004,
               HDR = function(H, intercept) H * -0.035),
    s4r = list(MDR = function(mu, intercept) intercept + mu * 0.009,
               HDR = function(H, intercept) H * 0.11),
    s5r = list(MDR = function(mu, intercept) intercept + mu * 0.015,
               HDR = function(H, intercept) rep(-5, length(H)))
  )
  out <- imap(scenarios, ~ make_som_theory_scenario(
    .y,
    .x$MDR,
    .x$HDR,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_audit_intermediates = write_audit_intermediates
  ))
  pv_audit_write_contract_csv(
    purrr::imap_dfr(out, function(result, scenario) {
      result$data %>%
        mutate(scenario = scenario, grid_n = result$grid_n, .before = 1) %>%
        transmute(scenario, mu, H, MDR, HDR, richness, range_mbh, cv_mbh, var_mbh, grid_n)
    }),
    "EC-021",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "SOM scenario grids"
  )
  # Stage 8R-NQ instrumentation: EC-021 metadata sidecar for units and scales only.
  pv_audit_write_contract_metadata_sidecar(
    purrr::map_dfr(names(out), function(scenario) {
      grid_n <- out[[scenario]]$grid_n
      tibble::tribble(
        ~column_name, ~column_role, ~unit_or_scale, ~scale_basis, ~formula_source, ~grid_n, ~metadata_note,
        "scenario", "scenario_label", "s2r-s5r deterministic scenario", "scenario list name", "run_som_s2_s5_R scenarios list", grid_n, "No primary EC-021 values are changed.",
        "mu", "input_axis", "original mu grid units", "seq(64, 2422, length.out = grid_n)", "make_som_theory_scenario", grid_n, "No primary EC-021 values are changed.",
        "H", "input_axis", "original H grid units", "seq(0.64, 144, length.out = grid_n)", "make_som_theory_scenario", grid_n, "No primary EC-021 values are changed.",
        "MDR", "component", "deterministic MDR contribution", "scenario-specific MDR function", "MDR_fun(mu, intercept)", grid_n, "No primary EC-021 values are changed.",
        "HDR", "component", "deterministic HDR contribution", "scenario-specific HDR function", "HDR_fun(H, intercept)", grid_n, "No primary EC-021 values are changed.",
        "richness", "response", "predicted richness", "intercept + MDR + HDR", "richness expression in make_som_theory_scenario", grid_n, "No primary EC-021 values are changed.",
        "range_mbh", "computed_measure", "mbh_range_gamma output", "0.01/0.99 gamma quantile range helper output", "mbh_range_gamma(mu, H)", grid_n, "No primary EC-021 values are changed.",
        "cv_mbh", "computed_measure", "mbh_cv_gamma output", "gamma coefficient-of-variation helper output", "mbh_cv_gamma(mu, H)", grid_n, "No primary EC-021 values are changed.",
        "var_mbh", "computed_measure", "mbh_var_gamma output", "gamma variance helper output", "mbh_var_gamma(mu, H)", grid_n, "No primary EC-021 values are changed.",
        "grid_n", "metadata", "row-generation grid size", "n passed to make_som_theory_scenario", "grid_n", grid_n, "No primary EC-021 values are changed."
      ) %>%
        mutate(scenario = scenario, .before = 1)
    }),
    "EC-021",
    "stage8rnq_patch_EC021_unit_metadata_r.csv",
    "unit_annotation",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-021 unit/scale metadata sidecar"
  )
  pv_audit_write_contract_csv(
    purrr::imap_dfr(out, function(result, scenario) {
      result$terms %>%
        transmute(
          scenario,
          model,
          term,
          estimate,
          std_error = std.error,
          statistic,
          p_value = p.value,
          grid_n = result$grid_n
        )
    }),
    "EC-022",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "SOM scenario model terms"
  )
	  # Stage 8R-NQ instrumentation: EC-022 metadata sidecar for model statistics only.
	  pv_audit_write_contract_metadata_sidecar(
	    purrr::imap_dfr(out, function(result, scenario) {
	      purrr::imap_dfr(result$models, function(model, model_name) {
	        # Stage 8R-NQ EC-022 sidecar: prevent the model column from being
	        # shadowed during tibble construction.
	        formula_obj <- stats::formula(model)
	        model_formula <- paste(deparse(formula_obj), collapse = " ")
	        predictor_column <- all.vars(formula_obj)[[2]]
	        tibble::tibble(
	          scenario = scenario,
	          model = model_name,
	          model_formula = model_formula,
	          fit_engine = "stats::lm",
	          predictor_column = predictor_column,
	          term_convention = "linear term plus I(predictor^2)",
	          df_source = "stats::df.residual(model)",
	          statistic_source = "broom::tidy statistic",
          std_error_source = "broom::tidy std.error",
          p_value_source = "broom::tidy p.value",
          model_family = "Gaussian linear model via stats::lm",
          grid_n = result$grid_n
        )
      })
    }),
    "EC-022",
    "stage8rnq_patch_EC022_model_statistic_metadata_r.csv",
    "model_statistic_metadata",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-022 model-statistic metadata sidecar"
  )
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    # This harmonized export, used only for provenance, mirrors the grid density
    # and coefficient term convention of the P&V Julia audit for EC-022.
    som_exact_grids <- purrr::imap(scenarios, function(spec, scenario) {
      pv_audit_som_theory_grid_julia_style(spec$MDR, spec$HDR, n = 500)
    })
    ec022_exact <- purrr::imap_dfr(som_exact_grids, function(grid, scenario) {
      pv_audit_som_polyfit_terms_julia_style(grid, scenario = scenario, grid_n = nrow(grid))
    })
    ec022_exact_path <- pv_audit_write_contract_provenance_csv(
      ec022_exact,
      "EC-022",
      "stage8lnq_exact_EC022_som_polyfit_terms_grid250000_r.csv",
      subdir = "tables",
      object_name = "EC-022_harmonized_som_polyfit_terms_grid250000",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Stage 8L-NQ exact-replication harmonized EC-022 output using Julia-style grid_n=250000"
    )
    pv_audit_write_contract_metadata_sidecar(
      tibble::tibble(
        output_name = "stage8lnq_exact_EC022_som_polyfit_terms_grid250000_r.csv",
        output_sha256 = pv_audit_sha256_or_na(ec022_exact_path),
        output_row_count = nrow(ec022_exact),
        output_column_count = ncol(ec022_exact),
        scenario_count = length(som_exact_grids),
        scenarios = paste(names(som_exact_grids), collapse = ";"),
        model_family = "Gaussian linear model via stats::lm",
        model_basis = "intercept + mbh + mbh^2",
        coefficient_terms = "intercept;mbh;mbh_squared",
        grid_axis_n = 500L,
        grid_n = unique(vapply(som_exact_grids, nrow, integer(1))),
        mbh_scale_convention = "range_mbh=mbh_range_gamma/100; cv_mbh=mbh_cv_gamma*10; var_mbh=mbh_var_gamma/10000",
        key_fields = "scenario,model,term",
        source_function = "pv_audit_som_polyfit_terms_julia_style",
        provenance_only = TRUE,
        manuscript_facing_output_changed = FALSE,
        metadata_note = "grid_n is total rows per scenario, matching Julia length(mu)=500*500; MBH predictors are scaled to Julia exact-replication units; no Julia output values are hard-coded."
      ),
      "EC-022",
      "stage8lnq_exact_EC022_som_polyfit_terms_grid250000_metadata_r.csv",
      "exact_replication_harmonization_metadata",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Stage 8L-NQ exact-replication EC-022 harmonized output metadata"
    )
  }
  out
}

# 02K R reimplementation of selected SOM S6 procedures from somfig.jl
run_som_s6_delta_location_R <- function(
  adf,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  if (!all(c("X_coord", "Y_coord") %in% names(adf))) {
    warning("Allouche data lacks X_coord/Y_coord; skipping SOM S6.")
    return(NULL)
  }
  s6_theme <- pv_visual_schematic_theme(base_size = 6.8) +
    theme(
      panel.border = element_rect(fill = NA, colour = "grey25", linewidth = 0.34),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.14),
      axis.title = element_text(size = 6.6),
      axis.text = element_text(size = 5.7, colour = "grey25"),
      plot.margin = margin(0.9, 0.9, 0.9, 0.9)
    )
  p1 <- ggplot(adf, aes(delta, X_coord / 100000)) +
    geom_point(alpha = 0.82, size = 0.72, colour = "#1F77B4") +
    labs(x = expression(delta), y = "X coord (e+5)") +
    s6_theme
  p2 <- ggplot(adf, aes(delta, Y_coord / 1000000)) +
    geom_point(alpha = 0.82, size = 0.72, colour = "#1F77B4") +
    labs(x = expression(delta), y = "Y coord (e+6)") +
    s6_theme
  pv_audit_write_contract_csv(
    adf %>%
      mutate(row_id = dplyr::row_number()) %>%
      transmute(row_id, delta, X_coord, Y_coord, X_scaled = X_coord / 100000, Y_scaled = Y_coord / 1000000),
    "EC-023",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "SOM S6 plotted dataframe"
  )
  save_png(p1 | p2, "SOM_s6_delta_location_R_equivalent.png", width = 4.85, height = 2.25)
  p1 | p2
}


coef_equivalence_test_row <- function(name, delta_model, comparison_model) {
  beta <- abs(coef(delta_model)[2]); se <- summary(delta_model)$coefficients[2, "Std. Error"]
  null <- abs(coef(comparison_model)[2]) / 2; df <- stats::df.residual(delta_model)
  t_stat <- (beta - null) / se
  p_two <- 2 * pt(-abs(t_stat), df = df)
  tibble(test = name, delta_abs_coef = beta, delta_se = se, comparison_abs_coef = abs(coef(comparison_model)[2]),
         null_half_comparison = null, t_statistic = t_stat, df = df, p_value_two_sided = p_two,
         abs_delta_over_comparison = beta / abs(coef(comparison_model)[2]))
}
run_delta_trend_test_R <- function(
  root = repo_root,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  elev <- read_elevation_heterogeneity(root) %>% mutate(
    rsvar = rescale01(elev_var), rssd = rescale01(sd), rsmean = rescale01(mean), rsiod = rescale01(iod),
    rscv = rescale01(cv), sqrtmean = sqrt(mean), sqrtmeandiv = sqrt(mean) / mean
  )
  crop <- read_crop_heterogeneity(root) %>% mutate(
    rsvar = rescale01(crop_var), rssd = rescale01(sd), rsmean = rescale01(mean), rsdelta = rescale01(delta),
    rscv = rescale01(cv), x = mean * (100 - mean), sqrtx = sqrt(mean * (100 - mean)), sqrtxdiv = sqrt(mean * (100 - mean)) / mean
  )
  rows <- bind_rows(
    coef_equivalence_test_row("elevation_var_vs_delta", lm(rsiod ~ mean, elev), lm(rsvar ~ mean, elev)),
    coef_equivalence_test_row("elevation_sd_vs_delta", lm(rsiod ~ sqrtmean, elev), lm(rssd ~ sqrtmean, elev)),
    coef_equivalence_test_row("elevation_cv_vs_delta", lm(rsiod ~ sqrtmeandiv, elev), lm(rscv ~ sqrtmeandiv, elev)),
    coef_equivalence_test_row("crop_var_vs_delta", lm(rsdelta ~ x, crop), lm(rsvar ~ x, crop)),
    coef_equivalence_test_row("crop_sd_vs_delta", lm(rsdelta ~ sqrtx, crop), lm(rssd ~ sqrtx, crop)),
    coef_equivalence_test_row("crop_cv_vs_delta", lm(rsdelta ~ sqrtxdiv, crop), lm(rscv ~ sqrtxdiv, crop))
  )
  pv_audit_write_contract_csv(
    bind_rows(
      elev %>%
        mutate(row_id = dplyr::row_number()) %>%
        transmute(
          dataset = "elevation",
          row_id,
          mean,
          rsvar,
          rssd,
          rsmean,
          rsiod,
          rsdelta = NA_real_,
          rscv,
          sqrtmean,
          sqrtmeandiv,
          x = NA_real_,
          sqrtx = NA_real_,
          sqrtxdiv = NA_real_
        ),
      crop %>%
        mutate(row_id = dplyr::row_number()) %>%
        transmute(
          dataset = "crop_cover",
          row_id,
          mean,
          rsvar,
          rssd,
          rsmean,
          rsiod = NA_real_,
          rsdelta,
          rscv,
          sqrtmean = NA_real_,
          sqrtmeandiv = NA_real_,
          x,
          sqrtx,
          sqrtxdiv
        )
    ),
    "EC-029",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Delta trend model inputs"
  )
  pv_audit_write_contract_csv(
    rows %>%
      transmute(
        test,
        delta_abs_coef,
        delta_se,
        comparison_abs_coef,
        null_half_comparison,
        t_statistic,
        df,
        p_value = p_value_two_sided,
        ratio = abs_delta_over_comparison
      ),
    "EC-030",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Delta trend test statistics"
  )
  # Stage 8R-NQ instrumentation: EC-030 metadata sidecar for test statistics only.
  pv_audit_write_contract_metadata_sidecar(
    tibble::tibble(
      test = rows$test,
      delta_model_formula = c(
        "rsiod ~ mean",
        "rsiod ~ sqrtmean",
        "rsiod ~ sqrtmeandiv",
        "rsdelta ~ x",
        "rsdelta ~ sqrtx",
        "rsdelta ~ sqrtxdiv"
      ),
      comparison_model_formula = c(
        "rsvar ~ mean",
        "rssd ~ sqrtmean",
        "rscv ~ sqrtmeandiv",
        "rsvar ~ x",
        "rssd ~ sqrtx",
        "rscv ~ sqrtxdiv"
      ),
      beta_rule = "abs(coef(delta_model)[2])",
      se_rule = "summary(delta_model)$coefficients[2, 'Std. Error']",
      null_rule = "abs(coef(comparison_model)[2]) / 2",
      df_rule = "stats::df.residual(delta_model)",
      t_statistic_rule = "(beta - null) / se",
      p_value_rule = "2 * pt(-abs(t_stat), df = df)",
      sidedness = "two_sided",
      ratio_rule = "beta / abs(coef(comparison_model)[2])"
    ),
    "EC-030",
    "stage8rnq_patch_EC030_test_statistic_metadata_r.csv",
    "model_statistic_metadata",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-030 test-statistic metadata sidecar"
  )
  write_csv(rows, safe_file(tab_dir, "delta_negligible_trend_tests_R.csv"))
  rows
}

# 02N 06_review2/2.1.jl -> review response figures
pv_review2_beta_density_panel <- function(p, q, panel_label, tag = NULL) {
  tag <- if (is.null(tag)) "" else tag
  x <- seq(0.001, 0.999, by = 0.001)
  mu <- p / (p + q)
  delta <- 1 / (p + q + 1)
  phi <- p + q
  df <- tibble(x = x, density = stats::dbeta(x, p, q))
  y_max <- max(df$density, na.rm = TRUE)
  label_x <- if (mu <= 0.5) 0.94 else 0.06
  label_hjust <- if (mu <= 0.5) 1 else 0

  ggplot(df, aes(x, density)) +
    geom_area(fill = "grey86", colour = NA) +
    geom_line(linewidth = 0.42, colour = "grey10") +
    geom_vline(xintercept = mu, linetype = "22", linewidth = 0.35, colour = "grey25") +
    annotate(
      "text",
      x = label_x,
      y = y_max * 0.92,
      hjust = label_hjust,
      vjust = 1,
      size = 2.1,
      label = paste0(
        "p = ", p, ", q = ", q, "\n",
        "mu = ", signif(mu, 3), ", delta = ", signif(delta, 3), "\n",
        "phi = ", phi
      )
    ) +
    scale_x_continuous(
      breaks = c(0, 1 / 3, 2 / 3, 1),
      labels = c("0", "1/3", "2/3", "1"),
      limits = c(0, 1)
    ) +
    labs(x = "x", y = "Density", title = panel_label) +
    pv_visual_schematic_theme(base_size = 6) +
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.02)) %>%
    pv_visual_panel_label(tag, size = 6.5)
}

pv_review2_delta_mu_panel <- function(kind = c("problem", "solution")) {
  kind <- match.arg(kind)
  mu_grid <- seq(0, 1, by = 0.001)
  if (identical(kind, "problem")) {
    curve_df <- tibble(mu = mu_grid, delta = mu_grid / (2 + mu_grid))
    point_df <- tibble(
      mu = c(1 / 3, 2 / 3),
      delta = c((1 / 3) / (2 + 1 / 3), (2 / 3) / (2 + 2 / 3)),
      label = c("B1", "B2")
    )
    title <- "Delta changes with mean"
  } else {
    phi <- 6
    curve_df <- tibble(mu = mu_grid, delta = rep(1 / (phi + 1), length(mu_grid)))
    point_df <- tibble(mu = c(1 / 3, 2 / 3), delta = rep(1 / (phi + 1), 2), label = c("B1", "B2"))
    title <- "Delta held constant"
  }

  ggplot(curve_df, aes(mu, delta)) +
    geom_line(linewidth = 0.46, colour = "#D55E00") +
    geom_point(data = point_df, aes(mu, delta), inherit.aes = FALSE, size = 1.6, colour = "black") +
    geom_text(data = point_df, aes(mu, delta, label = label), inherit.aes = FALSE, nudge_y = -0.018, size = 2.1) +
    scale_x_continuous(
      breaks = c(0, 1 / 3, 2 / 3, 1),
      labels = c("0", "1/3", "2/3", "1"),
      limits = c(0, 1)
    ) +
    coord_cartesian(ylim = c(0, 0.32), expand = FALSE) +
    labs(x = "mu", y = "delta", title = title) +
    pv_visual_schematic_theme(base_size = 6) +
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5))
}

pv_review2_formula_panel <- function(kind = c("problem", "solution")) {
  kind <- match.arg(kind)
  if (identical(kind, "problem")) {
    lines <- c("delta = mu / (p + mu)", "if p = 2,", "q in (0, infinity)")
  } else {
    lines <- c("delta = mu / (p + mu)", "if phi = 6,", "mu in (0, 1)")
  }
  tibble(x = 0, y = 0) %>%
    ggplot(aes(x, y)) +
    annotate("segment", x = 0.08, xend = 0.34, y = 0.78, yend = 0.78, colour = "#D55E00", linewidth = 0.55) +
    annotate("text", x = 0.08, y = 0.57, hjust = 0, label = paste(lines, collapse = "\n"), size = 2.4) +
    annotate("text", x = 0.08, y = 0.20, hjust = 0, label = "Relationship guide", fontface = "bold", size = 2.0) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_size = 6)
}

pv_review2_full_problem_solution_plot <- function(kind = c("problem", "solution")) {
  kind <- match.arg(kind)
  if (identical(kind, "problem")) {
    left <- pv_review2_beta_density_panel(2, 4, "B1", "a")
    right <- pv_review2_beta_density_panel(2, 1, "B2", "b")
    title <- "Review 2.1 problem full R equivalent"
  } else {
    left <- pv_review2_beta_density_panel(2, 4, "B1", "a")
    right <- pv_review2_beta_density_panel(4, 2, "B2", "b")
    title <- "Review 2.1 solution full R equivalent"
  }
  bottom <- (
    pv_visual_panel_label(pv_review2_delta_mu_panel(kind), "c", size = 6.5) |
      pv_review2_formula_panel(kind)
  ) + patchwork::plot_layout(widths = c(0.62, 0.38))

  ((left | right) / bottom) +
    patchwork::plot_layout(heights = c(1, 1.05)) +
    patchwork::plot_annotation(title = title)
}

run_review2_figures_R <- function(
  nsim = REVIEW2_SIM_N,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_audit_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  x <- seq(0.001, 0.999, by = 0.001)
  make_beta_density_plot <- function(p1, q1, p2, q2, filename, title) {
    df <- bind_rows(
      tibble(x = x, density = dbeta(x, p1, q1), panel = paste0("p=", p1, ", q=", q1), p = p1, q = q1),
      tibble(x = x, density = dbeta(x, p2, q2), panel = paste0("p=", p2, ", q=", q2), p = p2, q = q2)
    )
    p <- ggplot(df, aes(x, density)) + geom_area(fill = "grey85") + geom_line() +
      facet_wrap(~panel) + labs(title = title, x = "x", y = "Density") + theme_bw()
    save_png(p, filename, width = 7, height = 4)
    list(plot = p, data = df)
  }
  problem <- make_beta_density_plot(2,4,2,1,"review2_2.1problem_R_equivalent.png", "Review 2.1 problem")
  solution <- make_beta_density_plot(2,4,4,2,"review2_2.1solution_R_equivalent.png", "Review 2.1 solution")
  p_problem <- problem$plot
  p_solution <- solution$plot
  p_problem_full <- pv_review2_full_problem_solution_plot("problem")
  p_solution_full <- pv_review2_full_problem_solution_plot("solution")
  save_png(p_problem_full, "review2_2.1problem_full_R_equivalent.png", width = 7, height = 6.2)
  save_png(p_solution_full, "review2_2.1solution_full_R_equivalent.png", width = 7, height = 6.2)
  if (pv_audit_enabled(audit_mode, write_audit_intermediates)) {
    pv_audit_write_contract_csv(
      bind_rows(
        mutate(problem$data, figure = "problem", .before = 1),
        mutate(solution$data, figure = "solution", .before = 1)
      ) %>%
        mutate(
          mu = p / (p + q),
          delta = 1 / (p + q + 1),
          phi = p + q,
          curve_x = x,
          curve_y = density
        ) %>%
        transmute(figure, panel, p, q, x, density, mu, delta, phi, curve_x, curve_y),
      "EC-031",
      audit_paths = audit_paths,
      audit_mode = audit_mode,
      write_intermediates = write_audit_intermediates,
      notes = "Review2 deterministic beta layer data"
    )
  }

  write_audit_seed_record(
    "review2_pq_and_muvar_simulations",
    12345,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = paste0("nsim=", nsim)
  )
  set.seed(12345)
  review2_rng_kind_after_seed <- paste(RNGkind(), collapse = ";")
  p <- rexp(nsim, rate = 1 / 10); q <- rexp(nsim, rate = 1 / 10)
  mu <- p / (p + q); delta <- 1 / (p + q + 1); varb <- (p * q) / ((p + q)^2 * (p + q + 1))
  pqdf <- tibble(mu = mu, delta = delta, var = varb)
  pv_audit_write_contract_csv(
    pv_audit_review2_summary(
      pqdf,
      simulation = "pq_exponential",
      seed_status = "set.seed(12345)",
      distribution = "p,q ~ rexp(rate = 1/10)",
      notes = paste0("nsim=", nsim)
    ),
    "EC-032",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Review2 p/q simulation summary"
  )
  # Stage 8R-NQ instrumentation: EC-032 metadata sidecar for stochastic settings and seeds only.
  pv_audit_write_contract_metadata_sidecar(
    tibble::tibble(
      simulation = "pq_exponential",
      seed_value = 12345L,
      seed_call = "set.seed(12345)",
      rng_kind = review2_rng_kind_after_seed,
      draw_order = "after set.seed and before EC-033 runif stream",
      draw_count = paste0("p=", nsim, ";q=", nsim),
      distribution = "p,q ~ rexp(rate = 1/10)",
      summary_policy = "pv_audit_review2_summary means and mu quantiles",
      nsim = nsim,
      stream_dependency = "starts the review2_pq_and_muvar_simulations seeded stream"
    ),
    "EC-032",
    "stage8rnq_patch_EC032_stochastic_seed_metadata_r.csv",
    "stochastic_seed_metadata",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-032 stochastic/seed metadata sidecar"
  )
  p_pq <- (ggplot(pqdf, aes(mu, var)) + geom_hex(bins = 120) + scale_fill_viridis_c(trans="log10") + geom_smooth(se=FALSE, color="red") + theme_bw() + labs(x="mu", y="sigma^2")) |
    (ggplot(pqdf, aes(mu, delta)) + geom_hex(bins = 120) + scale_fill_viridis_c(trans="log10") + geom_smooth(se=FALSE, color="red") + theme_bw() + labs(x="mu", y="delta"))
  save_png(p_pq, "review2_2.1pqsolution_R_equivalent.png", width = 9, height = 4)

  mu2 <- runif(nsim); maxvar <- mu2 * (1 - mu2); var2 <- runif(nsim, 0, maxvar); delta2 <- var2 / (mu2 * (1 - mu2))
  mudf <- tibble(mu = mu2, var = var2, delta = delta2)
  pv_audit_write_contract_csv(
    pv_audit_review2_summary(
      mudf,
      simulation = "mu_variance_uniform",
      seed_status = "set.seed(12345)",
      notes = paste0("nsim=", nsim, "; mu~runif; var~runif(0, mu*(1-mu))")
    ),
    "EC-033",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Review2 mu/variance simulation summary"
  )
  # Stage 8R-NQ instrumentation: EC-033 metadata sidecar for stochastic settings and seeds only.
  pv_audit_write_contract_metadata_sidecar(
    tibble::tibble(
      simulation = "mu_variance_uniform",
      seed_value = 12345L,
      seed_call = "set.seed(12345) earlier in run_review2_figures_R",
      rng_kind = review2_rng_kind_after_seed,
      rng_stream_dependency = "continues after EC-032 rexp draws",
      draw_order = "mu2 runif then var2 conditional runif",
      draw_count = paste0("mu2=", nsim, ";var2=", nsim),
      distribution = "mu~runif; var~runif(0, mu*(1-mu))",
      summary_policy = "pv_audit_review2_summary means and mu quantiles",
      nsim = nsim
    ),
    "EC-033",
    "stage8rnq_patch_EC033_stochastic_seed_metadata_r.csv",
    "stochastic_seed_metadata",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates,
    notes = "Stage 8R-NQ EC-033 stochastic/seed metadata sidecar"
  )
  pv_audit_write_review2_ec034_julia_layers(
    pqdf,
    mudf,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_audit_intermediates
  )
  p_muvar <- (ggplot(mudf, aes(mu, var)) + geom_hex(bins = 120) + geom_smooth(se=FALSE, color="red") + theme_bw() + labs(x="mu", y="sigma^2")) |
    (ggplot(mudf, aes(mu, delta)) + geom_hex(bins = 120) + geom_smooth(se=FALSE, color="red") + theme_bw() + labs(x="mu", y="delta"))
  save_png(p_muvar, "review2_2.1muvarsolution_R_equivalent.png", width = 9, height = 4)

  list(problem = p_problem, solution = p_solution, problem_full = p_problem_full, solution_full = p_solution_full, pq = pqdf, muvar = mudf)
}

# ============================================================ #
# 02O OPTIONAL R REIMPLEMENTATION RUNNER ####
# ============================================================ #
run_all_R_replication <- function(root = repo_root) {
  check_repo_root(root)
  audit_paths <- pv_audit_paths_if_enabled(create = TRUE)
  old_output_dirs <- NULL
  if (!is.null(audit_paths)) {
    old_output_dirs <- activate_pv_audit_output_dirs(audit_paths)
    on.exit(restore_pv_audit_output_dirs(old_output_dirs), add = TRUE)
    write_audit_session_info(audit_paths = audit_paths)
  }

  message("Pellett & Valbuena input folder: ", root)
  message("Optional output folder: ", out_dir)

  message("Running the selected Figure 1 concept procedures in R...")
  fig1 <- run_figure1_concept_R(audit_paths = audit_paths)

  message("Running the selected Figure 2 mean-biased heterogeneity procedures in R...")
  fig2 <- run_figure2_MBH_comb_R(root, audit_paths = audit_paths)

  message("Running tests of mean-biased heterogeneity measures...")
  mbh_tests <- run_mean_biased_hypothesis_tests_R(root, audit_paths = audit_paths)

  message("Preparing the Catalonia data and fitting the selected Allouche models...")
  adf <- prepare_allouche_data_R(root, audit_paths = audit_paths)
  allouche <- run_allouche_analysis_R(adf, audit_paths = audit_paths)
  table1_rescaled <- write_table1_rescaled_model_terms(allouche$allouche_model_df)

  message("Running the selected Figure 3 corrected-HDR procedures in R...")
  fig3 <- run_figure3_corrected_HDR_R(allouche, root, audit_paths = audit_paths)

  message("Reimplementing selected supplementary-material Figure S1 analyses in R...")
  s1 <- run_som_s1_fixedCV_R(audit_paths = audit_paths)

  message("Reimplementing selected supplementary-material Figures S2–S5 analyses in R...")
  s2s5 <- run_som_s2_s5_R(audit_paths = audit_paths)

  message("Reimplementing selected supplementary-material Figure S6 analyses in R...")
  s6 <- run_som_s6_delta_location_R(adf, audit_paths = audit_paths)

  message("Running MacArthur foliage-height reanalysis...")
  mac <- run_macarthur_analysis_R(root, audit_paths = audit_paths)

  message("Running tests for negligible trends in δ...")
  delta_trend <- run_delta_trend_test_R(root, audit_paths = audit_paths)

  message("Running the selected review-response simulations and figures...")
  review2 <- run_review2_figures_R(audit_paths = audit_paths)

  if (RUN_EXPENSIVE_GLOBAL_RESAMPLING) {
    message("RUN_EXPENSIVE_GLOBAL_RESAMPLING is TRUE, but these partial R ports do not complete global sampling.")
  }

  manifest <- tibble(
    component = c("Figure1", "Figure2", "MBH hypothesis tests", "Allouche analysis", "Figure3", "SOM S1", "SOM S2-S5", "SOM S6", "MacArthur", "delta trend", "review2"),
    status = "completed_R_port",
    note = c(
      "R figure equivalent, not pixel-perfect CairoMakie",
      "R figure equivalent from analytic and empirical data",
      "Direct statistical port of empirical_hypo_test.jl",
      "Direct statistical port of prepare_data.jl and analysis.jl",
      "R figure equivalent using same model logic",
      "R figure equivalent",
      "R figure equivalents and polynomial summaries",
      "Direct scatter-plot port",
      "Direct statistical port; entropy estimator checked against Julia run",
      "Direct statistical port of 1-fig2t.jl",
      "R figure equivalents; random simulations use same n by default"
    )
  )
  write_csv(manifest, safe_file(out_dir, "R_replication_manifest.csv"))

  message("Optional tables written to: ", tab_dir)
  message("Optional figures written to: ", fig_dir)
  list(
    manifest = manifest, figure1 = fig1, figure2 = fig2, mbh_tests = mbh_tests,
    allouche = allouche, table1_rescaled = table1_rescaled, figure3 = fig3, s1 = s1, s2s5 = s2s5, s6 = s6,
    macarthur = mac, delta_trend = delta_trend, review2 = review2
  )
}
