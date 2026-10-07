#' Plot ClockSHAP waterfall for a single sample
#'
#' Create a horizontal waterfall plot showing age-adjusted ClockSHAP
#' contributions for a single sample, highlighting the top contributing
#' features and aggregating the remainder as "Other".
#'
#' The plot visualizes how feature-level contributions bridge the expected
#' age (conditioned on chronological age) to the predicted age.
#'
#' @param x A `clockshap` object.
#' @param sample A sample index (integer) or row name (character).
#' @param top_n Number of top absolute contributions to display.
#' @param show_age_labels Logical; show the expected/predicted age text labels.
#' @param show_delta_label Logical; show the deviation label above the arrow.
#' @param show_effect_labels Logical; show the contribution values inside bars.
#' @param effect_label_mode Label placement strategy for effect labels. `"auto"`
#'   places labels inside sufficiently large bars and outside small bars.
#'   `"inside"` places all labels inside, `"outside"` places all labels outside,
#'   and `"hide_small"` hides labels for small bars.
#' @param effect_label_min_frac Numeric scalar giving the minimum bar width
#'   (as a fraction of x-axis span) required for inside placement.
#' @param effect_label_outside_nudge_frac Numeric scalar giving outside-label
#'   horizontal offset (as a fraction of x-axis span).
#' @param show_effect_label_connectors Logical; draw connector segments from bar
#'   ends to outside labels.
#' @param x_axis_title Character scalar for the x-axis title.
#' @param y_axis_title Character scalar for the y-axis title. Use `NULL` for no
#'   title.
#' @param deviation_label Text template for the deviation annotation. Use
#'   `{value}` as a placeholder for the numeric value.
#' @param expected_age_label Text template for the expected age annotation. Use
#'   `{value}` as a placeholder for the numeric value.
#' @param predicted_age_label Text template for the predicted age annotation.
#'   Use `{value}` as a placeholder for the numeric value.
#' @param effect_label Text template for effect labels shown inside bars. Use
#'   `{value}` as a placeholder for the numeric value.
#' @param digits Number of decimal places used for all numeric labels.
#' @param top_label_digits Deprecated compatibility argument. Use `digits`.
#' @param effect_label_digits Deprecated compatibility argument. Use `digits`.
#'
#' @return A ggplot object.
#'
#' @details Labels use reconciled rounding from
#'   [clockshap_waterfall_values()], so displayed contributions sum to the
#'   displayed deviation and displayed expected age plus deviation equals
#'   displayed predicted age. To preserve these identities, an endpoint label
#'   may differ from independently rounding its stored value by one unit in the
#'   last displayed decimal place.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_text annotate
#' @importFrom ggplot2 scale_y_continuous scale_x_continuous scale_color_manual
#' @importFrom ggplot2 coord_cartesian theme_minimal theme labs
#' @importFrom ggplot2 element_text element_blank margin expansion
#' @importFrom dplyr arrange mutate desc
#' @importFrom tibble tibble
#' @importFrom scales pretty_breaks
#' @importFrom magrittr %>%
#' @importFrom grid unit arrow
#' @importFrom utils head
#'
#' @export
plot_clockshap_waterfall <- function(
    x,
    sample,
    top_n = 10,
    show_age_labels = TRUE,
    show_delta_label = TRUE,
    show_effect_labels = TRUE,
    effect_label_mode = "auto",
    effect_label_min_frac = 0.03,
    effect_label_outside_nudge_frac = 0.02,
    show_effect_label_connectors = TRUE,
    x_axis_title = "Age (years)",
    y_axis_title = NULL,
    deviation_label = "Deviation = {value} yr",
    expected_age_label = "Expected Age\n{value} yr",
    predicted_age_label = "Predicted Age\n{value} yr",
    effect_label = "{value}",
    top_label_digits = NULL,
    effect_label_digits = NULL,
    digits = 1
) {

  ## ------------------------------------------------------------
  ## Validation
  ## ------------------------------------------------------------
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }

  top_n <- .validate_waterfall_top_n(top_n)
  effect_label_mode <- match.arg(
    effect_label_mode,
    c("auto", "inside", "outside", "hide_small")
  )
  if (!is.numeric(effect_label_min_frac) || length(effect_label_min_frac) != 1 ||
      !is.finite(effect_label_min_frac) || effect_label_min_frac < 0) {
    stop("`effect_label_min_frac` must be a single non-negative finite number.",
         call. = FALSE)
  }
  if (!is.numeric(effect_label_outside_nudge_frac) ||
      length(effect_label_outside_nudge_frac) != 1 ||
      !is.finite(effect_label_outside_nudge_frac) ||
      effect_label_outside_nudge_frac <= 0) {
    stop("`effect_label_outside_nudge_frac` must be a single positive finite number.",
         call. = FALSE)
  }
  if (!is.logical(show_effect_label_connectors) ||
      length(show_effect_label_connectors) != 1 ||
      is.na(show_effect_label_connectors)) {
    stop("`show_effect_label_connectors` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(x_axis_title) || length(x_axis_title) != 1) {
    stop("`x_axis_title` must be a single character string.", call. = FALSE)
  }
  if (!is.null(y_axis_title) &&
      (!is.character(y_axis_title) || length(y_axis_title) != 1)) {
    stop("`y_axis_title` must be NULL or a single character string.",
         call. = FALSE)
  }
  if (!is.character(deviation_label) || length(deviation_label) != 1) {
    stop("`deviation_label` must be a single character string.",
         call. = FALSE)
  }
  if (!is.character(expected_age_label) ||
      length(expected_age_label) != 1) {
    stop("`expected_age_label` must be a single character string.",
         call. = FALSE)
  }
  if (!is.character(predicted_age_label) ||
      length(predicted_age_label) != 1) {
    stop("`predicted_age_label` must be a single character string.",
         call. = FALSE)
  }
  if (!is.character(effect_label) || length(effect_label) != 1) {
    stop("`effect_label` must be a single character string.",
         call. = FALSE)
  }
  digits_missing <- missing(digits)
  if (!digits_missing) {
    digits <- .validate_waterfall_digits(digits)
  }
  legacy_digits <- c(
    top_label_digits = top_label_digits,
    effect_label_digits = effect_label_digits
  )
  legacy_digits <- legacy_digits[!vapply(legacy_digits, is.null, logical(1))]
  if (length(legacy_digits) > 0) {
    legacy_digits <- vapply(names(legacy_digits), function(argument) {
      .validate_waterfall_digits(legacy_digits[[argument]], argument)
    }, integer(1))
    warning(
      "`top_label_digits` and `effect_label_digits` are deprecated; use `digits`.",
      call. = FALSE
    )
    if (length(unique(legacy_digits)) > 1) {
      stop("Deprecated precision arguments must have the same value.",
           call. = FALSE)
    }
    if (!digits_missing && !identical(digits, legacy_digits[[1]])) {
      stop("`digits` conflicts with a deprecated precision argument.",
           call. = FALSE)
    }
    if (digits_missing) {
      digits <- legacy_digits[[1]]
    }
  }
  if (digits_missing && length(legacy_digits) == 0) {
    digits <- .validate_waterfall_digits(digits)
  }

  unsigned_format <- paste0("%.", digits, "f")

  expected_age_label <- gsub("\\\\n", "\n", expected_age_label)
  predicted_age_label <- gsub("\\\\n", "\n", predicted_age_label)

  phi <- clockshap_phi(x)

  i <- .resolve_waterfall_sample(phi, sample)

  ## ------------------------------------------------------------
  ## Extract and order effects
  ## ------------------------------------------------------------
  effects <- .waterfall_contributions(phi[i, ], top_n)
  display <- clockshap_waterfall_values(x, sample, top_n, digits)

  ## ------------------------------------------------------------
  ## Expected and predicted ages
  ## ------------------------------------------------------------
  exp_i  <- clockshap_expected(x)[i]
  pred_i <- clockshap_predicted(x)[i]
  if (!is.finite(exp_i) || !is.finite(pred_i)) {
    stop("Selected sample has non-finite expected or predicted age; cannot plot.",
         call. = FALSE)
  }

  ## ------------------------------------------------------------
  ## Cumulative positions for waterfall
  ## ------------------------------------------------------------
  start_vals <- c(exp_i, exp_i + cumsum(head(effects, -1)))
  end_vals   <- start_vals + effects

  wf <- tibble(
    Feature = names(effects),
    Start   = start_vals,
    End     = end_vals,
    Mid     = (start_vals + end_vals) / 2,
    Eff     = effects,
    DisplayEff = display$contributions,
    Sign    = effects > 0,
    IsOther = names(effects) == "Other"
  )

  ## ------------------------------------------------------------
  ## Arrange for plotting
  ## ------------------------------------------------------------
  wf_plot <- wf %>%
    arrange(Feature == "Other", desc(abs(Eff))) %>%
    mutate(
      y = rev(seq_len(nrow(wf))),
      fill_key = interaction(Sign, IsOther)
    )

  if (!isTRUE(all.equal(
    unname(wf_plot$Eff),
    unname(wf_plot$End - wf_plot$Start),
    tolerance = 1e-12
  ))) {
    stop("Waterfall invariant violated: bar effects do not match bar lengths.",
         call. = FALSE)
  }

  ## vertical positioning
  y_min <- min(wf_plot$y)
  y_max <- max(wf_plot$y)
  y_base <- y_min - 0.5
  y_arrow <- y_max + 1.2
  y_label <- y_arrow + 0.8

  ## x-range padding
  x_min <- min(c(wf_plot$Start, wf_plot$End, exp_i, pred_i), na.rm = TRUE)
  x_max <- max(c(wf_plot$Start, wf_plot$End, exp_i, pred_i), na.rm = TRUE)

  x_span <- x_max - x_min
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- 1
  }
  eff_values <- .format_waterfall_signed(wf_plot$DisplayEff, digits)
  wf_plot$eff_label <- vapply(
    eff_values,
    function(v) sub("{value}", v, effect_label, fixed = TRUE),
    character(1)
  )

  small_cutoff <- effect_label_min_frac * x_span
  is_small <- abs(wf_plot$Eff) < small_cutoff
  label_place <- rep("inside", nrow(wf_plot))
  if (effect_label_mode == "outside") {
    label_place[] <- "outside"
  } else if (effect_label_mode == "hide_small") {
    label_place[is_small] <- "hidden"
  } else if (effect_label_mode == "auto") {
    label_place[is_small] <- "outside"
  }
  wf_plot$label_place <- label_place

  outside_nudge <- effect_label_outside_nudge_frac * x_span
  connector_gap <- 0.005 * x_span
  wf_plot$label_x <- ifelse(
    wf_plot$Eff >= 0,
    wf_plot$End + outside_nudge,
    wf_plot$End - outside_nudge
  )
  wf_plot$label_hjust <- ifelse(wf_plot$Eff >= 0, 0, 1)
  wf_plot$connector_xend <- ifelse(
    wf_plot$Eff >= 0,
    wf_plot$label_x - connector_gap,
    wf_plot$label_x + connector_gap
  )

  x_pad_left  <- 0.04 * x_span
  x_pad_right <- 0.08 * x_span
  if (any(wf_plot$label_place == "outside")) {
    if (any(wf_plot$label_place == "outside" & wf_plot$Eff < 0)) {
      x_pad_left <- max(x_pad_left, outside_nudge + 0.14 * x_span)
    }
    if (any(wf_plot$label_place == "outside" & wf_plot$Eff >= 0)) {
      x_pad_right <- max(x_pad_right, outside_nudge + 0.14 * x_span)
    }
  }

  x_breaks <- pretty_breaks(n = 6)(
    c(x_min - x_pad_left, x_max + x_pad_right)
  )
  x_limits <- c(x_min - x_pad_left, x_max + x_pad_right)
  x_breaks <- x_breaks[x_breaks >= x_limits[1] & x_breaks <= x_limits[2]]

  ## ------------------------------------------------------------
  ## Plot (base layers)
  ## ------------------------------------------------------------
  p <- ggplot(wf_plot) +

    ## custom vertical grid lines
    geom_segment(
      data = data.frame(
        x = x_breaks,
        xend = x_breaks,
        y = y_base,
        yend = y_arrow
      ),
      aes(
        x = x, xend = xend,
        y = y, yend = yend
      ),
      inherit.aes = FALSE,
      colour = "grey92",
      linewidth = 0.4
    ) +

    ## expected age line
    annotate(
      "segment",
      x = exp_i, xend = exp_i,
      y = y_base, yend = y_arrow - 0.25,
      linetype = "dotted",
      colour = "grey50"
    ) +

    ## predicted age line
    annotate(
      "segment",
      x = pred_i, xend = pred_i,
      y = y_base, yend = y_arrow - 0.25,
      linetype = "dotted",
      colour = "grey50"
    ) +

    ## arrow connecting expected -> predicted
    annotate(
      "segment",
      x = exp_i, xend = pred_i, y = y_arrow, yend = y_arrow,
      colour = "grey30",
      linewidth = 0.9,
      arrow = arrow(
        ends = "both",
        type = "closed",
        length = unit(0.06, "inches")
      )
    ) +

    ## waterfall bars (two layers: grey base + colored overlay)
    geom_segment(
      data = wf_plot,
      aes(x = Start, xend = End, y = y, yend = y),
      linewidth = 8.5,
      lineend = "butt",
      colour = "grey30"
    ) +
    geom_segment(
      data = wf_plot,
      aes(x = Start, xend = End, y = y, yend = y, colour = fill_key),
      linewidth = 7,
      lineend = "butt"
    ) +

    ## colors
    scale_color_manual(
      values = c(
        "TRUE.FALSE"  = "#d73027",
        "FALSE.FALSE" = "#4575b4",
        "TRUE.TRUE"   = "#b2182b",
        "FALSE.TRUE"  = "#2166ac"
      ),
      guide = "none"
    ) +

    ## axes
    scale_y_continuous(
      breaks = wf_plot$y,
      labels = wf_plot$Feature,
      expand = expansion(mult = c(0.05, 0.18))
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks
    ) +

    labs(x = x_axis_title, y = y_axis_title) +
    coord_cartesian(clip = "off") +

    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10, face = "bold", margin = margin(r = 2)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(t = 10, r = 12, b = 8, l = 8)
    )

  ## ------------------------------------------------------------
  ## Optional annotation layers
  ## ------------------------------------------------------------

  if (isTRUE(show_delta_label)) {
    p <- p + annotate(
      "text",
      x = (exp_i + pred_i) / 2,
      y = y_label,
      label = sub(
        "{value}",
        .format_waterfall_signed(display$deviation, digits),
        deviation_label,
        fixed = TRUE
      ),
      size = 4.5,
      fontface = "bold",
      colour = "grey20"
    )
  }

  if (isTRUE(show_age_labels)) {
    p <- p +
      annotate(
        "text",
        x = exp_i,
        y = y_label,
        label = sub("{value}", sprintf(unsigned_format, display$expected), expected_age_label,
                    fixed = TRUE),
        hjust = 0.5, vjust = 0,
        size = 4,
        colour = "grey20"
      ) +
      annotate(
        "text",
        x = pred_i,
        y = y_label,
        label = sub("{value}", sprintf(unsigned_format, display$predicted), predicted_age_label,
                    fixed = TRUE),
        hjust = 0.5, vjust = 0,
        size = 4,
        colour = "grey20"
      )
  }

  if (isTRUE(show_effect_labels)) {
    wf_inside <- wf_plot[wf_plot$label_place == "inside", , drop = FALSE]
    wf_outside <- wf_plot[wf_plot$label_place == "outside", , drop = FALSE]

    if (nrow(wf_inside) > 0) {
      p <- p + geom_text(
        data = wf_inside,
        aes(x = Mid, y = y),
        label = wf_inside$eff_label,
        colour = "white",
        size = 3.6,
        fontface = "bold"
      )
    }

    if (nrow(wf_outside) > 0) {
      if (isTRUE(show_effect_label_connectors)) {
        p <- p + geom_segment(
          data = wf_outside,
          aes(x = End, xend = connector_xend, y = y, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.5,
          colour = "grey40"
        )
      }

      p <- p + geom_text(
        data = wf_outside,
        aes(x = label_x, y = y, hjust = label_hjust),
        inherit.aes = FALSE,
        label = wf_outside$eff_label,
        colour = "grey15",
        size = 3.6,
        fontface = "bold"
      )
    }
  }

  return(p)
}
