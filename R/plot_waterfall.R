#' Plot ClockSHAP waterfall for a single sample
#'
#' Create a horizontal waterfall plot showing age-adjusted ClockSHAP
#' contributions for a single sample, highlighting the top contributing
#' features and aggregating the remainder as "Other".
#'
#' The plot visualizes how feature-level contributions bridge the expected
#' age (conditioned on chronological age) to the predicted biological age.
#'
#' @param x A `clockshap` object.
#' @param sample A sample index (integer) or row name (character).
#' @param top_n Number of top absolute contributions to display.
#' @param show_age_labels Logical; show the expected/predicted age text labels.
#' @param show_delta_label Logical; show the Delta Age label above the arrow.
#' @param show_effect_labels Logical; show the contribution values inside bars.
#'
#' @return A ggplot object.
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
    show_effect_labels = TRUE
) {

  ## ------------------------------------------------------------
  ## Validation
  ## ------------------------------------------------------------
  if (!inherits(x, "clockshap")) {
    stop("`x` must be a clockshap object.", call. = FALSE)
  }

  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  phi <- clockshap_phi(x)

  ## resolve sample index
  if (is.character(sample)) {
    if (is.null(rownames(phi)) || !sample %in% rownames(phi)) {
      stop("Sample name not found in ClockSHAP object.", call. = FALSE)
    }
    i <- match(sample, rownames(phi))
  } else {
    i <- as.integer(sample)
    if (is.na(i) || i < 1 || i > nrow(phi)) {
      stop("Sample index out of bounds.", call. = FALSE)
    }
  }

  ## ------------------------------------------------------------
  ## Extract and order effects
  ## ------------------------------------------------------------
  phi_all <- phi[i, ]
  phi_all <- phi_all[order(-abs(phi_all))]

  top_eff   <- head(phi_all, top_n)
  other_eff <- sum(phi_all) - sum(top_eff)

  effects <- c(top_eff, Other = other_eff)

  ## ------------------------------------------------------------
  ## Expected and predicted ages
  ## ------------------------------------------------------------
  exp_i  <- clockshap_expected(x)[i]
  pred_i <- clockshap_predicted(x)[i]

  exp_i_r  <- round(exp_i, 1)
  pred_i_r <- round(pred_i, 1)
  delta_i  <- pred_i_r - exp_i_r

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

  ## vertical positioning
  y_arrow <- max(wf_plot$y) + 1.2
  y_label <- y_arrow + 0.8

  ## x-range padding
  x_min <- min(c(wf_plot$Start, wf_plot$End, exp_i, pred_i), na.rm = TRUE)
  x_max <- max(c(wf_plot$Start, wf_plot$End, exp_i, pred_i), na.rm = TRUE)

  x_pad_left  <- 0.04 * (x_max - x_min)
  x_pad_right <- 0.08 * (x_max - x_min)

  x_breaks <- pretty_breaks(n = 6)(
    c(x_min - x_pad_left, x_max + x_pad_right)
  )

  ## ------------------------------------------------------------
  ## Plot (base layers)
  ## ------------------------------------------------------------
  p <- ggplot(wf_plot) +

    ## custom vertical grid lines
    geom_segment(
      data = data.frame(x = x_breaks),
      aes(
        x = x, xend = x,
        y = min(wf_plot$y) - 0.5,
        yend = y_arrow
      ),
      inherit.aes = FALSE,
      colour = "grey92",
      size = 0.4
    ) +

    ## expected age line
    geom_segment(
      aes(
        x = exp_i, xend = exp_i,
        y = min(wf_plot$y) - 0.5,
        yend = y_arrow - 0.25
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      colour = "grey50"
    ) +

    ## predicted age line
    geom_segment(
      aes(
        x = pred_i, xend = pred_i,
        y = min(wf_plot$y) - 0.5,
        yend = y_arrow - 0.25
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      colour = "grey50"
    ) +

    ## arrow connecting expected -> predicted
    geom_segment(
      aes(x = exp_i, xend = pred_i, y = y_arrow, yend = y_arrow),
      inherit.aes = FALSE,
      colour = "grey30",
      size = 0.9,
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
      size = 8.5,
      lineend = "butt",
      colour = "grey30"
    ) +
    geom_segment(
      data = wf_plot,
      aes(x = Start, xend = End, y = y, yend = y, colour = fill_key),
      size = 7,
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
      limits = c(x_min - x_pad_left, x_max + x_pad_right),
      breaks = x_breaks
    ) +

    labs(x = "Age (years)", y = NULL) +
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
      label = as.expression(
        bquote(Delta * "Age = " * .(sprintf("%+.1f", delta_i)) * " yr")
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
        label = sprintf("Expected Age\n%.1f yr", exp_i_r),
        hjust = 0.5, vjust = 0,
        size = 4,
        colour = "grey20"
      ) +
      annotate(
        "text",
        x = pred_i,
        y = y_label,
        label = sprintf("Predicted Age\n%.1f yr", pred_i_r),
        hjust = 0.5, vjust = 0,
        size = 4,
        colour = "grey20"
      )
  }

  if (isTRUE(show_effect_labels)) {
    p <- p + geom_text(
      data = wf_plot,
      aes(x = Mid, y = y, label = sprintf("%+.1f", Eff)),
      colour = "white",
      size = 3.6,
      fontface = "bold"
    )
  }

  return(p)
}
