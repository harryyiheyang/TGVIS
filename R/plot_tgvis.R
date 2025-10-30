#' Plot GWAS Results with Credible-Set Annotation
#'
#' This function visualizes GWAS results (Z-score, -log10(p), or cs.pratt)
#' with fine-mapping annotations derived from a credible-set summary table.
#' It supports SNP/Gene distinction via shapes and uses ggsci scientific color palettes.
#'
#' @param gwas_df A data frame containing at least: \code{SNP}, \code{Zscore}, \code{p}, \code{CHR}, and \code{BP}.
#' @param summary_df A data frame with fine-mapping summary, containing columns:
#' \code{variable}, \code{cs}, \code{cs.pip}, \code{cs.pratt}, \code{xqtl}, \code{type}.
#' @param cs.pratt_thres Numeric. Pratt index threshold for labeling and dashed line (default: 0.1).
#' @param y Character. One of \code{"z"}, \code{"p"}, or \code{"cs.pratt"} to determine y-axis variable.
#' @param palette Character. ggsci color palette name. Default is \code{"locuszoom"}.
#'
#' @return A ggplot object.
#'
#' @export
#' @importFrom dplyr mutate select distinct group_by summarise filter left_join case_when
#' @importFrom tidyr unnest
#' @importFrom stringr str_split str_trim
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot geom_point labs theme_bw theme guides guide_legend scale_shape_identity
#' @importFrom grid unit arrow
#' @importFrom stats na.omit
#' @import ggsci
#'
plot_tgvis <- function(
    gwas_df,
    summary_df,
    cs.pratt_thres = 0.1,
    y = c("z","p","cs.pratt"),
    palette = c(
      "locuszoom","npg","lancet","jco","nejm","d3",
      "ucscgb","aaas","igv","cosmic","uchicago","tron",
      "rickandmorty","futurama","simpsons"
    )
){
  y <- match.arg(y)
  palette <- match.arg(palette)

  if (!requireNamespace("ggsci", quietly = TRUE)) {
    stop('Package "ggsci" is required. Please install it: install.packages("ggsci")')
  }
  scale_fill_fun  <- get(paste0("scale_fill_",  palette), envir = asNamespace("ggsci"))
  scale_color_fun <- get(paste0("scale_color_", palette), envir = asNamespace("ggsci"))

  req_gwas <- c("SNP","Zscore","p","CHR","BP")
  miss_gwas <- setdiff(req_gwas, names(gwas_df))
  if (length(miss_gwas) > 0) stop("GWAS is missing columns: ", paste(miss_gwas, collapse = ", "))

  req_sum <- c("variable","cs","cs.pip","cs.pratt","xqtl","type")
  miss_sum <- setdiff(req_sum, names(summary_df))
  if (length(miss_sum) > 0) stop("summary is missing columns: ", paste(miss_sum, collapse = ", "))

  gwas <- dplyr::mutate(gwas_df, BP_Mb = BP / 1e6)

  map_df <- summary_df %>%
    dplyr::mutate(
      xqtl_raw = xqtl,
      xqtl = stringr::str_split(xqtl_raw, "~"),
      type_code = dplyr::case_when(
        tolower(type) == "snp"  ~ 21L,
        tolower(type) == "gene" ~ 22L,
        TRUE ~ NA_integer_
      )
    ) %>%
    tidyr::unnest(xqtl) %>%
    dplyr::mutate(xqtl = stringr::str_trim(xqtl)) %>%
    dplyr::select(SNP = xqtl, cs, variable, cs.pratt, type_code) %>%
    dplyr::distinct()

  snp_group_all <- map_df %>%
    dplyr::group_by(SNP) %>%
    dplyr::summarise(
      cs_all  = paste(sort(unique(cs)), collapse = ","),
      cs_primary = min(cs),
      type_code_primary = {
        mx <- suppressWarnings(max(type_code, na.rm = TRUE))
        if (is.infinite(mx)) NA_integer_ else as.integer(mx)
      },
      cspratt_max = {
        mx <- suppressWarnings(max(cs.pratt, na.rm = TRUE))
        if (is.infinite(mx)) NA_real_ else as.numeric(mx)
      },
      .groups = "drop"
    )

  lab_df <- map_df %>%
    dplyr::filter(cs.pratt >= cs.pratt_thres) %>%
    dplyr::group_by(SNP) %>%
    dplyr::summarise(label = paste(unique(variable), collapse = "\n"), .groups = "drop")

  gwas_anno <- gwas %>%
    dplyr::left_join(snp_group_all, by = "SNP") %>%
    dplyr::left_join(lab_df, by = "SNP") %>%
    dplyr::mutate(
      cs_primary = as.factor(cs_primary),
      yval = dplyr::case_when(
        y == "z" ~ Zscore,
        y == "p" ~ -log10(p),
        y == "cs.pratt" ~ cspratt_max
      )
    )

  y_lab <- switch(y,
                  "z" = "Z-score",
                  "p" = expression(-log[10](p)),
                  "cs.pratt" = "Pratt Index of Credible Set"
  )

  present_shapes <- sort(unique(stats::na.omit(gwas_anno$type_code_primary)))
  shape_breaks   <- present_shapes
  shape_labels   <- c("SNP","Gene")[match(shape_breaks, c(21, 22))]

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = gwas_anno,
      ggplot2::aes(x = BP_Mb, y = yval),
      color = "grey80", size = 3
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(gwas_anno, !is.na(cs_primary), !is.na(type_code_primary)),
      ggplot2::aes(x = BP_Mb, y = yval, fill = cs_primary, color = cs_primary, shape = type_code_primary),
      size = 5, stroke = 0.2
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(gwas_anno, !is.na(label)),
      ggplot2::aes(x = BP_Mb, y = yval, label = label),
      size = 3,
      min.segment.length = 0.05,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.size = 0.3,
      segment.color = "grey30",
      arrow = grid::arrow(length = grid::unit(0.01, "npc")),
      max.overlaps = Inf,
      seed = 1
    ) +
    scale_fill_fun(name = "Credible set") +
    scale_color_fun(name = "Credible set") +
    ggplot2::labs(
      x = paste0("BP in MB ", "(CHR=", unique(gwas_df$CHR), ")"),
      y = y_lab,
      title = if (y == "cs.pratt") "Locus Plot of Fine-Mapping" else "Locus Zoom Plot of GWAS"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12),
      legend.spacing.x = grid::unit(0.3, "cm"),
      legend.justification = "center",
      legend.box.just = "center",
      legend.box.spacing = grid::unit(0, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(nrow = 1, override.aes = list(shape = 21, size = 3)),
      color = "none"
    )

  if (length(shape_breaks) > 0) {
    p <- p +
      ggplot2::scale_shape_identity(
        guide  = "legend",
        name   = "Type",
        breaks = shape_breaks,
        labels = shape_labels
      ) +
      ggplot2::guides(
        shape = ggplot2::guide_legend(nrow = 1, override.aes = list(fill = "grey60", color = "black", size = 3))
      )
  } else {
    p <- p + ggplot2::guides(shape = "none")
  }

  if (y == "cs.pratt") {
    p <- p +
      ggplot2::geom_hline(yintercept = cs.pratt_thres, linetype = "dashed", size = 0.4, color = "grey40") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2))
  }

  return(p)
}
