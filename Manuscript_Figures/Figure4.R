suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

infile <- "../Data/Fig4_Data.txt"     # <- your per-replicate data file

df <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("replicate","tool","empirical_FDR","recall") %in% names(df)))
df$empirical_FDR <- as.numeric(df$empirical_FDR)
df$recall        <- as.numeric(df$recall)

nominal <- 0.05

## ---- 2. classify each row: individual tool vs consensus vote rule ----------
individual_tools <- c("edgeR","DESeq2","ALDEx2","metaSeq","ADAPT","ANCOMBC2","MaAsLin3")

classify <- function(tool) {
  if (tool %in% individual_tools) return("Individual tool")
  if (tool %in% c("Union","Intersect") || grepl("^Vote_", tool)) return("Consensus rule")
  NA_character_
}
df$type <- vapply(df$tool, classify, character(1))

## Harmonise vote labels across 6- and 7-tool replicates onto a common
## consensus ordering. Map every consensus rule to its k-of-N fraction, then
## to a shared display label on the 7-tool scale (nearest rung).
consensus_frac <- function(tool) {
  if (tool == "Union")     return(1/7)          # union = at least 1
  if (tool == "Intersect") return(7/7)          # intersect = all
  m <- regmatches(tool, regexec("^Vote_(\\d+)of(\\d+)$", tool))[[1]]
  if (length(m) == 3) return(as.numeric(m[2]) / as.numeric(m[3]))
  NA_real_
}
df$frac <- ifelse(df$type == "Consensus rule", vapply(df$tool, consensus_frac, numeric(1)), NA_real_)

## shared consensus display label (order from Union at top to Intersect)
frac_to_label <- function(f) {
  if (is.na(f)) return(NA_character_)
  # nearest k on a 7 scale
  k <- round(f * 7)
  if (k <= 1) return("Union (1 of 7)")
  if (k >= 7) return("Intersect (7 of 7)")
  sprintf("%d of 7", k)
}
df$consensus_label <- vapply(df$frac, frac_to_label, character(1))

## unified label used on the y-axis of panels b/c: tool name for individual,
## consensus_label for consensus rules
df$label <- ifelse(df$type == "Individual tool", df$tool, df$consensus_label)

## ---- 3. per-label summary across replicates --------------------------------
summ <- df %>%
  group_by(type, label) %>%
  summarise(
    mean_FDR    = mean(empirical_FDR, na.rm = TRUE),
    sd_FDR      = sd(empirical_FDR,   na.rm = TRUE),
    mean_recall = mean(recall,        na.rm = TRUE),
    sd_recall   = sd(recall,          na.rm = TRUE),
    pct_controlled = 100 * mean(empirical_FDR <= nominal, na.rm = TRUE),
    n_rep       = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(sd_FDR = ifelse(is.na(sd_FDR), 0, sd_FDR),
         sd_recall = ifelse(is.na(sd_recall), 0, sd_recall))

## y-axis order: consensus rules Union->Intersect, then individual tools
consensus_order <- c("Union (1 of 7)","2 of 7","3 of 7","4 of 7","5 of 7","6 of 7","Intersect (7 of 7)")
tool_order      <- individual_tools
lvl <- c(rev(tool_order), rev(consensus_order))   # bottom -> top for ggplot
summ$label <- factor(summ$label, levels = lvl)

## colours: orange = consensus, blue = individual (matching your reference)
cols <- c("Consensus rule" = "#E8890C", "Individual tool" = "#2b7fff")

## ---- panel a: FDR vs power trade-off ---------------------------------------
pa <- ggplot(summ, aes(mean_FDR, mean_recall, colour = type)) +
  geom_vline(xintercept = nominal, linetype = "dotted", colour = "grey50") +
  geom_path(data = dplyr::filter(summ, type == "Consensus rule") %>%
              dplyr::arrange(match(label, rev(lvl))),
            aes(group = 1), colour = "#E8890C", linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.6, max.overlaps = Inf,
                           show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "Empirical FDR (nominal 0.05)", y = "Power (recall)",
       title = "a  FDR / power trade-off",
       subtitle = "orange = consensus rules in vote order; blue = individual tools") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

## ---- panel b: FDR control at nominal q < 0.05 ------------------------------
pb <- ggplot(summ, aes(mean_FDR, label, colour = type)) +
  geom_vline(xintercept = nominal, linetype = "dotted", colour = "grey50") +
  geom_errorbarh(aes(xmin = pmax(0, mean_FDR - sd_FDR), xmax = mean_FDR + sd_FDR),
                 height = 0.3, linewidth = 0.4) +
  geom_point(size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "Empirical FDR (mean +/- SD)", y = NULL,
       title = "b  FDR control at nominal q < 0.05") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

## ---- panel c: reliability of control ---------------------------------------
pc <- ggplot(summ, aes(pct_controlled, label, colour = type)) +
  geom_vline(xintercept = 95, linetype = "dotted", colour = "grey50") +
  geom_point(size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(x = "% of replicates with FDR <= nominal", y = NULL,
       title = "c  Reliability of control",
       subtitle = "dotted line = 95% of replicates") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

## ---- assemble --------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- pa / (pb | pc) + patchwork::plot_layout(heights = c(1, 1))
  ggsave("FDR_summary_figure.pdf", fig, width = 10, height = 10, device = cairo_pdf)
} else {
  ggsave("FDR_summary_a.pdf", pa, width = 7, height = 5, device = cairo_pdf)
}

write.table(summ, "FDR_summary_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote FDR_summary_figure.pdf and FDR_summary_table.tsv")
print(summ)

  
