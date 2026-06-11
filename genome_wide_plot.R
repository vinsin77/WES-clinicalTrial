library(dplyr)
library(ggplot2)

scores <- read.table("scores.gistic",
                     header    = TRUE,
                     sep       = "\t",
                     stringsAsFactors = FALSE
)

colnames(scores) <- c("Type", "Chromosome", "Start", "End",
                      "log10q", "Gscore", "Amplitude", "Frequency")

scores <- scores %>%
    mutate(
        Chromosome = as.character(Chromosome),
        plot_score = ifelse(Type == "Amp", log10q, -log10q)
    ) %>%
    filter(!is.na(Chromosome))

chr_lengths <- data.frame(
    Chromosome = as.character(1:22),
    length = c(248956422, 242193529, 198295559, 190214555,
               181538259, 170805979, 159345973, 145138636,
               138394717, 133797422, 135086622, 133275309,
               114364328, 107043718, 101991189, 90338345,
               83257441,  80373285,  58617616,  64444167,
               46709983,  50818468)
) %>%
    mutate(
        offset = cumsum(as.numeric(lag(length, default = 0))),
        mid    = offset + length / 2
    )

scores <- scores %>%
    left_join(chr_lengths, by = "Chromosome") %>%
    mutate(
        pos_start = Start + offset,
        pos_end   = End   + offset,
        pos_mid   = (pos_start + pos_end) / 2,
        Chromosome = factor(Chromosome, levels = as.character(1:22))
    )
scores %>%
    filter(Type == "Amp") %>%
    group_by(Chromosome) %>%
    slice_max(log10q, n = 1) %>%
    ungroup() %>%
    mutate(Chromosome_char = as.character(Chromosome)) %>%
    left_join(cytobands_for_join, by = "Chromosome_char", relationship = "many-to-many") %>%
    filter(Start >= start, Start <= end) %>%
    names()
lesions <- read.table(
    "all_lesions.conf_95.txt",
    header    = TRUE,
    sep       = "\t",
    fill      = TRUE,
    quote     = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
lesions <- lesions[, nzchar(names(lesions)) & !is.na(names(lesions))]

sig_peaks <- lesions %>%
    filter(`q values` < 0.05) %>%
    mutate(Type = ifelse(grepl("Amplification", `Unique Name`), "Amp", "Del"))

chr_shade <- chr_lengths %>%
    mutate(shade = ifelse(as.integer(Chromosome) %% 2 == 0, "even", "odd"))

total_genome <- sum(chr_lengths$length)
min_width    <- total_genome * 0.003

scores <- scores %>%
    mutate(
        bar_start = pmax(0, pos_mid - pmax((pos_end - pos_start) / 2, min_width / 2)),
        bar_end   = pos_mid + pmax((pos_end - pos_start) / 2, min_width / 2)
    )

ylim_val <- max(abs(scores$log10q), na.rm = TRUE) * 1.1

# cytobands
cytobands <- read.table(
    "cytoBand.txt",
    header = FALSE, sep = "\t",
    col.names = c("chrom", "start", "end", "name", "gieStain")
) %>%
    filter(chrom %in% paste0("chr", 1:22)) %>%
    mutate(
        Chromosome = as.character(gsub("chr", "", chrom)),
        Chromosome = factor(Chromosome, levels = as.character(1:22))
    ) %>%
    left_join(chr_lengths, by = "Chromosome") %>%
    mutate(
        pos_start = start + offset,
        pos_end   = end   + offset,
        fill_color = case_when(
            gieStain == "gneg"    ~ "white",
            gieStain == "gpos25"  ~ "grey75",
            gieStain == "gpos50"  ~ "grey50",
            gieStain == "gpos75"  ~ "grey25",
            gieStain == "gpos100" ~ "black",
            gieStain == "acen"    ~ "#C06060",
            gieStain == "gvar"    ~ "grey70",
            gieStain == "stalk"   ~ "#C06060",
            TRUE                  ~ "white"
        )
    )

# Peak labels via proper join
cytobands_for_join <- cytobands %>%
    mutate(Chromosome_char = as.character(Chromosome))

peak_labels_amp <- scores %>%
    filter(Type == "Amp") %>%
    group_by(Chromosome) %>%
    slice_max(log10q, n = 1) %>%
    ungroup() %>%
    mutate(Chromosome_char = as.character(Chromosome)) %>%
    left_join(cytobands_for_join, by = "Chromosome_char", relationship = "many-to-many") %>%
    filter(Start >= start, Start <= end) %>%
    mutate(
        cytoband_name = paste0(Chromosome_char, name),
        label_y       = log10q + ylim_val * 0.03
    ) %>%
    select(Chromosome = Chromosome.x, pos_mid, log10q, cytoband_name, label_y) %>%
    distinct(Chromosome, .keep_all = TRUE)

peak_labels_del <- scores %>%
    filter(Type == "Del") %>%
    group_by(Chromosome) %>%
    slice_max(log10q, n = 1) %>%
    ungroup() %>%
    mutate(Chromosome_char = as.character(Chromosome)) %>%
    left_join(cytobands_for_join, by = "Chromosome_char", relationship = "many-to-many") %>%
    filter(Start >= start, Start <= end) %>%
    mutate(
        cytoband_name = paste0(Chromosome_char, name),
        label_y       = -(log10q + ylim_val * 0.03)
    ) %>%
    select(Chromosome = Chromosome.x, pos_mid, log10q, cytoband_name, label_y) %>%
    distinct(Chromosome, .keep_all = TRUE)

# Ideogram strip
ideogram_y_min <- -ylim_val * 1.15
ideogram_y_max <- -ylim_val * 1.05

p <- ggplot() +
    geom_rect(data = filter(chr_shade, shade == "even"),
              aes(xmin = offset, xmax = offset + length,
                  ymin = -ylim_val, ymax = ylim_val),
              fill = "gray80", alpha = 1, show.legend = FALSE) +
    
    geom_rect(data = filter(scores, Type == "Amp"),
              aes(xmin = bar_start, xmax = bar_end,
                  ymin = 0, ymax = log10q),
              fill = "#C0392B", alpha = 0.85) +
    
    geom_rect(data = filter(scores, Type == "Del"),
              aes(xmin = bar_start, xmax = bar_end,
                  ymin = -log10q, ymax = 0),
              fill = "#2471A3", alpha = 0.85) +
    
    geom_text(data = filter(peak_labels_amp, log10q > 1.301),
              aes(x = pos_mid, y = label_y, label = cytoband_name),
              angle = 90, hjust = 0, size = 2.5, color = "#C0392B") +
    
    geom_text(data = filter(peak_labels_del, log10q > 1.301),
              aes(x = pos_mid, y = label_y, label = cytoband_name),
              angle = 90, hjust = 1, size = 2.5, color = "#2471A3") +
    
    geom_rect(data = cytobands,
              aes(xmin = pos_start, xmax = pos_end,
                  ymin = ideogram_y_min, ymax = ideogram_y_max),
              fill = cytobands$fill_color, color = NA) +
    geom_rect(data = chr_lengths,
              aes(xmin = offset, xmax = offset + length,
                  ymin = ideogram_y_min, ymax = ideogram_y_max),
              fill = NA, color = "black", linewidth = 0.3) +
    
    geom_hline(yintercept = 0,      color = "black",  linewidth = 0.4) +
    geom_hline(yintercept =  1.301, color = "grey40", linewidth = 0.3, linetype = "dashed") +
    geom_hline(yintercept = -1.301, color = "grey40", linewidth = 0.3, linetype = "dashed") +
    
    scale_x_continuous(
        breaks = chr_lengths$mid,
        labels = as.character(1:22),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        name   = expression(-log[10](q)),
        breaks = seq(-60, 60, by = 20)
    ) +
    coord_cartesian(ylim = c(-ylim_val * 1.2, ylim_val * 1.2)) +
    labs(x = "Chromosome") +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x     = element_text(size = 8),
        panel.grid      = element_blank(),
        axis.line.x     = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    )
p
