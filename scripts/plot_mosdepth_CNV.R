library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(stringr)
library(dplyr)
library(cowplot)

chrlist = read_csv("genome/chrom_nums.csv",
                   col_names = c("Chr", "CHR","Length"),
                   col_types = "cii")
# filter by length perhaps
chrlistshort = chrlist %>% filter(Length >= 100000)
mosdepthdir = "coverage/mosdepth"
num_per_plot = 8
windows = c(10000,50000)
alternating_colors = rep(c("red", "black"), times = length(chrlist$Chr))



plot_strain <- function(strain, data,repeats,genes) {
  l = data %>% filter(Strain == strain)
  r = repeats
  g = genes
  Title = sprintf("Read coverage plot for %s", strain)
  p <-
    ggplot(l, aes(
      x = pos,
      y = NormDepth,
      color = factor(CHR)
    )) +
    scale_colour_manual(values = alternating_colors) +
    geom_point(
      alpha = 0.9,
      size = 0.8,
      shape = 16,
      show.legend = FALSE
    ) +
    labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(l$CHR)
    ) +
    scale_y_continuous(name = "Normalized Read Depth",
                       expand = c(0, 0),
                       limits = c(0, 3)) +
    theme_classic()
  
  pR <-
    ggplot(l, aes(
      x = pos,
      y = NormDepth,
      color = factor(CHR)
    )) +
    scale_colour_manual(values = alternating_colors) +
    geom_line(
      alpha = 0.7,
      size = 0.5,
      show.legend = FALSE
    ) +
    labs(title = sprintf("Repeat density plot for %s", strain), xlab = "Position", y = "Repeat density") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(l$CHR)
    ) +
    scale_y_continuous(name = "Repeat density",
                       expand = c(0, 0),
                       limits = c(0, 1)) +
    theme_classic()
  
  pG <-
    ggplot(g, aes(
      x = pos,
      y = Covfrac,
      color = factor(CHR)
    )) +
    scale_colour_manual(values = alternating_colors) +
    geom_line(
      alpha = 0.7,
      size = 0.5,
      show.legend = FALSE
    ) +
    labs(title = sprintf("Gene density plot for %s", strain), xlab = "Position", y = "Gene density") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(l$CHR)
    ) +
    scale_y_continuous(name = "Gene density",
                       expand = c(0, 0),
                       limits = c(0, 1)) +
    theme_classic()
  comboPlot = plot_grid(p,pR,pG,ncol=1,align = "v")
  #+ guides(fill = guide_legend(keywidth = 3,keyheight = 1))
}

plot_chrs <- function(chrom, data,repeats, genes) {
  Title = sprintf("Chr%s depth of coverage", chrom)
  l = data %>% filter(CHR == chrom)
  r = repeats %>% filter(CHR == chrom)
  g = genes %>% filter(CHR == chrom)

  p <-
    ggplot(l, aes(
      x = Start,
      y = NormDepth,
      color = factor(Strain)
    )) +
    geom_point(
      alpha = 0.7,
      size = 0.8,
      shape = 16,
      show.legend = FALSE
    ) + # scale_color_brewer(palette='RdYlBu',type='seq') +
    labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
    scale_x_continuous(expand = c(0,	0), name = "Position") +
    scale_y_continuous(name = "Normalized Read Depth",
                       expand = c(0, 0),
                       limits = c(0, 5)) +
    theme_classic()

  pR <-
    ggplot(r, aes(
      x = Start,
      y = Covfrac,
    )) +
    geom_line(
      alpha = 0.7,
      size = 0.5,
      color="darkblue",
      show.legend = FALSE
    ) + # scale_color_brewer(palette='RdYlBu',type='seq') +
    labs(title = sprintf("Chr%s repeat density",chrom), xlab = "Position", y = "Repeat density") +
    scale_x_continuous(expand = c(0,	0), name = "Position") +
    scale_y_continuous(name = "Repeat density",
                       expand = c(0, 0),
                       limits = c(0, 1)) +
    theme_classic()
  
  pG <-
    ggplot(g, aes(
      x = Start,
      y = Covfrac,
    )) +
    geom_line(
      alpha = 0.7,
      size = 0.5,
      color="darkred",
      show.legend = FALSE
    ) + 
    labs(title = sprintf("Chr%s gene density",chrom), xlab = "Position", y = "Gene density") +
    scale_x_continuous(expand = c(0,	0), name = "Position") +
    scale_y_continuous(name = "Gene density",
                       expand = c(0, 0),
                       limits = c(0, 1)) +
    theme_classic()
  comboPlot = plot_grid(p,pR,pG,ncol=1,align = "v")
}

for (window in windows) {
  ChromDepths = tibble(strain = c(), ChromDepth = c())
  genewindows <- read_tsv(sprintf("features/gene_windows_%d.bed",window),
                          col_names=c("Chr","Start","End","Count","Covbases","Windowsize","Covfrac"))
  
  repeatwindows <- read_tsv(sprintf("features/repeat_windows_%d.bed",window),
                            col_names=c("Chr","Start","End","Count","Covbases","Windowsize","Covfrac"))
  
  inpattern = sprintf(".%sbp.regions.bed.gz$", window)
  file_list <- list.files(path = mosdepthdir, pattern = inpattern)
  bedwindows <- data.frame()
  for (i in 1:length(file_list)) {   # make me a function for lapply someday?
    infile = sprintf("%s/%s", mosdepthdir, file_list[i])
    strain = str_replace(file_list[i], inpattern, "")
    t = read_tsv(
      infile,
      col_names = c("Chr", "Start", "End", "Depth"),
      col_types = "ciid"
    )
    t$Strain = c(strain)
    medianDepth = median(t$Depth)
    t$NormDepth = t$Depth / medianDepth
    t = t %>% inner_join(., chrlist)
    ChromDepths <-
      bind_rows(
        ChromDepths,
        t %>% group_by(Strain, CHR) %>% summarize(meddepth = sprintf("%.2f",median(Depth)), .groups =
                                                    "keep")
      )
    bedwindows <- bind_rows(bedwindows, t)
  }
  length(bedwindows$Strain)
  unique(bedwindows$Strain)
  length(unique(bedwindows$Strain))
  unique(bedwindows$CHR)
  d = bedwindows %>% filter(CHR %in% chrlistshort$CHR) %>% arrange(CHR, Start, Strain)
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$Start, d$CHR, length))
  d$pos = NA
  rd = repeatwindows %>% inner_join(., chrlist) %>% filter(Chr %in% chrlistshort$Chr) %>% arrange(CHR, Start)
  rd$index = rep.int(seq_along(unique(rd$CHR)), times = tapply(rd$Start, rd$CHR, length))
  rd$pos = NA
  rg = genewindows  %>% inner_join(., chrlist) %>% filter(Chr %in% chrlistshort$Chr) %>% arrange(CHR, Start)
  rg$index = rep.int(seq_along(unique(rg$CHR)), times = tapply(rg$Start, rg$CHR, length))
  rg$pos = NA
  
  nchr = length(unique(chrlistshort$CHR))
  lastbase = 0
  ticks = NULL
  minor = vector(, 8)
  
  for (i in 1:nchr) {
    if (i == 1) {
      d[d$index == i, ]$pos = d[d$index == i, ]$Start
      rg[rg$index == i, ]$pos = rg[rg$index == i, ]$Start
      rd[rd$index == i, ]$pos = rd[rd$index == i, ]$Start
    } else {
      ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
      lastbase = lastbase + max(d[d$index == (i - 1), "Start"])
      minor[i] = lastbase
      
      d[d$index == i, "Start"] = d[d$index == i, "Start"] - min(d[d$index == i, "Start"]) + 1
      rd[rd$index == i, "Start"] = rd[rd$index == i, "Start"] - min(rd[rd$index == i, "Start"]) + 1
      rg[rg$index == i, "Start"] = rg[rg$index == i, "Start"] - min(rg[rg$index == i, "Start"]) + 1
      d[d$index == i, "End"] = lastbase
      rd[rd$index == i, "End"] = lastbase
      rg[rg$index == i, "End"] = lastbase
      
      d[d$index == i, "pos"] = d[d$index == i, "Start"] + lastbase
      rd[rd$index == i, "pos"] = rd[rd$index == i, "Start"] + lastbase
      rg[rg$index == i, "pos"] = rg[rg$index == i, "Start"] + lastbase
    }
  }
  ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
  # ticks
  minorB <- tapply(d$End, d$index, max, probs = 0.5)
  Title = "Depth of sequence coverage"

  p <- ggplot(d, aes(x = pos, y = NormDepth, color = Strain)) +
    geom_vline(
      mapping = NULL,
      xintercept = minorB,
      alpha = 0.5,
      size = 0.1,
      colour = "grey15"
    ) +
    geom_point(
      alpha = 0.8,
      size = 0.4,
      shape = 16,
      show.legend = FALSE,
    ) +
    labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(d$CHR)
    ) +
    scale_y_continuous(name = "Normalized Read Depth", limits = c(0, 5)) +
    theme_classic()
  
  pRP <- ggplot(rd, aes(x = pos, y = Covfrac)) +
    geom_vline(
      mapping = NULL,
      xintercept = minorB,
      alpha = 0.5,
      size = 0.1,
      colour = "grey15"
    ) +
    geom_line(
      alpha = 0.8,
      size = 0.5,
      show.legend = FALSE, color = "darkblue",
    ) +
    labs(title = "Repeat density", xlab = "Position", y = "Repeat density") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(rd$CHR)
    ) +
    scale_y_continuous(name = "Repeat density", limits = c(0, 1)) +
    theme_classic()
  
  pGene <- ggplot(rg, aes(x = pos, y = Covfrac)) +
    geom_vline(
      mapping = NULL,
      xintercept = minorB,
      alpha = 0.5,
      size = 0.1,
      colour = "grey15"
    ) +
    geom_line(
      alpha = 0.8,
      size = 0.5,
      show.legend = FALSE,
      color = "darkred",
    ) +
    labs(title = "Gene Density", xlab = "Position", y = "Gene density") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(rg$CHR)
    ) +
    scale_y_continuous(name = "Gene density", limits = c(0, 1)) +
    theme_classic()
  
  pdffile = sprintf("plots/genomewide_coverage_%dkb.pdf", window / 1000)
  ggsave(p,
         file = pdffile,
         width = 16,
         height = 4)
  comboPlot = plot_grid(p,pRP,pGene,ncol=1,align = "v")
  pdffile = sprintf("plots/genomewide_multidata_%dkb.pdf", window / 1000)
  ggsave(comboPlot,
         file = pdffile,
         width = 16,
         height = 4)
  
  f <- ggplot(repeatwindows, aes(x = pos, y = NormDepth, color = Strain)) +
    geom_vline(
      mapping = NULL,
      xintercept = minorB,
      alpha = 0.5,
      size = 0.1,
      colour = "grey15"
    ) +
    geom_point(
      alpha = 0.8,
      size = 0.4,
      shape = 16,
      show.legend = FALSE,
    ) +
    labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
    scale_x_continuous(
      name = "Chromosome",
      expand = c(0, 0),
      breaks = ticks,
      labels = unique(d$CHR)
    ) +
    scale_y_continuous(name = "Normalized Read Depth", limits = c(0, 5)) +
    theme_classic()
  # break these up
  countset = floor(length(file_list) / num_per_plot)
  if ( ! (countset * num_per_plot) == length(file_list) ) {
    countset = countset + 1
  }
  strainlist <- unique(d$Strain)

  for (n in 1:countset) {
    start=(n-1)*num_per_plot + 1
    end=n*num_per_plot
    end = min(length(file_list),end)

    dsplit = d %>% filter(Strain %in% strainlist[start:end])
    p <-
      ggplot(dsplit, aes(x = pos, y = NormDepth, color = Strain)) + geom_vline(
        mapping = NULL,
        xintercept = minorB,
        alpha = 0.5,
        size = 0.1,
        colour = "grey15"
      ) +
      geom_point(
        alpha = 0.8,
        size = 0.8,
        shape = 16,
        show.legend = TRUE
      ) +
      labs(title = Title,
           xlab = "Position",
           y = "Normalized Read Depth") +
      scale_x_continuous(
        name = "Chromosome",
        expand = c(0, 0),
        breaks = ticks,
        labels = unique(d$CHR)
      ) +
      scale_y_continuous(name = "Normalized Read Depth",
                         expand = c(0, 0),
                         limits = c(0, 5)) +
      scale_color_brewer(palette = 'Set1') +
      theme_classic() + guides(fill = guide_legend(keywidth = 3, keyheight = 1)) + theme(legend.position = "top")
    comb = plot_grid(p,pRP,pGene,ncol=1)
    pdffile = sprintf("plots/genomewide_coverage_%dkb_split%d.pdf",
                      window / 1000,
                      n)
    ggsave(comb,
           file = pdffile,
           width = 24,
           height =10)
  }
  
  plts <- lapply(unique(d$Strain), plot_strain, data = d, repeats=rd, genes=rg)
  
  strains = unique(d$Strain)
  for (i in 1:length(strains)) {
    pdffile = sprintf("plots/StrainPlot_%dkb.%s.pdf", window / 1000, strains[[i]])
    ggsave(plot = plts[[i]],
           file = pdffile,
           width = 16, height=2)
  }
  
  plts <- lapply(1:nchr, plot_chrs, data = d, repeats=rd,genes=rg)
  for (i in 1:nchr) {
    pdffile = sprintf("plots/ChrPlot_%dkb.Chr%s.pdf", window / 1000, i)
    ggsave(plot = plts[[i]], file = pdffile,width=20)
  }
  
  cd <-
    ChromDepths %>% pivot_wider(names_from = CHR,
                                values_from = meddepth,
                                names_prefix = "Chr")
  write_tsv(cd, sprintf("coverage/Chromosome_Depths_%dkb_windows.tsv",window / 1000))
}

