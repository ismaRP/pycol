#!/usr/bin/env Rscript

library(readr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(ggplot2, warn.conflicts = F)
library(magrittr, warn.conflicts = F)

library(getopt, warn.conflicts = F)
library(MALDIutils, warn.conflicts = F)
library(ggpubr, warn.conflicts = F)

options(warn=-1)

optspec = matrix(c(
  'peaks',      'p', 1, 'character',
  'markers',    'm', 1, 'character',
  'metadata',   'd', 1, 'character',
  'min_frac',   'r', 1, 'double',
  'outfolder',  'o', 1, 'character',
  'outmarkers', 'f', 1, 'character'
), byrow=TRUE, ncol=4)
opt = getopt(optspec)

if (is.null(opt$min_frac)) opt$min_frac=0.1


peaks_file = opt$peaks
markers_file = opt$markers
metadata_file = opt$metadata
min_frac = opt$min_frac

metadata = read_csv(metadata_file, show_col_types = F)

markers = read_csv(markers_file, show_col_types = F) %>% # filter(multiseq==FALSE)
  arrange(mass1) %>%
  mutate(missed.cleaves = as.integer(missed.cleaves))


peaks = readRDS(peaks_file)
peaks = peaks[paste0(metadata$sample_name, '_', metadata$replicate)]

mask = metadata$taxid != -1
metadata = metadata[mask, ]
peaks = peaks[mask]

nmarkers = nrow(markers)
npeaks = length(peaks)
message(
  sprintf('Matching %d markers against %d peaks spectra objects ...',
          nmarkers, npeaks))
message('Calculating fraction of spectra that each marker matches ...')
markers = pept_fly(markers, peaks, peaksby = metadata$taxid, tolerance=0.002,
                   match_tol = 300, match_tol_unit = 'ppm', aug_deam=T)
message('DONE\n')


message('Making plots...')
message('\tMass vs error (ppm)')
mass_error_plot = ggplot(markers) +
  geom_point(aes(x=mass1, y=error_ppm, color=multiseq), size=2, alpha=0.5) +
  geom_smooth(aes(x=mass1, y=error_ppm)) +
  xlab('m/z') + ylab('Error ppm') +
  scale_color_discrete('name' = 'Other peptides within 0.1 Da') +
  geom_hline(aes(yintercept=0)) +
  theme_pubclean() +
  theme(axis.line.x.bottom = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0.5), 'cm'),
        panel.grid.major.x = element_line(linetype = 'dotted', linewidth=0.5, color='grey'),
        legend.title=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), axis.text=element_text(size=8)
  )

message('\tFraction of appearance vs error (ppm')
frac_error_plot = ggplot(markers) +
  geom_point(aes(x=malditof_frac, y=error_ppm, color=multiseq), size=2, alpha=0.5) +
  geom_vline(aes(xintercept=x), color='red', data=data.frame(x=min_frac)) +
  xlab('Fraction of samples') +
  scale_color_discrete('name' = 'Other peptides within 0.1 Da') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(position = 'right', sec.axis = dup_axis()) +
  theme_pubclean() +
  theme(axis.line.x.bottom = element_blank(),
        axis.line.y.left = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
        panel.grid.major.x = element_line(linetype = 'dotted', linewidth=0.5, color='grey'),
        legend.title=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), axis.text=element_text(size=8)
  )

message('\tCombining and saving ...')
p = ggarrange(mass_error_plot, frac_error_plot, ncol=2, nrow=1, common.legend = T)
ggsave('comb_error.pdf', p, 'pdf', path=opt$outfolder, width=7.5, height=3.5)
message('\tDONE\n')

message('\tFraction by multisequence cluster')
frac_mc_plot = markers %>% filter(malditof_frac>0.1) %>%
  ggplot() +
  geom_jitter(aes(x=multiseq, y=malditof_frac, color=multiseq), size=3, alpha=0.8) +
  xlab('Multiseq') + ylab('Fraction of samples') +
  scale_color_discrete('name' = 'Multiple sequences\n with same mass')
message('\tSaving...')
ggsave('frac_mc_plot.pdf', frac_mc_plot, 'pdf', path=opt$outfolder, width=8, height=5)
message('\tDONE\n')

message('\tError (ppm) by multisequence cluster')
error_mc_plot = markers %>% filter(malditof_frac>0.1) %>%
  ggplot() +
  geom_jitter(aes(x=multiseq, y=error_ppm, color=multiseq), size=3, alpha=0.8) +
  xlab('Multiseq') + ylab('Error ppm') +
  scale_color_discrete('name' = 'Has other peptides\n within 0.1 Da')
message('\tSaving...')
ggsave('error_mc_plot.pdf', error_mc_plot, 'pdf', path=opt$outfolder, width=8, height=5)
message('\tDONE\n')

message(sprintf('Filtering markers that have a match in more than %f of the spectra ...', min_frac))
markers_zooms = markers %>%
  filter(malditof_maxfrac > min_frac)
nmarkers_filt = nrow(markers_zooms)
message(sprintf('\t%d markers left', nmarkers_filt))
message('DONE\n')

message(sprintf('Saving markers in %s ...', opt$outmarkers))
write_csv(markers_zooms, opt$outmarkers)
message('DONE\n')
