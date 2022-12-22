#!/usr/bin/env Rscript

library(readr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(magrittr, warn.conflicts = F)

library(getopt, warn.conflicts = F)
library(MALDIutils, warn.conflicts = F)

options(warn=-1)

optspec = matrix(c(
  'spectra',    's', 1, 'character',
  'metadata',   'd', 1, 'character',
  'markers',    'm', 1, 'character',
  'hws',        'h', 1, 'integer',
  'outfolder',  'o', 1, 'character',
  'outbase',    'b', 2, 'character',
  'ncores',     't', 2, 'integer',
  'iocores',    'j', 2, 'integer',
  'nchunks',    'c', 2, 'integer'
), byrow=TRUE, ncol=4)

opt = getopt(optspec)


if (is.null(opt$nchunks)) opt$nchunks=3
if (is.null(opt$ncores)) opt$ncores=6
if (is.null(opt$iocores)) opt$iocores=1

metadata_file = opt$metadata
markers_file = opt$markers
spectra_folder = opt$spectra

metadata = read_csv(metadata_file, show_col_types = F)
markers_zooms = read_csv(markers_file, show_col_types = F)

nspectra = nrow(metadata)
nmarkers = nrow(markers_zooms)

message(sprintf('Aligning %d markers to %d spectra ...', nmarkers, nspectra))
cor_data = get_max_cor(markers_zooms, spectra_folder, metadata,
                       opt$nchunks, opt$ncores, opt$iocores, vch=5,
                       myby = 0.01, gauss=0.2, laglim=0.6,
                       method='SavitzkyGolay', halfWindowSize=opt$hws)
message('DONE\n')

message('Combining markers X spectra with correlation and lag...')
pos_metadata = match(
  paste0(cor_data$sample_name, '_', cor_data$replicate),
  paste0(metadata$sample_name, '_', metadata$replicate))
pos_markers = match(cor_data$pept_id, markers_zooms$pept_id)

cor_data = bind_cols(
  cor_data,
  markers_zooms[pos_markers,] %>% rename(marker_taxid = taxid) %>%
    select(-pept_id),
  metadata[pos_metadata,] %>% rename(sample_taxid = taxid) %>%
    select(-sample_name, -replicate))
message('DONE\n')

message('Saving object in %s...', file.path(opt$outfolder, paste0('cordata_',opt$outbase,'.rds')))
saveRDS(cor_data, file.path(opt$outfolder, paste0('cordata_',opt$outbase,'.rds')))
message('DONE\n')
