#!/usr/bin/env Rscript

library(data.table)
library(dtplyr)
library(readr)
library(dplyr)
library(stringr)

setDTthreads($task.cpus)

prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

files <- Sys.glob('input/*')
if ( length(files) == 0 ) stop('No input files')

# One file at a time: memory use is set by the largest file, not by the number of samples.
# TPM is computed per file, so each sample must be in exactly one file.
out <- gzfile(paste0(prefix, '.counts.tsv.gz'), 'wb')
first <- TRUE
seen <- character()

for ( file in files ) {
    sample <- str_remove(tail(names(fread(file, sep = '\t', skip = 1, nrows = 0)), 1), '[.]sorted[.]bam\$')
    if ( sample %in% seen ) stop(paste0('Sample ', sample, ' is in more than one input file'))
    seen <- c(seen, sample)
    fread(file, sep = '\t', skip = 1) %>%
        melt(measure.vars = c(ncol(.)), variable.name = 'sample', value.name = 'count') %>%
        lazy_dt() %>%
        filter(count > 0) %>%
        mutate(
            sample = str_remove(sample, '[.]sorted[.]bam\$'),
            r = count/Length
        ) %>%
        rename(orf = Geneid, chr = Chr, start = Start, end = End, strand = Strand, length = Length) %>%
        group_by(sample) %>%
        # Rounded so tables independently computed elsewhere from the same underlying
        # counts (e.g. after a caller-consolidation step using a different R backend)
        # stay byte-comparable -- data.table/dtplyr and plain dplyr can otherwise round
        # the last significant digit of the same mathematical value differently.
        mutate(tpm = round(r/sum(r) * 1e6, 6)) %>% ungroup() %>%
        select(-r) %>%
        as_tibble() %>%
        write_tsv(out, col_names = first)
    first <- FALSE
}

close(out)

writeLines(
    c(
        "\"$task.process\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    dplyr: ", packageVersion('dplyr')),
        paste0("    readr: ", packageVersion('readr')),
        paste0("    stringr: ", packageVersion('stringr')),
        paste0("    dtplyr: ", packageVersion('dtplyr')),
        paste0("    data.table: ", packageVersion('data.table'))
    ),
    "versions.yml"
)
