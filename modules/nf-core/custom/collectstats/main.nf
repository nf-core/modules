// Safely quote a Groovy value as a single-quoted R string literal. Used everywhere
// below that a sample ID or other user-controlled value gets embedded into generated
// R source -- without this, a value containing a quote could break the generated R
// syntax, or worse.
def rq(v) {
    return "'" + v.toString().replace('\\', '\\\\').replace("'", "\\'") + "'"
}

process CUSTOM_COLLECTSTATS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9d1a8065dcebe35c513db5eb4a5237795aa4bae4756086f3c4319a820a8a647c/data'
:         'community.wave.seqera.io/library/custom_collectstats:c1f477e6251c36bb' }"

    input:
    tuple val(meta), val(samples_meta), path(trimlogs), path(bblogs), path(idxstats), path(fcs), path(mergetab)

    output:
    tuple val(meta), path("${outfile}"), emit: overall_stats
    path "versions.yml"                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile    = "${prefix}.overall_stats.tsv.gz"

    def read_trimlogs = ""
    if ( trimlogs ) {
        def se_trimlogs = samples_meta
            .findAll { s -> s.single_end }
            .collect { s -> "${rq(s.id)}, ${rq(s.id + '.*_trimming_report.txt')}, 1," }
        def pe_trimlogs = samples_meta
            .findAll { s -> ! s.single_end }
            .collect { s -> "${rq(s.id)}, ${rq(s.id + '_1.*_trimming_report.txt')}, 2," }

        read_trimlogs = """
        trimming <- tribble(
            ~sample, ~tlogglob, ~mult,
            ${se_trimlogs.join("\n")}
            ${pe_trimlogs.join("\n")}
        ) %>%
            mutate(
                d = map(
                    tlogglob,
                    function(g) {
                        # Deliberately using Sys.glob() + readLines() rather than shelling
                        # out to grep/sed via pipe(sprintf(...)): the glob and file content
                        # are only ever passed as R function arguments here, never
                        # interpolated into a shell command string.
                        f <- Sys.glob(g)[1]
                        line <- grep('Reads written \\\\(passing filters\\\\)', readLines(f), value = TRUE)
                        n <- as.integer(gsub(',', '', sub(' .*', '', sub('.*: *', '', line))))
                        tibble(n_post_trimming = n)
                    }
                )
            ) %>%
            unnest(d) %>%
            transmute(sample, n_post_trimming = n_post_trimming * mult)
        """
    } else {
        read_trimlogs = "trimming <- tibble(sample = character(), n_post_trimming = integer())"
    }

    if ( mergetab ) {
        read_mergetab = """
        mergetab <- read_tsv("${mergetab}", show_col_types = FALSE)
        """
    } else {
        read_mergetab = """
        mergetab <- tibble(sample = character())
        """
    }

    """
    #!/usr/bin/env Rscript

    # All read counts joined into this table (n_post_trimming, n_non_contaminated,
    # idxs_n_mapped/idxs_n_unmapped, and the featureCounts columns) are expected to use
    # the same unit -- individual reads, with both mates of a pair counted separately,
    # never fragments/pairs as a single unit -- so the columns stay directly comparable.
    # This depends on how the upstream trimming/decontamination/alignment/counting steps
    # are invoked (e.g. whether a paired-end-aware counting tool is asked to count reads
    # or fragments), not on anything this module itself controls.

    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(stringr)

    start    <- tibble(sample = c(${samples_meta.collect { s -> rq(s.id) }.join(', ')}))

    ${read_trimlogs}

    # Reading each idxstats file directly instead of shelling out to grep, to avoid
    # passing a sample-derived filename through a shell command string.
    idxs <- tibble(fname = Sys.glob('*.idxstats')) %>%
        mutate(
            sample = str_remove(basename(fname), '.idxstats'),
            d = map(
                fname,
                \\(f) {
                    read_tsv(f, col_names = c('chr', 'length', 'idxs_n_mapped', 'idxs_n_unmapped'), col_types = 'ciii') %>%
                        filter(chr != '*') %>%
                        summarise(idxs_n_mapped = sum(idxs_n_mapped), idxs_n_unmapped = sum(idxs_n_unmapped))
                }
            )
        ) %>%
        unnest(d) %>%
        select(-fname)

    # Column types given by name rather than position: this only depends on the
    # feature-counts input's columns existing and being correctly typed, not also on
    # their order, so a future column reorder upstream won't silently misparse here.
    counts <- read_tsv(
        c('${fcs.join("', '")}'),
        id = 'fname',
        col_types = cols(
            .default = col_guess(),
            sample   = col_character(),
            count    = col_integer()
        )
    ) %>%
        mutate(feature = str_replace(fname, '^[^.]+\\\\.([^.]+)\\\\..*', '\\\\1')) %>%
        group_by(sample, feature) %>%
        summarise(count = sum(count), .groups = 'drop') %>%
        pivot_wider(names_from = feature, values_from = count, values_fill = 0)

    # Reading each bbduk log directly (instead of shelling out to grep/sed) avoids
    # passing a sample-derived filename through a shell command string.
    bbduk <- tibble(sample = character(), n_non_contaminated = integer())
    for ( f in Sys.glob('*.bbduk.log') ) {
        s <- str_remove(f, '.bbduk.log')
        line <- grep('Result:', readLines(f), value = TRUE)
        n <- as.integer(sub(' reads.*', '', sub('Result:[[:space:]]*', '', line)))
        bbduk <- bbduk %>% union(tibble(n_non_contaminated = n, sample = s))
    }
    if ( nrow(bbduk) == 0 ) bbduk <- bbduk %>% select(sample)

    # Add in stats from taxonomy and function, if provided
    ${read_mergetab}

    # Write output
    start %>%
        left_join(trimming, by = join_by(sample)) %>%
        left_join(bbduk, by = join_by(sample)) %>%
        left_join(idxs, by = join_by(sample)) %>%
        left_join(counts, by = join_by(sample)) %>%
        left_join(mergetab, by = join_by(sample)) %>%
        arrange(sample) %>%
        write_tsv("${outfile}")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    purrr: ", packageVersion('purrr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile    = "${prefix}.overall_stats.tsv.gz"
    """
    cat /dev/null | gzip -c > ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        readr: 2.0.0
        purrr: 0.3.4
        tidyr: 1.1.3
        stringr: 1.4.0
    END_VERSIONS
    """
}
