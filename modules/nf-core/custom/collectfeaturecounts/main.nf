process CUSTOM_COLLECTFEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/143121bf17b3f4e9539d8cd17aaaa428aae54ce827471c4fc6399465263efa2e/data' :
        'community.wave.seqera.io/library/r-base_r-dplyr_r-readr_r-purrr_pruned:0f879b99d6a89834' }"

    input:
    tuple val(meta), path(inputfiles)

    output:
    tuple val(meta), path("${prefix}.counts.tsv.gz"), emit: counts
    path "versions.yml"                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_list = inputfiles instanceof List ? inputfiles : [inputfiles]
    def r_files = file_list.collect { "'${it}'" }.join(', ')
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    setDTthreads($task.cpus)

    tibble(f = c($r_files)) %>%
        mutate(
            d = purrr::map(
                f,
                function(file) {
                    fread(file, sep = '\\t', skip = 1) %>%
                        melt(measure.vars = c(ncol(.)), variable.name = 'sample', value.name = 'count') %>%
                        lazy_dt() %>%
                        filter(count > 0) %>%
                        mutate(
                            sample = str_remove(sample, '.sorted.bam'),
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
                        as_tibble()
                }
            )
        ) %>%
        tidyr::unnest(d) %>%
        select(-f) %>%
        write_tsv("${prefix}.counts.tsv.gz")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    stringr: ", packageVersion('stringr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')),
            paste0("    data.table: ", packageVersion('data.table'))
        ),
        "versions.yml"
    )
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv
    gzip ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        readr: 2.0.0
        stringr: 1.4.0
        dtplyr: 1.1.0
        data.table: 1.14.0
    END_VERSIONS
    """
}
