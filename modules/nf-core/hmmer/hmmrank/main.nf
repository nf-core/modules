process HMMER_HMMRANK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'quay.io/biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(tblouts), path(domtblouts)    // HMMER_HMMSEARCH.out.target_summary and, optionally, out.domain_summary

    output:
    tuple val(meta), path("*.hmmrank.tsv.gz"), emit: hmmrank
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def domtbl_list = domtblouts ? "c('${domtblouts.join("','")}')" : ''
    def read_domtblouts = domtblouts ? """
    # Read the domtblout files and reduce each hit's domains to one row per coordinate set: the
    # outer bounds of the hit, plus the size of the union of its domains. A hit is keyed on the
    # query name as well as the file, since one HMM file may hold several models and the whole
    # file is then searched in one go, putting several models in the same table.

    islands <- function(d, set) {
        d %>%
            select(accno, profile, query, f = !!sym(paste0(set, '_from')), t = !!sym(paste0(set, '_to'))) %>%
            arrange(accno, profile, query, f) %>%
            group_by(accno, profile, query) %>%
            # Sorted by start, so anything overlapping a row lies before it. cummax(t) is how far
            # right the rows up to here reach; lag() shifts that so each row sees how far
            # everything before it reached. A start past that opens a new island, and since
            # starts are sorted no later domain can fill the gap either. Coordinates are
            # inclusive, so a domain starting exactly one position beyond the previous reach
            # continues it rather than opening a gap: hence the + 1L. Note the default must be
            # 0L, not 0: lag() casts it to the column type and errors on a lossy double cast.
            mutate(island = cumsum(f > lag(cummax(t), default = 0L) + 1L)) %>%
            group_by(accno, profile, query, island) %>%
            summarise(s = min(f), e = max(t), .groups = 'drop_last') %>%
            summarise(
                from = min(s), to = max(e), len = sum(e - s + 1), n_islands = n(),
                .groups = 'drop'
            ) %>%
            rename_with(~ paste0(set, '_', .x), c(from, to, len, n_islands))
    }

    domains <- read_fwf(
        ${domtbl_list}, fwf_cols(content = c(1, NA)),
        col_types = cols(content = col_character()), comment = '#', id = 'fname'
    ) %>%
        filter(! str_detect(content, '^ *#')) %>%
        separate(
            content,
            c(
                'accno', 'd0', 'tlen', 'query', 'd1', 'qlen', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10',
                'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'd11', 'rest'
            ),
            '\\\\s+', extra = 'merge', convert = FALSE
        ) %>%
        transmute(
            profile = basename(fname) %>% str_remove('^${prefix}\\\\.') %>% str_remove('\\\\.domtbl\\\\.gz\$'),
            accno, query,
            across(c(tlen, qlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to), as.integer)
        )

    domain_coords <- domains %>%
        distinct(accno, profile, query, tlen, qlen) %>%
        left_join(islands(domains, 'hmm'), by = c('accno', 'profile', 'query')) %>%
        left_join(islands(domains, 'ali'), by = c('accno', 'profile', 'query')) %>%
        left_join(islands(domains, 'env'), by = c('accno', 'profile', 'query'))
""" : ''
    def join_domtblouts = domtblouts ? "        left_join(domain_coords, by = c('accno', 'profile', 'profile_desc' = 'query')) %>%\n" : ''

    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
${read_domtblouts}
    # Read all the tblout files

    read_fwf(c('${tblouts.join("','")}'), fwf_cols(content = c(1, NA)), col_types = cols(content = col_character()), comment='#', id = 'fname') %>%
        filter(! str_detect(content, '^ *#')) %>%
        separate(
            content,
            c('accno', 't0', 'profile_desc', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'),
            '\\\\s+',  extra='merge', convert = FALSE
        ) %>%
        transmute(profile = basename(fname) %>% str_remove('^${prefix}\\\\.') %>% str_remove('\\\\.tbl\\\\.gz\$'), accno, profile_desc, evalue = as.double(evalue), score = as.double(score)) %>%
        # Group and calculate a rank based on score and evalue; let ties be resolved by profile in alphabetical order
        group_by(accno) %>%
        arrange(desc(score), evalue, profile) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
${join_domtblouts}        write_tsv('${prefix}.hmmrank.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    r-base: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    r-tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def coord_columns = domtblouts ? '\ttlen\tqlen\thmm_from\thmm_to\thmm_len\thmm_n_islands\tali_from\tali_to\tali_len\tali_n_islands\tenv_from\tenv_to\tenv_len\tenv_n_islands' : ''
    """
    echo 'profile\taccno\tprofile_desc\tevalue\tscore\trank${coord_columns}'  > ${prefix}.hmmrank.tsv
    gzip ${prefix}.hmmrank.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
