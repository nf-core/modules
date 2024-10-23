process HMMER_HMMRANK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(tblouts)      // HMMER_HMMSEARCH.out.target_summary

    output:
    tuple val(meta), path("*.hmmrank.tsv.gz"), emit: hmmrank
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # Read all the tblout files

    read_fwf(c('${tblouts.join("','")}'), fwf_cols(content = c(1, NA)), col_types = cols(content = col_character()), comment='#', id = 'fname') %>%
        filter(! str_detect(content, '^ *#')) %>%
        separate(
            content,
            c('accno', 't0', 'profile_desc', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'),
            '\\\\s+',  extra='merge', convert = FALSE
        ) %>%
        transmute(profile = basename(fname) %>% str_remove('${prefix}\\\\.') %>% str_remove('.tbl.gz'), accno, profile_desc, evalue = as.double(evalue), score = as.double(score)) %>%
        # Group and calculate a rank based on score and evalue; let ties be resolved by profile in alphabetical order
        group_by(accno) %>%
        arrange(desc(score), evalue, profile) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        write_tsv('${prefix}.hmmrank.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo \"profile\taccno\tprofile_desc\tevalue\tscore\trank\" | gzip -c > ${prefix}.hmmrank.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.0
        tidyverse: 2.0
    END_VERSIONS
    """
}
