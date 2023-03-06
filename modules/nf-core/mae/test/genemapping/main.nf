process GENE_NAME_MAPPING {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each gene_map_data

    output:
        tuple val(annotation), path("output.tsv")     , emit: result

    shell:
        annotation   = gene_map_data.annotation
        gtf          = gene_map_data.gtf
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(rtracklayer)
                library(data.table)
                library(magrittr)
                library(tidyr)
            })

            gtf_dt <- import("!{gtf}") %>% as.data.table
            if (!"gene_name" %in% colnames(gtf_dt)) {
            gtf_dt[gene_name := gene_id]
            }
            if('gene_biotype' %in% colnames(gtf_dt))
            setnames(gtf_dt, 'gene_biotype', 'gene_type')
            gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type)]

            # make gene_names unique
            gtf_dt[, N := 1:.N, by = gene_name] # warning message
            gtf_dt[, gene_name_orig := gene_name]
            gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
            gtf_dt[, N := NULL]

            fwrite(gtf_dt, "output.tsv", na = NA)
        '''
}
