process PREPROCESSGENEANNOTATION {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
//      tuple val(meta),
        each gtfData

    output:
        tuple val(annotation), path("txdb.db"), \
            path("gene_name_map.tsv"), path("count_ranges.Rds") , emit: result
        path "versions.yml"                                     , emit: versions

    script:
        gtf = gtfData.path
        annotation = gtfData.version
        """
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(GenomicFeatures)
                library(data.table)
                library(rtracklayer)
                library(magrittr)
            })

            ## Create txdb
            txdb <- makeTxDbFromGFF("$gtf")
            txdb <- keepStandardChromosomes(txdb)
            saveDb(txdb, "txdb.db")

            # save count ranges
            count_ranges <- exonsBy(txdb, by = "gene")
            saveRDS(count_ranges, "count_ranges.Rds")

            # Get a gene annotation table
            gtf_dt <- import("$gtf") %>% as.data.table
            if (!"gene_name" %in% colnames(gtf_dt)) {
            gtf_dt[, gene_name := gene_id]
            }
            gtf_dt <- gtf_dt[type == 'gene']

            if('gene_biotype' %in% colnames(gtf_dt))
            setnames(gtf_dt, 'gene_biotype', 'gene_type')

            # Subset to the following columns only
            columns <- c('seqnames', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type', 'gene_status')
            columns <- intersect(columns, colnames(gtf_dt))
            gtf_dt <- gtf_dt[, columns, with = FALSE]

            # make gene_names unique
            gtf_dt[, N := 1:.N, by = gene_name] # warning message
            gtf_dt[, gene_name_orig := gene_name]
            gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
            gtf_dt[, N := NULL]

            fwrite(gtf_dt, "gene_name_map.tsv", na = NA)

            # run the version part
            cat(file="versions.yml", "${task.process}:\naberrantexpression: 1.3.0")
        """
}
