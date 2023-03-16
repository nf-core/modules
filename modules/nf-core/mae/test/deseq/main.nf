process DESEQ {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each allelic_counts
        val (params)

    output:
        tuple val(vcf), val(rnaID), path("output.Rds")     , emit: result

    shell:
        vcf     = allelic_counts.vcf
        rnaID   = allelic_counts.rnaID

        counts  = allelic_counts.counts

        addAF = params.mae.addAF
        maxAF = params.mae.maxAF
        genome_assembly = params.genomeAssembly
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(stringr)
                library(tMAE)
            })

            message("Started with deseq")

            # Read mae counts
            mae_counts <- fread("!{counts}", fill=TRUE)
            mae_counts <- mae_counts[contig != '']
            mae_counts[, position := as.numeric(position)]

            # Sort by chr
            mae_counts <- mae_counts[!grep("Default|opcode", contig)]
            mae_counts[, contig := factor(contig,
                            levels = unique(str_sort(mae_counts$contig, numeric = TRUE)))]


            print("Running DESeq...")
            # Function from tMAE pkg
            rmae <- DESeq4MAE(mae_counts) ## negative binomial test for allelic counts

            ### Add AF information from gnomAD
            if ("!{addAF}" == "true") {
                print("Adding gnomAD allele frequencies...")

                # obtain the assembly from the config
                genome_assembly <- "!{genome_assembly}"
                rmae <- add_gnomAD_AF(rmae, genome_assembly = genome_assembly,
                    max_af_cutoff = !{maxAF}, populations = c("AF", "AF_afr", "AF_amr", "AF_eas", "AF_nfe"))
            } else {
                rmae[, rare := NA]
            }

            saveRDS(rmae, "output.Rds")
        '''
}
