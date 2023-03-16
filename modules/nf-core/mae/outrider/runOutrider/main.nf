process RUNOUTRIDER {
//    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each filter_data
        val(implementation)
        val(maxTestedDimensionProportion)

    output:
        tuple val(preprocess), path("ods.Rds")      , emit: result
        path "versions.yml"                         , emit: versions

    shell:
        preprocess = filter_data.preprocess
        ods_unfitted = filter_data.odsUnfit
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(OUTRIDER)
                library(SummarizedExperiment)
                library(ggplot2)
                library(data.table)
                library(dplyr)
                library(magrittr)
                library(tools)
            })

            ods <- readRDS("!{ods_unfitted}")
            implementation <- "!{implementation}"
            mp <- !{maxTestedDimensionProportion}

            ## subset filtered
            ods <- ods[mcols(ods)$passedFilter,]

            # add gene ranges to rowData
            gr <- unlist(endoapply(rowRanges(ods), range))
            if(length(gr) > 0){
                rd <- rowData(ods)
                rowRanges(ods) <- gr
                rowData(ods) <- rd
            }

            ods <- estimateSizeFactors(ods)

            ## find optimal encoding dimension
            a <- 5
            b <- min(ncol(ods), nrow(ods)) / mp   # N/3

            maxSteps <- 15
            if(mp < 4){
                maxSteps <- 20
            }

            Nsteps <- min(maxSteps, b)   # Do at most 20 steps or N/3
            # Do unique in case 2 were repeated
            pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
            ods <- findEncodingDim(ods, params = pars_q, implementation = implementation)

            ## fit OUTRIDER
            ods <- OUTRIDER(ods, implementation = implementation)
            message("outrider fitting finished")

            saveRDS(ods, "ods.Rds")

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
