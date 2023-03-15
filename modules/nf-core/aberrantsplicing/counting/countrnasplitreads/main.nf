process COUNTRNA_SPLITREADS_SAMPLEWISE {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData
        val genomeAssembly
        val recount
        val keepNonStandardChrs

    output:
        tuple val(groups), path("FRASER_output")     , emit: result

    shell:
        groups   = groupData.group
        sampleID = groupData.sampleID
        output   = groupData.output

        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(rmarkdown)
                library(knitr)
                library(devtools)
                library(yaml)
                library(BBmisc)
                library(GenomicAlignments)
                library(tidyr)
                library(data.table)
                library(dplyr)
                library(plotly)
                library(DelayedMatrixStats)
                library(FRASER)
                library(rhdf5)
                library(purrr)
                library(base)
                library(BSgenome)
            })

            dataset    <- "!{groups}"
            genomeAssembly <- "!{genomeAssembly}"

            # Read FRASER object
            file.copy("!{output}", ".", recursive=TRUE)
            fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))

            # Get sample id from wildcard
            sample_id <- "!{sampleID}"

            # If data is not strand specific, add genome info
            genome <- NULL
            if(strandSpecific(fds) == 0){
                genome <- getBSgenome(genomeAssembly)
            }

            # Count splitReads for a given sample id
            sample_result <- countSplitReads(sampleID = sample_id,
                                            fds = fds,
                                            recount = "!{recount}",
                                            keepNonStandardChromosomes = "!{keepNonStandardChrs}",
                                            genome = genome)

            message(date(), ": ", dataset, ", ", sample_id,
                    " no. splice junctions (split counts) = ", length(sample_result))

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
