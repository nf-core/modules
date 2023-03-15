process AUTOENCODER {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData
        val params

    output:
        tuple val(groups), path("FRASER_output"),   emit: result

    shell:
        groups      = groupData.group
        FR_output   = groupData.output

        implementation = params.implementation
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
            })

            options("FRASER.maxSamplesNoHDF5"=1)
            options("FRASER.maxJunctionsNoHDF5"=-1)

            h5disableFileLocking()

            # copy FRASER_output
            file.copy("!{FR_output}", ".", recursive = TRUE)

            dataset    <- "!{groups}"
            workingDir <- "FRASER_output"

            register(MulticoreParam(1))
            # Limit number of threads for DelayedArray operations
            setAutoBPPARAM(MulticoreParam(1))

            fds <- loadFraserDataSet(dir=workingDir, name=dataset)

            # Fit autoencoder
            # run it for every type
            implementation <- "!{implementation}"

            for(type in psiTypes){
                currentType(fds) <- type
                q <- bestQ(fds, type)
                verbose(fds) <- 3   # Add verbosity to the FRASER object
                fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
                fds <- saveFraserDataSet(fds)
            }
        '''
}
