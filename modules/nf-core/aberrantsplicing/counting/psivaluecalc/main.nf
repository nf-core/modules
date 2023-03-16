process PSI_VALUE_CALC {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData

    output:
        tuple val(groups), path("FRASER_output"), \
            path("FRASER_output/savedObjects/raw-local-*/theta.h5")    , emit: result

    shell:
        groups       = groupData.group
        FR_output    = groupData.output
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

            dataset    <- "!{groups}"

            # Read FRASER object
            file.copy("!{FR_output}", ".", recursive = TRUE)
            fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))

            # Calculating PSI values
            fds <- calculatePSIValues(fds)

            # FRASER object after PSI value calculation
            fds <- saveFraserDataSet(fds)

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
