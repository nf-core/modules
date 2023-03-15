process STATSAE {
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
        maxTestedDimensionProportion = params.maxTestedDimensionProportion
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
            
            # Load Zscores data
            fds <- loadFraserDataSet(dir=workingDir, name=dataset)

            # Calculate stats
            for (type in psiTypes) {
                # Zscores
                fds <- calculateZscore(fds, type=type)
                # Pvalues
                fds <- calculatePvalues(fds, type=type)
                # Adjust Pvalues
                fds <- calculatePadjValues(fds, type=type)
            }

            fds <- saveFraserDataSet(fds)
        '''
}
