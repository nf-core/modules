process FITHYPER {
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


            # Disabled, c-sandu
            # if ("random_seed" %in% names(params)){
            #    rseed <- snakemake@config$random_seed
            #    if(isTRUE(rseed)){
            #        set.seed(42)
            #    } else if (is.numeric(rseed)){
            #        set.seed(as.integer(rseed))
            #    }
            # }

            #+ input
            dataset    <- "!{groups}"
            workingDir <- "FRASER_output"

            register(MulticoreParam(1))
            # Limit number of threads for DelayedArray operations
            setAutoBPPARAM(MulticoreParam(1))

            # Load PSI data
            fds <- loadFraserDataSet(dir=workingDir, name=dataset)

            # Run hyper parameter optimization
            implementation <- "!{implementation}"
            mp <- !{maxTestedDimensionProportion}

            # Get range for latent space dimension
            a <- 2 
            b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

            maxSteps <- 12
            if(mp < 6){
            maxSteps <- 15
            }

            Nsteps <- min(maxSteps, b)
            pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique

            for(type in psiTypes){
                message(date(), ": ", type)
                fds <- optimHyperParams(fds, type=type, 
                                        implementation=implementation,
                                        q_param=pars_q,
                                        plot = FALSE)
                fds <- saveFraserDataSet(fds)
            }
            fds <- saveFraserDataSet(fds)
        '''
}
