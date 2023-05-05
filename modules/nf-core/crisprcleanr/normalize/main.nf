process CRISPRCLEANR_NORMALIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::r-crisprcleanr=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-crisprcleanr:3.0.0--r42hdfd78af_1':
        'biocontainers/r-crisprcleanr:3.0.0--r42hdfd78af_1' }"

    input:
    tuple val(meta), path(count_file), path(library_file)
    val(min_reads)
    val(min_targeted_genes)

    output:
    tuple val(meta), path("*_norm_table.tsv"), emit: norm_count_file
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(CRISPRcleanR)
    library <- read.delim('${library_file}', header=T,sep="\t")
    row.names(library) <- library[["CODE"]]
    normANDfcs <- ccr.NormfoldChanges('${count_file}',saveToFig = FALSE,min_reads=${min_reads},EXPname='${meta.id}', libraryAnnotation=library,display=FALSE)
    gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs[["logFCs"]],library)
    correctedFCs <- ccr.GWclean(gwSortedFCs,display=FALSE,label='${meta.id}')
    correctedCounts <- ccr.correctCounts('${meta.id}',
                            normANDfcs[["norm_counts"]],
                            correctedFCs,
                            library,
                            minTargetedGenes=${min_targeted_genes},
                            OutDir='./')

    write.table(correctedCounts, file=paste0("${prefix}","_norm_table.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

    version_file_path <- "versions.yml"
    version_crisprcleanr <- paste(unlist(packageVersion("CRISPRcleanR")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    crisprcleanr: ", f, sep = "")
    writeLines(version_crisprcleanr, f)
    close(f)

    """
}
