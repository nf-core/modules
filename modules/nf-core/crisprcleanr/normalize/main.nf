process CRISPRCLEANR_NORMALIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-crisprcleanr:3.0.0--r42hdfd78af_1':
        'biocontainers/r-crisprcleanr:3.0.0--r42hdfd78af_1' }"

    input:
    tuple val(meta), path(count_file), path(library_file)
    val(min_reads)
    val(min_targeted_genes)

    output:
    tuple val(meta), path("*_norm_table.tsv"), emit: norm_count_file
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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

    version_crisprcleanr <- as.character(packageVersion("CRISPRcleanR"))
    version_rbase <- paste(R.version[['major']],R.version[['minor']], sep = ".")
    writeLines(c(
            '"${task.process}":',
            paste('    r-base:', version_rbase),
            paste('    crisprcleanr:', version_crisprcleanr)
    ), 'versions.yml')

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_norm_table.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: R --version | sed '1!d; s/.*version //; s/ .*//'
        crisprcleanr: Rscript -e 'packageVersion("CRISPRcleanR")'
    END_VERSIONS
    """

}
