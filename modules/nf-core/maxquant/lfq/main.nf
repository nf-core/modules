process MAXQUANT_LFQ {
    tag "$meta.id"
    label 'process_long'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/maxquant:2.0.3.0--py310hdfd78af_1' :
        'biocontainers/maxquant:2.0.3.0--py310hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta), path(paramfile)
    path(raw)

    output:
    tuple val(meta), path("${prefix}/*.txt"), emit: maxquant_txt
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    sed \"s_<numThreads>.*_<numThreads>$task.cpus</numThreads>_\" ${paramfile} > mqpar_changed.xml
    sed -i \"s|PLACEHOLDER|\$PWD/|g\" mqpar_changed.xml

    mkdir ${prefix}
    maxquant \\
        ${args} \\
        mqpar_changed.xml
    mv combined/txt/*.txt ${prefix}/

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            maxquant: \$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \")
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch '${prefix}/Oxidation (M)Sites.txt'
    touch ${prefix}/allPeptides.txt
    touch ${prefix}/evidence.txt
    touch ${prefix}/matchedFeatures.txt
    touch ${prefix}/modificationSpecificPeptides.txt
    touch ${prefix}/ms3Scans.txt
    touch ${prefix}/msScans.txt
    touch ${prefix}/msms.txt
    touch ${prefix}/msmsScans.txt
    touch ${prefix}/mzRange.txt
    touch ${prefix}/parameters.txt
    touch ${prefix}/peptides.txt
    touch ${prefix}/proteinGroups.txt
    touch ${prefix}/summary.txt

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            maxquant: \$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \")
    END_VERSIONS
    """
}
