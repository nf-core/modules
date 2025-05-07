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
    tuple val(meta), path("results/*.txt"), emit: maxquant_txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sed \"s_<numThreads>.*_<numThreads>$task.cpus</numThreads>_\" ${paramfile} > mqpar_changed.xml
    sed -i \"s|PLACEHOLDER|\$PWD/|g\" mqpar_changed.xml

    mkdir temp results
    maxquant \\
        ${args} \\
        mqpar_changed.xml
    mv combined/txt/*.txt results/

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            maxquant: \$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \")
    END_VERSIONS
    """

    stub:
    """
    mkdir results
    touch 'results/Oxidation (M)Sites.txt'
    touch results/allPeptides.txt
    touch results/evidence.txt
    touch results/matchedFeatures.txt
    touch results/modificationSpecificPeptides.txt
    touch results/ms3Scans.txt
    touch results/msScans.txt
    touch results/msms.txt
    touch results/msmsScans.txt
    touch results/mzRange.txt
    touch results/parameters.txt
    touch results/peptides.txt
    touch results/proteinGroups.txt
    touch results/summary.txt

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            maxquant: \$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \")
    END_VERSIONS
    """
}
