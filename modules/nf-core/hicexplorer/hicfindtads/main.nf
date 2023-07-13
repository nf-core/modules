process HICEXPLORER_HICFINDTADS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hicexplorer=3.7.2 numpy=1.23.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1':
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(matrix)
    val resolution

    output:
    tuple val(meta), val(resolution), path("*hicfindtads*")      , emit:results
    tuple val(meta), val(resolution), path("*domains.bed")       , emit:tads
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${resolution}"
    ## add defaults for required parameters.
    def bin_size = meta.bin? meta.bin.toInteger() : resolution
    def mDepth = [bin_size * 3, bin_size * 10]
    def mDepthArg = ['--minDepth', '--maxDepth']
    args = args.tokenize()
    for(i=0; i<2; i++){
        idx = args.indexOf(mDepthArg[i])
        if(idx>=0){
            mDepth[i] = args[idx + 1]
            args.remove(idx+1)
            args.remove(idx)
        }
    }
    def correctForMultipleTesting = 'fdr'
    idx = args.indexOf('correctForMultipleTesting')
    if(idx>=0){
        correctForMultipleTesting = args[idx + 1]
        args.remove(idx+1)
        args.remove(idx)
    }
    args = args.join(' ')
    """
    hicFindTADs \\
        ${args} \\
        --matrix $matrix \\
        --outPrefix ${prefix}_hicfindtads \\
        --correctForMultipleTesting $correctForMultipleTesting \\
        --minDepth ${mDepth[0]} \\
        --maxDepth ${mDepth[1]} \\
        --step $resolution \\
        --numberOfProcessors ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicFindTADs --version 2>&1 | sed 's/hicFindTADs //')
    END_VERSIONS
    """
}
