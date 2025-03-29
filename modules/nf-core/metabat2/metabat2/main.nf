process METABAT2_METABAT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.17--hd498684_0' :
        'biocontainers/metabat2:2.17--hd498684_0' }"

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    tuple val(meta), path("*.tooShort.fa.gz")                    , optional:true, emit: tooshort
    tuple val(meta), path("*.lowDepth.fa.gz")                    , optional:true, emit: lowdepth
    tuple val(meta), path("*.unbinned.fa.gz")                    , optional:true, emit: unbinned
    tuple val(meta), path("*.tsv.gz")                            , optional:true, emit: membership
    tuple val(meta), path("*[!lowDepth|tooShort|unbinned].fa.gz"), optional:true, emit: fasta
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def clean_depth      = depth.toString() - ~/\.gz$/
    def decompress_depth = (depth && depth.toString() != clean_depth) ? "gzip -d -f $depth" : ""
    def depth_input      = depth ? "-a ${clean_depth}" : ""
    """
    $decompress_depth

    metabat2 \\
        ${args} \\
        -i $fasta \\
        ${depth_input} \\
        -t $task.cpus \\
        --saveCls \\
        -o ${prefix}

    gzip -cn ${prefix} > ${prefix}.tsv.gz
    find . -name "*.fa" -type f | xargs -t -n 1 bgzip -@ ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def args             = task.ext.args   ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def decompress_depth = depth           ? "gzip -d -f $depth"    : ""
    def depth_file       = depth           ? "-a ${depth.baseName}" : ""
    """
    echo "" | gzip -c > ${prefix}.1.fa.gz
    echo "" | gzip -c > ${prefix}.1.tooShort.fa.gz
    echo "" | gzip -c > ${prefix}.1.lowDepth.fa.gz
    echo "" | gzip -c > ${prefix}.1.unbinned.fa.gz
    echo "" | gzip -c > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}
