process MINIA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community.wave.seqera.io/library/minia:3.2.6--a3c6327b13b0c75a' :
        'community.wave.seqera.io/library/minia:3.2.6--a3c6327b13b0c75a' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contigs.fa.gz'), emit: contigs
    tuple val(meta), path('*.unitigs.fa.gz'), emit: unitigs
    tuple val(meta), path('*.h5')           , emit: h5
    tuple val(meta), path("*-minia.log")    , emit: log
    tuple val("${task.process}"), val("minia"), eval("minia -v | grep Minia | sed 's/Minia version //g'"), topic: versions, emit: versions_minia

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_list = reads.join(",")
    """
    echo "${read_list}" | sed 's/,/\\n/g' > input_files.txt
    minia \\
        $args \\
        -nb-cores $task.cpus \\
        -in input_files.txt \\
        -out $prefix > ${prefix}-minia.log 2>&1

    if [ -f ${prefix}.contigs.fa ]; then
        gzip -cn ${prefix}.contigs.fa > ${prefix}.contigs.fa.gz
    fi
    if [ -f ${prefix}.unitigs.fa ]; then
        gzip -cn ${prefix}.unitigs.fa > ${prefix}.unitigs.fa.gz
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.contigs.fa.gz
    echo "" | gzip > ${prefix}.unitigs.fa.gz
    touch ${prefix}.h5
    touch ${prefix}-minia.log
    """
}
