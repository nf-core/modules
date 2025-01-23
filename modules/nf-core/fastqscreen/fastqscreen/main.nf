process FASTQSCREEN_FASTQSCREEN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc53eee7ca23c32220a9662fbb63c67769756544b6d74a1ee85cf439ea79a7ee/data' :
        'community.wave.seqera.io/library/fastq-screen_perl-gdgraph:5c1786a5d5bc1309'}"

    input:
    tuple val(meta), path(reads, arity: '1..2')
    tuple val(ref_names), path(ref_dirs, name:"ref*"), val(ref_basenames), val(ref_aligners)

    output:
    tuple val(meta), path("*.txt")     , emit: txt
    tuple val(meta), path("*.png")     , emit: png  , optional: true
    tuple val(meta), path("*.html")    , emit: html
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def config_content = ref_names.withIndex().collect { name, i -> "DATABASE ${name} ./${ref_dirs[i]}/${ref_basenames[i]} ${ref_aligners[i]}" }.join('\n')
    """
    echo '${config_content}' > fastq_screen.conf

    fastq_screen \\
        --conf fastq_screen.conf \\
        --threads ${task.cpus} \\
        $reads \\
        $args

    mv *_screen.txt ${prefix}_screen.txt
    mv *_screen.html ${prefix}_screen.html
    mv *_screen.png ${prefix}_screen.png

    fastq_screen_version=\$(fastq_screen --version 2>&1 | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    echo "\\\"${task.process}\\\":" > versions.yml
    echo "    fastqscreen: \$fastq_screen_version" >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_screen.html
    touch ${prefix}_screen.png
    touch ${prefix}_screen.txt

    fastq_screen_version=\$(fastq_screen --version 2>&1 | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    echo "\\\"${task.process}\\\":" > versions.yml
    echo "    fastqscreen: \$fastq_screen_version" >> versions.yml
    """

}
