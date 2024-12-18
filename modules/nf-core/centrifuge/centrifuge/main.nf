process CENTRIFUGE_CENTRIFUGE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4.1--hdcf5f25_1' :
        'biocontainers/centrifuge:1.0.4.1--hdcf5f25_1' }"

    input:
    tuple val(meta), path(reads)
    path db
    val save_unaligned
    val save_aligned

    output:
    tuple val(meta), path("*report.txt")                 , emit: report
    tuple val(meta), path("*results.txt")                , emit: results
    tuple val(meta), path("*.{sam,tab}")                 , optional: true, emit: sam
    tuple val(meta), path("*.mapped.fastq{,.1,.2}.gz")   , optional: true, emit: fastq_mapped
    tuple val(meta), path("*.unmapped.fastq{,.1,.2}.gz") , optional: true, emit: fastq_unmapped
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-U ${reads}" :  "-1 ${reads[0]} -2 ${reads[1]}"
    def unaligned = ''
    def aligned = ''
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        aligned = save_aligned ? "--al-gz ${prefix}.mapped.fastq.gz" : ''
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        aligned = save_aligned ? "--al-conc-gz ${prefix}.mapped.fastq.gz" : ''
    }
    """
    ## we add "-no-name ._" to ensure silly Mac OSX metafiles files aren't included
    db_name=`find -L ${db} -name "*.1.cf" -not -name "._*"  | sed 's/\\.1.cf\$//'`

    ## make a directory for placing the pipe files in somewhere other than default /tmp
    ## otherwise get pipefile name clashes when multiple centrifuge runs on same node
    ## use /tmp at the same time
    mkdir ./temp

    centrifuge \\
        -x \$db_name \\
        --temp-directory ./temp \\
        -p $task.cpus \\
        $paired \\
        --report-file ${prefix}.report.txt \\
        -S ${prefix}.results.txt \\
        $unaligned \\
        $aligned \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version  | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-U ${reads}" :  "-1 ${reads[0]} -2 ${reads[1]}"
    def unaligned = ''
    def aligned = ''
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        aligned = save_aligned ? "--al-gz ${prefix}.mapped.fastq.gz" : ''
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        aligned = save_aligned ? "--al-conc-gz ${prefix}.mapped.fastq.gz" : ''
    }
    """
    touch ${prefix}.report.txt
    touch ${prefix}.results.txt
    touch ${prefix}.sam
    echo | gzip -n > ${prefix}.unmapped.fastq.gz
    echo | gzip -n > ${prefix}.mapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version  | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}
