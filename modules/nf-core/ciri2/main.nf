
process CIRI2 {
    tag "$meta.id"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciri2:2.0.6--pl5321hdfd78af_0':
        'quay.io/biocontainers/ciri2:2.0.6--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(sam) // Required input sam_file
    path  fasta     // optional reference fasta file
    path  annotation        // optional GTF or GFF3 annotation file
    path  ref_dir      // optional reference directory path

    output:

    tuple val(meta), path("*.txt"), emit: circrna
    tuple val(meta), path("*.txt.log"), emit: log, optional: true
    tuple val(meta), path("CIRIerror.log"), emit: error_log, optional: true

    tuple val("${task.process}"), val('CIRI2.pl'), eval("CIRI2.pl --help | sed -n 's/^Version:\\s*//p'"), topic: versions, emit: versions_ciri2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def anno = annotation ? "--anno ${annotation}" : ""
    def reference = fasta ? "--ref_file ${fasta}" : "--ref_dir ${ref_dir}"

    """
    CIRI2.pl \\
        --in $sam \\
        --out ${prefix}.txt \\
        $reference \\
        $anno \\
        --thread_num $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.txt
    """
}
