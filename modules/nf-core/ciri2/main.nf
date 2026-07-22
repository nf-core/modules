
process CIRI2 {
    tag "$meta.id"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciri2:2.0.6--pl5321hdfd78af_0':
        'quay.io/biocontainers/ciri2:2.0.6--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(sam) // Required -I sam_file
    path  fasta     // optional -F reference fasta file
    path  annotation        // optional -A GTF or GFF3 annotation file
    path  ref_dir      // optional -R reference directory path

    output:

    tuple val(meta), path("*.txt"), emit: circrna
    tuple val(meta), path("*.txt.log"), emit: logfile, optional: true
    tuple val(meta), path("CIRIerror.log"), emit: error_log, optional: true


    tuple val("${task.process}"), val('CIRI2.pl'), val("2.0.6"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def anno = annotation ? "-A ${annotation}" : ""
    def reference = fasta ? "-F ${fasta}" : "-R ${ref_dir}"

    """
    CIRI2.pl \\
        -I $sam \\
        -O ${prefix}.txt \\
        $reference \\
        $anno \\
        -T $task.cpus \\
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
