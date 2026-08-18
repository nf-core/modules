process CIRI2 {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/69/693c63dc2bc67a4dc24baeee0fd6e1b3a480022c379b4ab6e33d77be917c504e/data':
        'community.wave.seqera.io/library/ciri2_samtools:2e503c1b320bcb93' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(annotation)
    tuple val(meta4), path(ref_dir)

    output:
    tuple val(meta), path("*.txt"), emit: circrna
    tuple val(meta), path("*.txt.log"), emit: log, optional: true
    tuple val(meta), path("CIRIerror.log"), emit: error_log, optional: true
    tuple val("${task.process}"), val('CIRI2.pl'), eval("CIRI2.pl --help | sed -n 's/^Version:\\s*//p'"), topic: versions, emit: versions_ciri2
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cram_ref = fasta ? "-T ${fasta}" : ""
    def anno = annotation ? "--anno ${annotation}" : ""
    def reference = fasta ? "--ref_file ${fasta}" : "--ref_dir ${ref_dir}"

    if (input.getExtension() == 'cram' && !fasta) {
      error "CIRI2 module: a reference fasta is required when input is a CRAM file"
  }

    def is_sam = input.getExtension() == 'sam'
    def sam_file = is_sam ? input.name : "${prefix}.sam"
    def convert_cmd = is_sam ? '' : "samtools view -h $cram_ref $input > ${sam_file}"
    def cleanup_cmd = is_sam ? '' : "rm ${sam_file}"
    """
    $convert_cmd
    CIRI2.pl \\
        --in ${sam_file} \\
        --out ${prefix}.txt \\
        $reference \\
        $anno \\
        --thread_num $task.cpus \\
        $args
    $cleanup_cmd
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}.txt
    """
}
