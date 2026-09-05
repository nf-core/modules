process STARSUITE_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://biodepot/star-suite@sha256:0ced4e4941b8a7347c1583598381e6abe25d798df09759c59ede0e6c6dc4e1d3' :
        'docker.io/biodepot/star-suite@sha256:0ced4e4941b8a7347c1583598381e6abe25d798df09759c59ede0e6c6dc4e1d3' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("starsuite"), emit: index
    tuple val("${task.process}"), val('starsuite'), eval('STAR --version'), emit: versions_starsuite, topic: versions
    tuple val("${task.process}"), val('star'), eval('STAR --upstream-version'), emit: versions_star, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "STARSUITE_GENOMEGENERATE does not support Conda. Please use Docker, Singularity, Apptainer, or Podman instead."
    }
    def args   = task.ext.args ?: ''
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir starsuite
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir starsuite/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --runThreadN $task.cpus \\
        $memory \\
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "STARSUITE_GENOMEGENERATE does not support Conda. Please use Docker, Singularity, Apptainer, or Podman instead."
    }
    """
    mkdir starsuite
    touch starsuite/Genome
    touch starsuite/Log.out
    touch starsuite/SA
    touch starsuite/SAindex
    touch starsuite/genomeParameters.txt
    touch starsuite/sjdbList.out.tab
    touch starsuite/transcriptome.fa
    """
}
