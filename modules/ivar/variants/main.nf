process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/ivar:1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  gff

    output:
    tuple val(meta), path("*.tsv")    , emit: tsv
    tuple val(meta), path("*.mpileup"), optional:true, emit: mpileup
    path "versions.yml"               , emit: versions

    script:
    def prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def save_mpileup = params.save_mpileup ? "tee ${prefix}.mpileup |" : ""
    def features     = params.gff ? "-g $gff" : ""
    """
    samtools mpileup \\
        $options.args2 \\
        --reference $fasta \\
        $bam | \\
        $save_mpileup  \\
        ivar variants \\
            $options.args \\
            $features \\
            -r $fasta \\
            -p $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//')
    END_VERSIONS
    """
}
