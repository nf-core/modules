process REPEATMASKER_RMOUTTOGFF3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92d028817d6ba2afc60f0d7550b7b510395b3b5426420b7e1c6666ee38ebd5a1/data'
        : 'community.wave.seqera.io/library/repeatmasker:4.1.5--3cc5e72bc3d6ed07'}"

    input:
    tuple val(meta), path(out)

    output:
    tuple val(meta), path("*.gff3") , emit: gff3
    tuple val("${task.process}"), val('repeatmasker'), eval("RepeatMasker -v | sed 's/RepeatMasker version //1'"), topic: versions, emit: versions_repeatmasker

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    rm_path=\$(dirname \$(realpath \$(which RepeatMasker)))

    PERL5LIB=\$rm_path rmOutToGFF3.pl \\
        $out \\
        > ${prefix}.gff3
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    """
}
