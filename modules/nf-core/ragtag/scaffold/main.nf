process RAGTAG_SCAFFOLD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/3836fefb293bbbf2be52f760b1833819ae30fe71cd3becf270a4bd5149e5b9b6/data'
        : 'community.wave.seqera.io/library/ragtag_gzip:5de7cc3ac89b3ac4' }"

    input:
    tuple val(meta), path(assembly, name: 'assembly/*')
    tuple val(meta2), path(reference, name: 'reference/*')
    tuple val(meta3), path(exclude)
    tuple val(meta4), path(skip), path(hard_skip)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: corrected_assembly
    tuple val(meta), path("*.agp"), emit: corrected_agp
    tuple val(meta), path("*.stats"), emit: corrected_stats
    tuple val(meta), path("*.paf.gz"), emit: ragtag_paf
    tuple val(meta), path("*.paf.log"), emit: ragtag_paf_log, optional: true
    tuple val(meta), path("*.delta.gz"), emit: ragtag_delta, optional: true
    tuple val(meta), path("*.delta.log"), emit: ragtag_delta_log, optional: true
    tuple val(meta), path("*.confidence.txt"), emit: confidence_txt
    tuple val("${task.process}"), val('ragtag'), eval("ragtag.py -v | sed 's/v//'"), emit: versions_ragtag, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_exclude = exclude ? "-e ${exclude}" : ""
    def arg_skip = skip ? "-j ${skip}" : ""
    def arg_hard_skip = hard_skip ? "-J ${hard_skip}" : ""
    """
    if [[ ${assembly} == *.gz ]]
    then
        zcat ${assembly} > assembly.fa
    else
        ln -s ${assembly} assembly.fa
    fi

    if [[ ${reference} == *.gz ]]
    then
        zcat ${reference} > reference.fa
    else
        ln -s ${reference} reference.fa
    fi

    ragtag.py scaffold reference.fa assembly.fa \\
        -o . \\
        -t ${task.cpus} \\
        ${arg_exclude} \\
        ${arg_skip} \\
        ${arg_hard_skip} \\
        ${args} \\
        2>| >( tee ${prefix}.stderr.log >&2 ) \\
        | tee ${prefix}.stdout.log

    for file in ragtag.scaffold.*; do
        mv \$file \${file/ragtag.scaffold/${prefix}}
    done

    gzip *.fasta 2>/dev/null || true
    gzip *.delta 2>/dev/null || true
    gzip *.paf 2>/dev/null || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fasta.gz
    touch ${prefix}.agp
    touch ${prefix}.stats
    echo "" | gzip > ${prefix}.asm.paf.gz
    echo "" | gzip > ${prefix}.asm.delta.gz
    touch ${prefix}.delta.log
    touch ${prefix}.paf.log
    touch ${prefix}.confidence.txt
    """
}
