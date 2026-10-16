process PANTERAGA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6e35ae32df658707535eea50f66f206e576ec8e949325cd2bd5a774a63876d33/data'
        : 'community.wave.seqera.io/library/panteraga:1.4.0--0ec78c92c1b49c00'}"

    input:
    tuple val(meta), path(genomes, stageAs: "genomes/*")
    val seed

    output:
    tuple val(meta), path("*-pantera-final-pass.fa"), emit: te_consensus_passed_fasta
    tuple val(meta), path("*-pantera-final.fa"), emit: te_consensus_fasta
    tuple val(meta), path("*-pantera-final.stats.tsv"), emit: te_stats
    tuple val(meta), path("alignments/*.maf"), emit: alignments
    tuple val(meta), path("annotations/*.bed"), emit: annotations
    tuple val(meta), path("*.1aln"), emit: onealn, optional: true
    tuple val(meta), path("*.log"), emit: panteraga_log, optional: true
    // xargs trims whitespace that was being emitted
    tuple val("${task.process}"), val('panteraga'), eval("panteraGA --version | sed 's/panteraGA //' | xargs"), topic: versions, emit: versions_panteraga

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seed_arg = args.contains("--seed") || args.contains("-r") ? "" : seed ? "--seed ${seed}" : ""
    """
    for file in genomes/*; do
        echo \${file} >> genomes.txt
    done

    panteraGA \\
        ${args} \\
        ${seed_arg} \\
        --genomes genomes.txt \\
        -T ${task.cpus} \\
        -b ${prefix} \\
        -o .

    mv pantera.log ${prefix}.pantera.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir alignments
    mkdir annotations
    touch ${prefix}-pantera-final-pass.fa
    touch ${prefix}-pantera-final.fa
    touch ${prefix}-pantera-final.stats.tsv

    for f in Copia_1.test.LTR.Copia.maf DNA_4.test.DNA.maf I.Jockey_1.test.LINE.I.Jockey.maf; do
        f=\${f/test/${prefix}}
        touch alignments/\${f}
    done

    for genome in genomes/*; do
        genome=\$(basename \${genome})
        touch annotations/\${genome}.bed
    done
    """
}
