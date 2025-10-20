process PICRUST2_PIPELINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6a/6a9b732fec49b33334dcff4f8875dcd552402ead5654443a54c7f4f79823df78/data'
        : 'community.wave.seqera.io/library/picrust2:2.6.2--a7c158f7c987b452'}"

    input:
    tuple val(meta), path(sequences), path(otu_table)

    output:
    tuple val(meta), path("${prefix}/")
    tuple val(meta), path("${prefix}/*_reduced.tre"), emit: trees
    tuple val(meta), path("${prefix}/*_metagenome_out/pred_metagenome_unstrat.tsv.gz"), emit: function_abundances
    tuple val(meta), path("${prefix}/pathways_out/path_abun_unstrat.tsv.gz"), emit: pathway_abundances
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    picrust2_pipeline.py \\
        --study_fasta ${sequences} \\
        --input ${otu_table} \\
        --processes ${task.cpus} \\
        --output ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picrust2: \$( picrust2_pipeline.py --version | sed 's/PICRUSt2 //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    mkdir -p ${prefix}/{EC_metagenome_out,pathways_out}
    touch ${prefix}/bac_reduced.tre
    echo '' | gzip -c > ${prefix}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
    echo '' | gzip -c > ${prefix}/pathways_out/path_abun_unstrat.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picrust2: \$( picrust2_pipeline.py --version | sed 's/PICRUSt2 //' )
    END_VERSIONS
    """
}
