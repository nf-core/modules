process KINSHIPREAD {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/578b1b874f2cfe4a99acbb2813b84653d81c2756ff1990d46b86a06d6e3b2170/data':
        'community.wave.seqera.io/library/kinship-read:2.1.1--b828062a2eb94ec6' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${prefix}_Read_Results.tsv"), emit: results
    tuple val(meta), path("${prefix}_meansP0_AncientDNA_normalized_READv2.tsv"), emit: normalised_means
    tuple val(meta), path("${prefix}_Read_Results.pdf"), emit: results_pdf, optional: true
    tuple val("${task.process}"), val('kinship-read'), eval("kinship-read --version | sed 's/v//'"), topic: versions, emit: versions_read

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_prefix = bed.toString().replaceAll(/\.bed$/, '')
    """
    kinship-read \\
        --input ${input_prefix} \\
        $args \\

    mv Read_Results.tsv ${prefix}_Read_Results.tsv
    mv meansP0_AncientDNA_normalized_READv2 ${prefix}_meansP0_AncientDNA_normalized_READv2.tsv
    if [[ -f Read_Results.pdf ]]; then
        mv Read_Results.pdf ${prefix}_Read_Results.pdf
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_prefix = bed.toString().replaceAll(/\.bed$/, '')
    """
    echo $input_prefix
    echo $args
    
    touch ${prefix}_Read_Results.tsv
    touch ${prefix}_meansP0_AncientDNA_normalized_READv2.tsv
    touch ${prefix}_Read_Results.pdf
    """
}
