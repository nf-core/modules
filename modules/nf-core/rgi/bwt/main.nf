process RGI_BWT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rgi:6.0.5--pyh05cac1d_0'
        : 'biocontainers/rgi:6.0.5--pyh05cac1d_0'}"

    input:
    tuple val(meta), path(reads, arity: '1..2')
    path card
    path wildcard

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt"), emit: tsv
    tuple val(meta), path("temp/"), emit: tmp
    tuple val("${task.process}"), val('rgi')         , eval("rgi main --version")  , emit: versions_rgi, topic: versions
    tuple val("${task.process}"), val('rgi-database'), eval("echo \$DB_VERSION")   , emit: versions_db , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // This customizes the command: rgi load
    def args = task.ext.args ?: ''
    // This customizes the command: rgi main
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_one = reads[0]
    def read_two = reads.size() > 1 ? reads[1] : null
    def load_wildcard = ""

    if (wildcard) {
        load_wildcard = """ \\
            --wildcard_annotation ${wildcard}/wildcard_database_v\$DB_VERSION.fasta \\
            --wildcard_annotation_all_models ${wildcard}/wildcard_database_v\$DB_VERSION\\_all.fasta \\
            --wildcard_index ${wildcard}/wildcard/index-for-model-sequences.txt \\
            --amr_kmers ${wildcard}/wildcard/all_amr_61mers.txt \\
            --kmer_database ${wildcard}/wildcard/61_kmer_db.json \\
            --kmer_size 61
        """
    }

    """
    DB_VERSION=\$(ls ${card}/card_database_*_all.fasta | sed "s/${card}\\/card_database_v\\([0-9].*[0-9]\\).*/\\1/")

    rgi \\
        load \\
        ${args} \\
        --local \\
        --card_json ${card}/card.json \\
        --debug \\
        --card_annotation ${card}/card_database_v\$DB_VERSION.fasta \\
        --card_annotation_all_models ${card}/card_database_v\$DB_VERSION\\_all.fasta \\
        ${load_wildcard}

    rgi \\
        bwt \\
        ${args} \\
        --local \\
        --threads ${task.cpus} \\
        --output_file ${prefix} \\
        --read_one ${read_one} \\
        ${ read_two ? "--read_two ${read_two}" : "" }


    mkdir temp/
    for FILE in *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json; do [[ -e \$FILE ]] && mv \$FILE temp/; done

    """

    stub:
    """
    mkdir -p temp
    touch test.json
    touch test.txt

    DB_VERSION=stub_version
    """
}
