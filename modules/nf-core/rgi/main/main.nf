process RGI_MAIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_0':
        'biocontainers/rgi:6.0.3--pyha8f3691_0' }"

    input:
    tuple val(meta), path(fasta)
    val(wildcard)
    path(database)
    
    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt") , emit: tsv
    tuple val(meta), path("temp/") , emit: tmp
    path("database_output/")       , emit: db_dir
    env RGI_VERSION                , emit: tool_version
    env DB_VERSION                 , emit: db_version
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def db_download_cmd = ""
    def db_process_cmd = ""
    def load_wildcard = ""

    if (!database) { // Download CARD DB if no path given
    db_download_cmd = "wget https://card.mcmaster.ca/latest/data;tar -xvf data ./card.json"
    db_process_cmd = "rgi card_annotation -i card.json"
    }

    // if (wildcard.equalsIgnoreCase("yes")) {
    //     load_wildcard = "--wildcard_annotation $database/wildcard_database_v\$DB_VERSION.fasta \\
    //     --wildcard_annotation_all_models $database/wildcard_database_v\$DB_VERSION\\_all.fasta \\
    //     --wildcard_index $database/wildcard/index-for-model-sequences.txt \\
    //     --wildcard_version \$DB_VERSION \\
    //     --amr_kmers $database/wildcard/all_amr_61mers.txt \\
    //     --kmer_database $database/wildcard/61_kmer_db.json \\
    //     --kmer_size 61"
    //     if (!database) { // Download both DBs if no path is given
    //         db_download_cmd = "wget https://card.mcmaster.ca/latest/data;tar -xvf data ./card.json;wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants;mkdir -p wildcard;tar -xjf wildcard_data.tar.bz2 -C wildcard;gunzip wildcard/*.gz"
    //         db_process_cmd = "rgi card_annotation -i card.json;rgi wildcard_annotation -i wildcard --card_json card.json -v \$DB_VERSION"
    //         database = "."
    // }

    """
    $db_download_cmd
    DB_VERSION=\$(ls $database/card_database_*_all.fasta | sed 's/card_database_v\\([0-9].*[0-9]\\).*/\\1/')
    $db_process_cmd
    
    rgi \\
        load \\
        $args \\
        --card_json $database/card.json \\
        --debug --local \\
        --card_annotation $database/card_database_v\$DB_VERSION.fasta \\
        --card_annotation_all_models $database/card_database_v\$DB_VERSION\\_all.fasta \\
        $load_wildcard

    rgi \\
        main \\
        $args \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    mkdir temp/
    mv *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json temp/

    mkdir database_output
    mv $database/* database_output

    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """
}
