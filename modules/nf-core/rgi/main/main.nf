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
    path("database_output/")       , emit: db_out
    env RGI_VERSION                , emit: tool_version
    env DB_VERSION                 , emit: db_version
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def db_download_cmd = ""
    def db_card_cmd = ""
    def db_wildcard_cmd = ""
    def load_wildcard = ""
    def copy_db_cmd = "cp $database/card* database_output"

    if (wildcard.equalsIgnoreCase("yes")) {
        if (!database) { // Download both DBs if no path is given
            db_download_cmd = "wget https://card.mcmaster.ca/latest/data;tar -xvf data ./card.json -C database_output;wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants;mkdir -p wildcard;tar -xjf wildcard_data.tar.bz2 -C wildcard;gunzip wildcard/*.gz"
            db_card_cmd = "rgi card_annotation -i database_output/card.json > card_annotation.log 2>&1"
            db_wildcard_cmd = "rgi wildcard_annotation -i wildcard --card_json database_output/card.json -v \$DB_VERSION > wildcard_annotation.log 2>&1"
            database = "."
        }
        load_wildcard = """ \\
            --wildcard_annotation database_output/wildcard_database_v\$DB_VERSION.fasta \\
            --wildcard_annotation_all_models database_output/wildcard_database_v\$DB_VERSION\\_all.fasta \\
            --wildcard_index database_output/wildcard/index-for-model-sequences.txt \\
            --amr_kmers database_output/wildcard/all_amr_61mers.txt \\
            --kmer_database database_output/wildcard/61_kmer_db.json \\
            --kmer_size 61
        """
        copy_db_cmd = "cp $database/card* database_output;cp -r $database/wildcard/ database_output"
    } else if (!database) { // Download CARD DB if no path given
        db_download_cmd = "wget https://card.mcmaster.ca/latest/data;tar -xvf data ./card.json -C database_output"
        db_card_cmd = "rgi card_annotation -i database_output/card.json"
        database = "."
        copy_db_cmd = "cp $database/card* database_output"
    }

    """
    mkdir database_output
    $db_download_cmd
    $db_card_cmd
    DB_VERSION=\$(ls $database/card_database_*_all.fasta | sed "s/$database\\/card_database_v\\([0-9].*[0-9]\\).*/\\1/")
    $db_wildcard_cmd
    $copy_db_cmd

    rgi \\
        load \\
        $args \\
        --card_json database_output/card.json \\
        --debug --local \\
        --card_annotation database_output/card_database_v\$DB_VERSION.fasta \\
        --card_annotation_all_models database_output/card_database_v\$DB_VERSION\\_all.fasta \\
        $load_wildcard

    rgi \\
        main \\
        $args \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    mkdir temp/
    mv *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json temp/

    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """

stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p temp
    touch test.json
    touch test.txt
    mkdir database_output

    RGI_VERSION=\$(rgi main --version)
    DB_VERSION=stub_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """
}
