process FASTQSCREEN_BUILDFROMINDEX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    val(genome_names)
    path(indexes), stageAs: "dir*"

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    dir = "FastQ_Screen_Genomes"
    folder = indexes.collect { it.toString() }
    database = [genome_names, folder].transpose()
    copy_indexes = folder.collect { "cp -r ${it} $dir/${it}"}.join(" && ")

    // Folder name and index (within folder) name could be different - use bash to look for index name
    config = database
        .collect { "########## ${it[0]} \nDATABASE ${it[0]} $dir/${it[1]}/${it[1] + '_to_be_replaced'}" }
        .join("\n\n")
        .replace("\n", "\\n")

    """
    mkdir $dir
    $copy_indexes

    echo "$config" > fastq_screen.conf
    sed 's/\\\\n/\\n/g' fastq_screen.conf > $dir/fastq_screen.conf
    echo "Replace folder name real index name"

    for f in ${folder.join(' ')}
    do
        TO_REPLACE=\${f}_to_be_replaced

        REPLACE_WITH=(\$f"/*.rev.1.bt2")
        REPLACE_WITH=\$(basename \$REPLACE_WITH)
        REPLACE_WITH=\${REPLACE_WITH%%.rev.1.bt2}

        sed "s/\$TO_REPLACE/\$REPLACE_WITH/g" $dir/fastq_screen.conf > new_conf.conf
        mv new_conf.conf $dir/fastq_screen.conf
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    dir = "FastQ_Screen_Genomes"
    """
    mkdir $dir
    touch $dir/fastq_screen.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """
}
