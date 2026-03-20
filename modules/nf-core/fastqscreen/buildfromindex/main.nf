process FASTQSCREEN_BUILDFROMINDEX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc53eee7ca23c32220a9662fbb63c67769756544b6d74a1ee85cf439ea79a7ee/data'
        : 'community.wave.seqera.io/library/fastq-screen_perl-gdgraph:5c1786a5d5bc1309'}"

    input:
    val genome_names
    path (indexes), stageAs: "dir*"

    output:
    path ("FastQ_Screen_Genomes"), emit: database
    tuple val("${task.process}"), val('fastqscreen'), eval('fastq_screen --version 2>&1 | sed "s/^.*FastQ Screen v//;"'), emit: versions_fastqscreen, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    dir = "FastQ_Screen_Genomes"
    folder = indexes.collect { index -> index.toString() }
    database = [genome_names, folder].transpose()
    copy_indexes = folder.collect { index -> "cp -r ${index} ${dir}/${index}" }.join(" && ")

    // Folder name and index (within folder) name could be different - use bash to look for index name
    config = database
        .collect { index -> "########## ${index[0]} \nDATABASE ${index[0]} ${dir}/${index[1]}/${index[1] + '_to_be_replaced'}" }
        .join("\n\n")
        .replace("\n", "\\n")

    """
    mkdir ${dir}
    ${copy_indexes}

    echo "${config}" > fastq_screen.conf
    sed 's/\\\\n/\\n/g' fastq_screen.conf > ${dir}/fastq_screen.conf
    echo "Replace folder name real index name"

    for f in ${folder.join(' ')}
    do
        TO_REPLACE=\${f}_to_be_replaced

        REPLACE_WITH=(\$f"/*.rev.1.bt2")
        REPLACE_WITH=\$(basename \$REPLACE_WITH)
        REPLACE_WITH=\${REPLACE_WITH%%.rev.1.bt2}

        sed "s/\$TO_REPLACE/\$REPLACE_WITH/g" ${dir}/fastq_screen.conf > new_conf.conf
        mv new_conf.conf ${dir}/fastq_screen.conf
    done
    """

    stub:
    dir = "FastQ_Screen_Genomes"
    """
    mkdir ${dir}
    touch ${dir}/fastq_screen.conf
    """
}
