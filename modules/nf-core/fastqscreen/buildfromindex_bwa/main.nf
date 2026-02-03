process FASTQSCREEN_BUILDFROMINDEX_BWA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    path(bwa_indices, stageAs: "index_?/*")  // List of [index_dir] from BWA_INDEX.out.index.map{ meta, dir -> dir }.collect()

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_dir = "FastQ_Screen_Genomes"
    """
    mkdir -p ${index_dir}

    # Copy all index files to the same directory
    for idx_dir in ${bwa_indices}; do
        cp -r "\$idx_dir/"* ${index_dir}/
    done

    # Build config by scanning copied directories
    echo "# FastQ Screen Configuration" > ${index_dir}/fastq_screen.conf
    for genome in ${index_dir}/*.amb; do
        genome_name=\$(basename "\$genome")
        echo "DATABASE\t\${genome_name}\t${index_dir}/\${genome_name%.amb}" >> ${index_dir}/fastq_screen.conf
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    index_dir = "FastQ_Screen_Genomes"
    """
    mkdir ${index_dir}
    touch ${index_dir}/fastq_screen.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """
}
