process FASTQSCREEN_BUILDFROMINDEX_BWA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    tuple val(meta_index), path(bwa_index)
    tuple val(meta_fasta), path(fasta)

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta_fasta.id}"
    def dir = "FastQ_Screen_Genomes"
    def fasta_name = fasta.name
    def fasta_ext = fasta.extension
    """
    mkdir -p ${dir}/${prefix}

    # Copy fasta file
    cp ${fasta} ${dir}/${prefix}/

    # Copy BWA index files
    # cp ${bwa_index}/* ${dir}/${prefix}/

    # Copy and optionally rename BWA index files, since
    # fastq-screen expects the filename of the index files to
    # include the full filename of the fasta (including its extensions)
    # E.g. genome.amb -> genome.fasta.amb
    for index_file in ${bwa_index}/*; do
        if [[ \${index_file} =~ "\\.${fasta_ext}\\." ]]; then
            cp \${index_file} ${dir}/${prefix}
        else
            cp \${index_file} ${dir}/${prefix}/${fasta_name}.\${index_file##*.}
        fi
    done

    # Verify index was created correctly and matches fasta
    if [ ! -f "${dir}/${prefix}/${fasta_name}.amb" ]; then
        echo "ERROR: BWA index filenames (${bwa_index}) do not match fasta name (${fasta_name})."
        exit 1
    fi

    # Generate config
    cat > ${dir}/fastq_screen.conf <<EOF
    DATABASE ${prefix} ${dir}/${prefix}/${fasta.name}
    EOF

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
