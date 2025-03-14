process GRAPHTYPER_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphtyper:2.7.2--h7d7f7ad_0':
        'biocontainers/graphtyper:2.7.2--h7d7f7ad_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_fai)
    path region_file  // can be empty if --region is supplied to task.ext.args

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_path_text = bam.join('\\n')
    def bai_path_text = bam.collect{"${it}.csi"}.join('\\n')
    def region_text = region_file.size() > 0 ? "--region_file ${region_file}" : ""
    if (region_file.size() == 0 && ! args.contains("region")) {
        error "GRAPHTYPER_GENOTYPE requires either a region file or a region specified using '--region' in ext.args"
    }
    """
    # Decompress reference file if needed and set file names
    if [[ $ref =~ \\.gz\$ ]]; then
        gzip -dc $ref > __my__reference__.fasta
    else
        ln -s $ref __my__reference__.fasta
    fi
    ln -s $ref_fai __my__reference__.fasta.fai

    # Make file of file names to pass BAM paths to graphtyper
    printf "$bam_path_text" > bam_list.txt
    printf "$bai_path_text" > bai_list.txt

    # Call graphtyper genotype
    graphtyper \\
        genotype \\
        __my__reference__.fasta \\
        $args \\
        --sams bam_list.txt \\
        --sams_index bai_list.txt \\
        --threads $task.cpus \\
        $region_text

    # Move result files into working directory for output
    find results -maxdepth 2 -name '*.vcf*' > output_paths.txt
    sed 's_results/__g' output_paths.txt | sed 's_/_-_g' > output_names.txt
    paste -d ' ' output_paths.txt output_names.txt | xargs -I {} echo "mv {}" > mv_commands.sh
    source mv_commands.sh

    # Save version information for graphtyper
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > region1.vcf.gz
    echo | gzip > region2.vcf.gz
    touch region1.vcf.gz.tbi
    touch region2.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """

}
