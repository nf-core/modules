process GLIMPSE2_PHASE {
    tag "$meta.id"
    label 'process_medium'

    beforeScript  """
    if cat /proc/cpuinfo | grep avx2 -q
    then
        echo "Feature AVX2 present on host"
    else
        echo "Feature AVX2 not present on host"
        exit 1
    fi
    """

    conda "bioconda::glimpse-bio=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.0--hf340a29_0':
        'biocontainers/glimpse-bio:2.0.0--hf340a29_0' }"

    input:
        tuple val(meta) , path(input), path(input_index), path(samples_file), val(input_region), val(output_region), path(reference), path(reference_index), path(map)
        tuple val(meta2), path(fasta_reference), path(fasta_reference_index)

    output:
        tuple val(meta), path("*.{vcf,bcf,bgen}"), emit: phased_variant
        tuple val(meta), path("*.txt.gz")        , emit: stats_coverage, optional: true
        path "versions.yml"                      , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def region = input_region    ? "${output_region.replace(":","_")}" : "${reference}"
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}_${region}"
    def suffix = task.ext.suffix ?: "bcf"

    def map_command           = map                 ? "--map $map"                    : ""
    def samples_file_command  = samples_file        ? "--samples-file $samples_file"  : ""
    def fasta_command         = fasta_reference     ? "--fasta $fasta_reference"      : ""
    def input_region_cmd      = input_region        ? "--input-region $input_region"  : ""
    def output_region_cmd     = output_region       ? "--output-region $output_region": ""

    def input_bam             = input.any { it.extension in ["cram","bam"]}

    """
    if $input_bam
    then
        ls -1 | grep '\\.cram\$\\|\\.bam\$' > all_bam.txt
        input_command="--bam-list all_bam.txt"
    else
        input_command="--input-gl $input"
    fi

    GLIMPSE2_phase \\
        $args \\
        \$input_command \\
        --reference $reference \\
        $map_command \\
        $fasta_command \\
        $samples_file_command \\
        $input_region_cmd \\
        $output_region_cmd \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
