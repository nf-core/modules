process GLIMPSE2_PHASE {
    tag "$meta.id"
    label 'process_medium'

    beforeScript  """
    if cat /proc/cpuinfo | grep avx2 -q
    then
        echo "Feature AVX2 present"
    else
        echo "Feature AVX2 not present on node"
        exit 1
    fi
    """

    conda "bioconda::glimpse-bio=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.0--hf340a29_0':
        'quay.io/biocontainers/glimpse-bio:2.0.0--hf340a29_0' }"

    input:
        tuple val(meta), path(input), path(input_index), val(input_region), val(output_region), path(reference), path(reference_index), path(map), path(samples_file)
        tuple val(meta2), path(fasta_reference), path(fasta_reference_index)
    output:
        tuple val(meta), path("*.{vcf,bcf,bgen}"), emit: phased_variant
        tuple val(meta), path("*.txt.gz")        , emit: stats_coverage, optional: true
        path "versions.yml"                      , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}_${input_region.replace(":","_")}"
    def suffix = task.ext.suffix ?: "bcf"

    def map_command           = map                 ? "--map $map"                    : ""
    def samples_file_command  = samples_file        ? "--samples-file $samples_file"  : ""
    def input_region_command  = input_region        ? "--input-region $input_region"  : ""
    def output_region_command = output_region       ? "--output-region $output_region": ""
    def fasta_command         = fasta_reference     ? "--fasta $fasta_reference"      : ""
    def input_bam             = input.any { it.toString().endsWith(".{bam,cram}")}

    """
    if $input_bam
    then
        printf "%s\\n" *.bam > all_bam.txt
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
        $input_region_command \\
        $output_region_command \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_split_reference --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
