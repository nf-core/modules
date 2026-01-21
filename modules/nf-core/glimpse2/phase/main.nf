process GLIMPSE2_PHASE {
    tag "${meta.id}"
    label 'process_medium'

    beforeScript """
    if cat /proc/cpuinfo | grep avx2 -q
    then
        echo "Feature AVX2 present on host"
    else
        echo "Feature AVX2 not present on host"
        exit 1
    fi
    """

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.1--h46b9e50_1'
        : 'biocontainers/glimpse-bio:2.0.1--h46b9e50_1'}"

    input:
    tuple val(meta), path(input, arity: '1..*'), path(input_index), path(bamlist), path(samples_file), val(input_region), val(output_region), path(reference), path(reference_index), path(map)
    tuple val(meta2), path(fasta_reference), path(fasta_reference_index)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bgen}"), emit: phased_variants
    tuple val(meta), path("*.txt.gz"), emit: stats_coverage, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def region = input_region ? "${output_region.replace(":", "_")}" : "${reference}"
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}_${region}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    def map_command = map ? "--map ${map}" : ""
    def samples_file_command = samples_file ? "--samples-file ${samples_file}" : ""
    def fasta_command = fasta_reference ? "--fasta ${fasta_reference}" : ""
    def input_region_cmd = input_region ? "--input-region ${input_region}" : ""
    def output_region_cmd = output_region ? "--output-region ${output_region}" : ""

    def input_type = input
        .collect { input_ ->
            input_.toString().endsWithAny("cram", "bam")
                ? "bam"
                : input_.toString().endsWithAny("vcf", "bcf", "vcf.gz")
                    ? "gl"
                    : input_.getExtension()
        }
        .unique()

    if (input_type.size() > 1 | !(input_type.contains("gl") | input_type.contains("bam"))) {
        error("Input files must be of the same type and either .bam/.cram or .vcf/.vcf.gz/.bcf format. Found: ${input_type}")
    }
    else {
        input_type = input_type[0]
    }
    if (input_type == "gl" & input.size() > 1) {
        error("Only one input .vcf/.vcf.gz/.bcf file can be provided")
    }
    def input_list = input.size() > 1

    """
    if [ -n "${bamlist}" ] ;
    then
        input_command="--bam-list ${bamlist}"
    elif ${input_list} ;
    then
        ls -1 | grep '\\.cram\$\\|\\.bam\$' | sort > all_bam.txt
        input_command="--bam-list all_bam.txt"
    else
        if [ "${input_type}" == "bam" ];
        then
            input_command="--bam-file ${input}"
        elif [ "${input_type}" == "gl" ];
        then
            input_command="--input-gl ${input}"
        else
            echo "Input file type not recognised"
            echo "${input_type}"
            exit 1
        fi
    fi

    GLIMPSE2_phase \\
        ${args} \\
        \$input_command \\
        --reference ${reference} \\
        ${map_command} \\
        ${fasta_command} \\
        ${samples_file_command} \\
        ${input_region_cmd} \\
        ${output_region_cmd} \\
        --thread ${task.cpus} \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """

    stub:
    def region = input_region ? "${output_region.replace(":", "_")}" : "${reference}"
    def prefix = task.ext.prefix ?: "${meta.id}_${region}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    def create_cmd = suffix.endsWith(".gz") ? "echo | gzip > ${prefix}.${suffix}" : "touch ${prefix}.${suffix}"
    """
    ${create_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
