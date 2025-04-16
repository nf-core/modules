process SHAPEIT5_PHASECOMMON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--hb60d31d_0':
        'biocontainers/shapeit5:5.1.1--hb60d31d_0'}"

    input:
        tuple val(meta) , path(input), path(input_index), path(pedigree), val(region)
        tuple val(meta2), path(reference), path(reference_index)
        tuple val(meta3), path(scaffold), path(scaffold_index)
        tuple val(meta4), path(map)

    output:
        tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: phased_variant
        path "versions.yml"                               , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"

    if ("$input" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def map_command       = map       ? "--map $map"             : ""
    def reference_command = reference ? "--reference $reference" : ""
    def scaffold_command  = scaffold  ? "--scaffold $scaffold"   : ""
    def pedigree_command  = pedigree  ? "--pedigree $pedigree"   : ""

    """
    SHAPEIT5_phase_common \\
        $args \\
        --input $input \\
        $map_command \\
        $reference_command \\
        $scaffold_command \\
        $pedigree_command \\
        --region $region \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: "\$(SHAPEIT5_phase_common | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """

    stub:
    def prefix     = task.ext.prefix        ?: "${meta.id}"
    def suffix     = task.ext.suffix        ?: "bcf"
    def create_cmd = suffix.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: "\$(SHAPEIT5_phase_common | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
