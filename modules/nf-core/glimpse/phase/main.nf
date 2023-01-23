process GLIMPSE_PHASE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::glimpse-bio=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--h2ce4488_2':
        'quay.io/biocontainers/glimpse-bio:1.1.1--hce55b13_1' }"

    input:
        tuple val(meta), path(input), path(inp_index)
        path(map)
        path(samples_file)
        tuple val(ref), path(reference), path(ref_index)
        path(region_file)
        val(input_region)
        val(output_region)

    output:
        tuple val(meta), path("*.{vcf,bcf}"), emit: phased_variant
        path "versions.yml"                 , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def map_command = map ? "--map $map" : ""
    def samples_file_command = samples_file ? "--samples-files $samples_file":""
    def input_region_command = input_region ? "--input-region $input_region":""
    def output_region_command = output_region ? "--output-region $output_region":""

    def by_file = region_file ? true : false

    if(by_file) {
        """
        while IFS="" read -r LINE || [ -n "\$LINE" ];
        do
            printf -v ID "%02d" \$(echo \$LINE | cut -d" " -f1)
            IRG=\$(echo \$LINE | cut -d" " -f3)
            ORG=\$(echo \$LINE | cut -d" " -f4)
            OUT=${prefix}_\${ID}.bcf

            GLIMPSE_phase \\
                $args \\
                --input $input \\
                --reference $reference \\
                $map_command \\
                $samples_file_command \\
                --input-region \${IRG} \\
                --output-region \${ORG} \\
                --thread $task.cpus \\
                --output \${OUT}

        done < $region_file

        version_glimpse=\$(GLIMPSE_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')
        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                glimpse: "\$version_glimpse"
        END_VERSIONS
        """
    }
    else {
        """
        GLIMPSE_phase \\
            $args \\
            --input $input \\
            --reference $reference \\
            $map_command \\
            $samples_file_command \\
            $input_region_command \\
            $output_region_command \\
            --thread $task.cpus \\
            --output ${prefix}.bcf

        version_glimpse=\$(GLIMPSE_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')
        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                glimpse: "\$version_glimpse"
        END_VERSIONS
        """
    }
}
