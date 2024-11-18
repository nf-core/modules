process BCFTOOLS_PLUGINSPLIT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(samples)
    path(groups)
    path(regions)
    path(targets)

    output:
    tuple val(meta), path("*/*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*/*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*/*.csi")                    , emit: csi, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: ""

    def samples_arg = samples ? "--samples-file ${samples}" : ""
    def groups_arg  = groups  ? "--groups-file ${groups}"   : ""
    def regions_arg = regions ? "--regions-file ${regions}" : ""
    def targets_arg = targets ? "--targets-file ${targets}" : ""

    """
    bcftools plugin split \\
        ${args} \\
        ${vcf} \\
        ${samples_arg} \\
        ${groups_arg} \\
        ${regions_arg} \\
        ${targets_arg} \\
        --output ${prefix}

    if [ -n "${suffix}" ]; then
        for file in ${prefix}/*; do
            # Extract the basename
            base_name=\$(basename "\$file")
            # Extract the part of the basename before the first dot
            name_before_dot="\${base_name%%.*}"
            # Extract the extension
            extension="\${base_name#\${name_before_dot}}"
            # Construct the new name
            new_name="\${name_before_dot}${suffix}\${extension}"
            mv "\$file" "${prefix}/\$new_name"
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: ""

    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                "vcf"
    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi") ? "tbi" :
                args.contains("--write-index=csi") || args.contains("-W=csi") ? "csi" :
                args.contains("--write-index") || args.contains("-W") ? "csi" :
                ""
    def determination_file = samples ?: targets
    def create_cmd = extension.matches("vcf|bcf") ? "touch " : "echo '' | gzip > "
    """
    mkdir -p ${prefix}

    cut -f 3 ${determination_file} | sed -e 's/\$/.${extension}/' > files.txt
    while IFS= read -r filename;
    do ${create_cmd} "${prefix}/\$filename";
    if [ -n "${index}" ]; then
        index_file=\$(sed -e 's/\$/.${index}/' <<< \$filename);
        touch ${prefix}/\$index_file;
    fi;
    done < files.txt

    if [ -n "${suffix}" ]; then
        for file in ${prefix}/*; do
            # Extract the basename
            base_name=\$(basename "\$file")
            # Extract the part of the basename before the first dot
            name_before_dot="\${base_name%%.*}"
            # Extract the extension
            extension="\${base_name#\${name_before_dot}}"
            # Construct the new name
            new_name="\${name_before_dot}${suffix}\${extension}"
            mv "\$file" "${prefix}/\$new_name"
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
