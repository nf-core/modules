process BCFTOOLS_PLUGINSPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf, stageAs: "input/*"), path(tbi, stageAs: "input/*"), path(samples), path(groups), path(regions), path(targets)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def samples_arg = samples ? "--samples-file ${samples}" : ""
    def groups_arg = groups ? "--groups-file ${groups}" : ""
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
        --output .
    """

    stub:
    def args = task.ext.args ?: ''

    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""
    def determination_file = samples ?: targets
    def create_cmd = extension.matches("vcf|bcf") ? "touch " : "echo '' | gzip > "
    """
    cut -f 3 ${determination_file} | sed -e 's/\$/.${extension}/' > files.txt
    while IFS= read -r filename;
        do ${create_cmd} "./\$filename";
        if [ -n "${index}" ]; then
            index_file=\$(sed -e 's/\$/.${index}/' <<< \$filename);
            touch ./\$index_file;
        fi;
    done < files.txt
    """
}
