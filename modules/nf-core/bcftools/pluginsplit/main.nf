process BCFTOOLS_PLUGINSPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/914b9b38bbfa22c458ceee6c578fd3d70ee1bffc33d3a064b0663ed4944e7872/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.25--130aa7603bf3bf5c'}"

    input:
    tuple val(meta), path(vcf, stageAs: "input/*"), path(tbi, stageAs: "input/*"), path(samples), path(groups), path(regions), path(targets)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}"), emit: index, optional: true
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
