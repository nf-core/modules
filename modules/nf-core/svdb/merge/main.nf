process SVDB_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff995c756aa8a3c0af13b1d054eacd536a11d35de5fa288dacf558bc21696968/data':
        'community.wave.seqera.io/library/bcftools_svdb:ae3b14d2d608fd81' }"

    input:
    tuple val(meta), path(vcfs)
    val(input_priority)
    val(sort_inputs)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*.csi")                    , emit: csi, optional: true
    tuple val("${task.process}"), val('svdb'), eval("svdb | sed -nE 's/.*SVDB-([0-9.]+).*/\\1/p'"), emit: versions_svdb, topic: versions
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Ensure priority list matches the number of VCFs if priority is provided
    if (input_priority && vcfs.collect().size() != input_priority.collect().size()) {
        error "If priority is used, one tag per VCF is needed"
    }

    def input = ""
    def prio = ""
    if (input_priority) {
        if (vcfs.collect().size() > 1 && sort_inputs) {
            // make vcf-priority pairs and sort on VCF name, so priority is also sorted the same
            def pairs = vcfs.indices.collect { index -> [vcfs[index], input_priority[index]] }
            pairs = pairs.sort { a, b -> a[0].name <=> b[0].name }
            vcfs = pairs.collect { vcf -> vcf[0] }
            priority = pairs.collect { pair -> pair[1] }
        } else {
            priority = input_priority
        }

        // Build inputs
        prio = "--priority ${input_priority.join(',')}"
        input = vcfs
            .withIndex()
            .collect { vcf, index -> "${vcf}:${priority[index]}" }
            .join(" ")

    } else {
        // if there's no priority input just sort the vcfs by name if possible
        input = (vcfs.collect().size() > 1 && sort_inputs) ? vcfs.sort { vcf_file -> vcf_file.name } : vcfs
    }

    def extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf.gz"
    """
    svdb \\
        --merge \\
        $args \\
        $prio \\
        --vcf $input |\\
        bcftools view \\
            $args2 \\
            --threads ${task.cpus} \\
            --output ${prefix}.${extension}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2   = task.ext.args2 ?: ''
    def extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf.gz"
    def index_type = args2.contains("--write-index=tbi") || args2.contains("-W=tbi") ? "tbi" :
                args2.contains("--write-index=csi") || args2.contains("-W=csi") ? "csi" :
                args2.contains("--write-index") || args2.contains("-W") ? "csi" :
                ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index_type.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index_type}" : ""
    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    """
}
