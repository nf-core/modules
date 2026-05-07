process HTSLIB_BGZIPTABIX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/53334d0d6ed42b863e50a0feb601505035ed5d5ae9fe6507cadceacc9f1545aa/data':
        'community.wave.seqera.io/library/htslib:1.23.1--45117a0a8dbaa21c' }"

    input:
    tuple val(meta), path(in_file)
    val action
    val make_index
    val out_ext

    output:
    tuple val(meta), path("${outfile}")          , emit: output
    tuple val(meta), path("${outfile}.{tbi,csi}"), emit: index , optional: true
    // all htslib tools have the same version, we use bgzip
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1! d; s/bgzip (htslib) //'"), topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    def allowed_actions = ["compress", "decompress"]
    if (action !in allowed_actions) {
        error "htslib/bgziptabix: Invalid action: ${action}. Allowed actions are: ${allowed_actions.join(', ')}"
    }

    if (action == "decompress" && make_index) {
        log.warn("htslib/bgziptabix: Cannot create index when decompressing. Ignoring make_index option.")
    }

    def args  = task.ext.args   ?: ''
    def args2 = task.ext.args2  ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    outfile   = action == "compress" ? (out_ext ? "${prefix}.${out_ext}.gz" : "${prefix}.gz") : (out_ext ? "${prefix}.${out_ext}" : "${prefix}")

    def compress_cmd = action == "compress" ? "bgzip -c ${args} -@ ${task.cpus}" : "cat"
    def bgzip_cmd    = action == "compress" ? "[ '\$(basename ${in_file})' != '\$(basename ${outfile})' ] && ln -s ${in_file} ${outfile} ;;" : "bgzip -d -c ${args} -@ ${task.cpus} ${in_file} > ${outfile} ;;"
    def tabix_cmd    = make_index ? "tabix -@ ${task.cpus} ${args2} -f ${outfile}" : ""
    """
    FILE_TYPE=\$(htsfile ${in_file})
    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            ${bgzip_cmd}
        *gzip-compressed*)
            [ "\$(basename ${in_file})" == "\$(basename ${outfile})" ] && echo "Input and output names cannot be the same" && exit 1
            zcat  ${in_file} | ${compress_cmd} > ${outfile} ;;
        *bzip2-compressed*)
            bzcat ${in_file} | ${compress_cmd} > ${outfile} ;;
        *XZ-compressed*)
            xzcat ${in_file} | ${compress_cmd} > ${outfile} ;;
        *)
            ${compress_cmd} ${in_file} > ${outfile} ;;
    esac

    ${tabix_cmd}
    """

    stub:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    outfile   = action == "compress" ? (out_ext ? "${prefix}.${out_ext}.gz" : "${prefix}.gz") : (out_ext ? "${prefix}.${out_ext}" : "${prefix}")

    def touch_cmd = action == "compress" ? "echo | bgzip -c" : "echo"
    def tabix_cmd = make_index ? "tabix -@ ${task.cpus} ${args2} -f ${outfile}" : ""
    """
    echo $args

    $touch_cmd > $outfile

    ${tabix_cmd}
    """
}
