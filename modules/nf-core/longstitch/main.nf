process LONGSTITCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longstitch:1.0.5--hdfd78af_0':
        'biocontainers/longstitch:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(reads)
    val command
    val span
    val genomesize
    val longmap

    output:
    tuple val(meta), path("*.tigmint-ntLink.fa"),                                           emit: tigmint_ntLink_fasta,                 optional: true
    tuple val(meta), path("*.tigmint-ntLink-arks.fa"),                                      emit: tigmint_ntLink_arcs_fasta,            optional: true
    tuple val(meta), path("*.ntLink-arks.fa"),                                              emit: ntLink_arcs_fasta,                    optional: true
    tuple val(meta), path("*k*.w*.z*.n*.scaffold.dot"),                                     emit: scaffold_dot,                         optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*.dist.gv"),                          emit: links_scaffolds_dist_gv,              optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_l*.assembly_correspondence.tsv"),   emit: links_assembly_correspondence_tsv,    optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_l*.gv"),                            emit: links_gv,                             optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_l*.log"),                           emit: links_log,                            optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_l*.scaffolds"),                     emit: links_scaffolds,                      optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_l*.scaffolds.fa"),                  emit: links_scaffolds_fa,                   optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_l*.*.tigpair_checkpoint.tsv"),      emit: links_checkpoint_tsv,                 optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_main.tsv"),                         emit: arcs_tsv,                             optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*_original.gv"),                      emit: arcs_gv,                              optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds_*.tigpair_checkpoint.tsv"),           emit: arcs_checkpoint_tsv,                  optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds.fa"),                                 emit: arcs_scaffolds_fa,                    optional: true
    tuple val(meta), path("*k*.w*.z*.ntLink.scaffolds.renamed.fa"),                         emit: arcs_scaffolds_renamed,               optional: true
    tuple val(meta), path("*k*.w*.z*.stitch.abyss-scaffold.fa"),                            emit: abyss_scaffolds_fa,                   optional: true
    tuple val(meta), path("*k*.w*.z*.stitch.path"),                                         emit: stitch_path,                          optional: true
    tuple val(meta), path("*k*.w*.z*.trimmed_scafs.agp"),                                   emit: trimmed_scaffolds_agp,                optional: true
    tuple val(meta), path("*k*.w*.z*.trimmed_scafs.fa"),                                    emit: trimmed_scaffolds_fasta,              optional: true
    tuple val(meta), path("*k*.w*.z*.trimmed_scafs.path"),                                  emit: trimmed_scaffolds_path,               optional: true
    tuple val(meta), path("*k*.w*.z*.trimmed_scafs.tsv"),                                   emit: trimmed_scaffolds_tsv,                optional: true
    tuple val(meta), path("*k*.w*.z*.verbose_mapping.tsv"),                                 emit: verbose_mapping_tsv,                  optional: true
    tuple val(meta), path("*.k*.w???.tsv"),                                                 emit: tsv,                                  optional: true
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // defs need to happen before ifs; https://github.com/nextflow-io/nextflow/issues/804
    def valid_commands = ["run", "tigmint-ntLink-arks", "tigmint-ntLink", "ntLink-arks"] // run is equivalent to tigmint-ntLink
    def valid_longmaps = [ "ont", "pb", "hifi" ]
    def longmap_val = longmap ? longmap : "ont"
    def arg_longmap = "longmap=${longmap_val}"
    def arg_span = span ? "span=${span}" : ""
    def arg_genomesize = genomesize ? "G=${genomesize}" : ""
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ( !valid_longmaps.contains(longmap_val) ) { error "Unrecognised longmap option. Options: ${valid_longmaps.join(', ')}" }
    if ( !span && !genomesize ) { error "longstitch requires either span or genomesize" }
    if ( !valid_commands.contains(command) )  { error "Unrecognised command to run longstitch. Options: ${valid_commands.join(', ')}" }

    """
    if [[ ${assembly} == *.gz ]]; then
        zcat ${assembly} | fold -w 120 > assembly.fa
    fi

    if [[ ${assembly} == *.fa || ${assembly} == *.fasta ]]; then
        cat ${assembly} | fold -w 120 > assembly.fa
    fi

    if [[ ${reads} == *.fa || ${reads} == *.fasta || ${reads} == *.fq || ${reads} == *.fastq ]]; then
        gzip -c ${reads} > ${reads}.gz
        ln -s ${reads}.gz reads.fq.gz
    fi

    if [[ ${reads} == *.gz ]]; then
        ln -s ${reads} reads.fq.gz
    fi

    longstitch \\
        ${command} \\
        draft=assembly \\
        reads=reads \\
        ${arg_longmap} \\
        ${arg_span} \\
        ${arg_genomesize} \\
        t=$task.cpus \\
        out_prefix=${prefix} \\
        ${args}

    # ["run", "tigmint-ntLink-arks", "tigmint-ntLink", "ntLink-arks"] // run is equivalent to tigmint-ntLink

    if [ ${command} == "run" || ${command} == "tigmint-ntLink" || ${command} == "tigmint-ntLink-arks" ]; then
        mv  *.tigmint-ntLink.longstitch-scaffolds.fa  ${prefix}.tigmint-ntLink.fa
        sed -i 's/\\(scaffold[0-9]*\\),.*/\\1/' ${prefix}.tigmint-ntLink.fa
    fi

    if [ ${command} == "tigmint-ntLink-arks" ]; then
        mv *.tigmint-ntLink-arks.longstitch-scaffolds.fa ${prefix}.tigmint-ntLink-arks.fa
        sed -i 's/\\(scaffold[0-9]*\\),.*/\\1/' ${prefix}.tigmint-ntLink-arks.fa
    fi

    if [ ${command} == "run" || ${command} == "ntLink-arks"]; then
        mv  *.ntLink-arks.longstitch-scaffolds.fa  ${prefix}.ntLink-arks.fa
        sed -i 's/\\(scaffold[0-9]*\\),.*/\\1/' ${prefix}.ntLink-arks.fa
    fi

    # Rename all outputs from tigmint, ntLinks
    for file in assembly.fa.*; do
        mv \$file \${file/assembly.fa./${prefix}.}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longstitch: \$(echo \$(longstitch | head -n1 | sed 's/LongStitch v//'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tigmint-ntLink.fa
    touch ${prefix}.tigmint-ntLink-arks.fa
    touch ${prefix}.ntLink-arks.fa
    touch ${prefix}.k21.w100.z1000.n_stub.scaffold.dot
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub.dist.gv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_l1.assembly_correspondence.tsv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_l1.gv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_l1.log
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_l1.scaffolds
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_l1.scaffolds.fa
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_l1.tigpair_checkpoint.tsv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_main.tsv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub_original.gv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds_stub.tigpair_checkpoint.tsv
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds.fa
    touch ${prefix}.k21.w100.z1000.ntLink.scaffolds.renamed.fa
    touch ${prefix}.k21.w100.z1000.stitch.abyss-scaffold.fa
    touch ${prefix}.k21.w100.z1000.stitch.path
    touch ${prefix}.k21.w100.z1000.trimmed_scafs.agp
    touch ${prefix}.k21.w100.z1000.trimmed_scafs.fa
    touch ${prefix}.k21.w100.z1000.trimmed_scafs.path
    touch ${prefix}.k21.w100.z1000.trimmed_scafs.tsv
    touch ${prefix}.k21.w100.z1000.verbose_mapping.tsv
    touch ${prefix}.k21.w100.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longstitch: \$(echo \$(longstitch | head -n1 | sed 's/LongStitch v//'))
    END_VERSIONS
    """
}
