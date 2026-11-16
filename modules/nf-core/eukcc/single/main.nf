process EUKCC_SINGLE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4644a12bcaad4ce5d083946976cbd6e0cfef498beeb183e9a9cfda149db49e7/data'
        : 'community.wave.seqera.io/library/epa-ng_eukcc_hmmer_metaeuk_pruned:537c5825af8f6f76'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(db)
    path(genomes, stageAs: 'genomes/*')

    output:
    tuple val(meta), path("${prefix}/eukcc.tsv")                , emit: tsv
    tuple val(meta), path("${prefix}/eukcc.log")                , emit: log
    tuple val(meta), path("${prefix}/scmg_marker_table.csv.gz") , emit: marker_table, optional: true
    tuple val(meta), path("${prefix}/workdir")                  , emit: workdir, optional: true
    tuple val("${task.process}"), val('eukcc'), eval("eukcc --version | sed 's/^EukCC version //'"), topic: versions, emit: versions_eukcc
    tuple val("${task.process}"), val('metaeuk'), eval("metaeuk version"), topic: versions, emit: versions_metaeuk
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -o '^# HMMER [0-9.]*' | cut -d' ' -f3"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('epa-ng'), eval("epa-ng --version | sed 's/^EPA-ng v//'"), topic: versions, emit: versions_epang
    tuple val("${task.process}"), val('pplacer'), eval("pplacer --version | sed 's/^v//; s/-.*//'"), topic: versions, emit: versions_pplacer
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | head -1 | sed 's/^gzip //'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def genomes_arg = genomes ? [genomes].flatten().collect { genome -> "--genomes \"${genome}\"" }.join(' ') : ''
    """
    eukcc \\
        single \\
        ${args} \\
        ${genomes_arg} \\
        --db "${db}" \\
        --threads ${task.cpus} \\
        --out "${prefix}" \\
        "${fasta}"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}"
    printf "fasta\\tcompleteness\\tcontamination\\tncbi_lng\\tncbi_lng_names\\n" > "${prefix}/eukcc.tsv"
    touch "${prefix}/eukcc.log"
    ${args.contains('--extra') ? "echo \"\" | gzip > \"${prefix}/scmg_marker_table.csv.gz\"" : ''}
    ${args.contains('--keep') ? "mkdir \"${prefix}/workdir\"" : ''}
    """
}
