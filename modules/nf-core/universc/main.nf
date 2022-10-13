process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
        conda (params.enable_conda ? "hcc::cellranger=3.0.2" : null)
    }
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    "tomkellygenetics/universc:1.2.3"
    //container "${ workflow.containerEngine == 'podman' && !task.ext.singularity_pull_docker_container ?
    //    "tomkellygenetics/universc:1.2.3"
    //container "${ workflow.containerEngine == 'docker' && !task.ext.singularity_pull_docker_container ?
    //    "tomkellygenetics/universc:1.2.3"
    container "tomkellygenetics/universc:1.2.3"

    input:
    tuple val(meta), path(reads)
    path  reference


    output:
    tuple val(meta), path("sample-${meta.id}/outs/*"), emit: outs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_arg = meta.samples.unique().join(",")
    def reference_name = reference.name
    def input_reads = meta.single_end ? "--file $reads" : "-R1 ${reads[0]} -R2 ${reads[1]}"
    """
    bash /universc/launch_universc.sh \\
        --technology '${meta.technology}' \\
        --id 'sample-${meta.id}' \\
        ${input_reads} \\
        --reference ${reference_name} \\
        --sample ${sample_arg} \\
        --jobmode "local" \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        --per-cell-data \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g'
        universc: bash /universc/launch_universc.sh --version | grep version | grep universc  |  's/^.* //g'
    END_VERSIONS
    """


    stub:
    """
    mkdir -p "sample-${meta.id}/outs/"
    touch sample-${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g'
        universc: bash /universc/launch_universc.sh --version | grep version | grep universc | sed 's/^.* //g'
    END_VERSIONS
    """
}

process CELLRANGER_COUNT_OS {
    tag "$meta.id"
    label 'process_medium'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
        conda (params.enable_conda ? "hcc::cellranger=3.0.2" : null)
    }
    container "tomkellygenetics/universc:1.2.3"

    input:
    tuple val(meta), path(reads)
    path  reference
    output:
    tuple val(meta), path("sample-${meta.id}/outs/*"), emit: outs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_arg = meta.samples.unique().join(",")
    def reference_name = reference.name
    """
    rm -rf ${reference_name}/star
    mkdir -p ${reference_name}/star
    unpigz -k ${reference_name}/genes/genes.gtf.gz

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN 1 \\
        --genomeDir ${reference_name}/star \\
        --genomeFastaFiles ${reference_name}/fasta/genome.fa \\
        --genomeSAindexNbases 6 \\
        --genomeChrBinNbits 15 \\
        --genomeSAsparseD 1 \\
        --limitGenomeGenerateRAM 17179869184 \\
        --sjdbGTFfile ${reference_name}/genes/genes.gtf
    sed -i "s/cellranger-cellranger-....../cellranger-3.0.2/g" ${reference_name}/reference.json

    cp ${reference_name}/fasta/genome.fa genome.fa
    cp ${reference_name}/genes/genes.gtf genes.gtf
    rm -rf ${reference_name}

    cellranger mkref \\
        --genome=${reference_name} \\
        --fasta=genome.fa \\
        --genes=genes.gtf

    cellranger \\
        count  \\
        --id='sample-${meta.id}' \\
        --fastqs=. \\
        --transcriptome=${reference_name} \\
        --sample=${sample_arg} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g'
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "sample-${meta.id}/outs/"
    touch sample-${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g'
    END_VERSIONS
    """
}
