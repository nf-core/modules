// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process SCSPLIT_COUNT {
    tag "$meta.id"
    label 'process_high'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/irenerobles93/scsplit:0.0.1"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //    'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), tuple val(sampleId), path(bam), path(bai), path(vcf)
    val tag
    val com
    val ref
    val alt

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def common_data = com != 'None' ? "--com $com" : ''
    def vcf_data = "-v $vcf"
    def bam_data = "-i $bam"
    def barcode_data = "-b $barcode"
    def tag_data = "--tag $tag"

    """
    git clone https://github.com/jon-xu/scSplit
    mkdir scsplit_${prefix}
    mkdir $out
    touch scsplit_${prefix}/params.csv
    echo -e "Argument,Value \n vcf,$vcf \n bam,$bam \n barcode,$barcode \n common_data,${common_data_name} \n num,${num} \n sub,${sub_yesno} \n ems,${ems} \n dbl,${dbl} \n vcf_known_data,${vcf_known_data_name}" >> scsplit_${sampleId}/params.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scsplit: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scsplit: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
