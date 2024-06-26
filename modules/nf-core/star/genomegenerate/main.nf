process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-18542ebf762268f62154251d379b96d001894a57:153528b341042daed3d265b3320b7efb5bf385a8-0' :
        'biocontainers/mulled-v2-18542ebf762268f62154251d379b96d001894a57:153528b341042daed3d265b3320b7efb5bf385a8-0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("star")  , emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args_list   = args.tokenize()
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    if (args_list.contains('--genomeSAindexNbases')) {
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            $include_gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        NUM_BASES=`awk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fai}`

        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            $include_gtf \\
            --runThreadN $task.cpus \\
            --genomeSAindexNbases \$NUM_BASES \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    }

    stub:
    if (gtf) {
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/exonGeTrInfo.tab
        touch star/exonInfo.tab
        touch star/geneInfo.tab
        touch star/genomeParameters.txt
        touch star/sjdbInfo.txt
        touch star/sjdbList.fromGTF.out.tab
        touch star/sjdbList.out.tab
        touch star/transcriptInfo.tab

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/genomeParameters.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    }
}
