process XTEA {
    tag "XTEA_$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xtea:0.1.9--hdfd78af_0':
        'biocontainers/xtea:0.1.9--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path xtea_resources_dir
    path genecode_gff3
    each xtea_mode

    output:
    tuple val(meta), path("path_work_folder/*/*/*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
 	echo ${prefix} > sample_id.txt
	echo "${prefix} ${bam}" > illumina_bam_list.txt

	mkdir path_work_folder
	
	xtea \\
        $args \\
		 -i sample_id.txt \\
		 -b illumina_bam_list.txt \\
		 -x null \\
		 -p path_work_folder \\
		 -o submit_jobs.sh \\
		 -l ${xtea_resources_dir} \\
		 -r ${fasta} \\
		 -g ${genecode_gff3} \\
		 --xtea /usr/local/bin \\ 
		 -y ${xtea_mode}
	
	sh path_work_folder/*/*/run_xTEA_pipeline.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xtea: \$(xtea --version |& awk '{print $2}' | sed -r 's/v//g)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ''> ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xtea: \$(xtea --version |& awk '{print $2}' | sed -r 's/v//g)
    END_VERSIONS
    """
}
