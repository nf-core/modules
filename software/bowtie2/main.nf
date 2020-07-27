nextflow.preview.dsl=2
params.genome = ''

process BOWTIE2 {
    // depending on the genome used one might want/need to adjust the memory settings.
    // For the E. coli test data this is probably not required

    // label 'bigMem'
    // label 'multiCore'

    publishDir "$outdir/bowtie2",
        mode: "copy", overwrite: true

    input:
        tuple val(name), path(reads)
        val (outdir)
        val (bowtie2_args)
        val (verbose)

    output:
        path "*bam",  emit: bam
        path "*stats.txt", emit: stats

    script:
        if (verbose){
            println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
        }

        cores = 4

        readString = ""

        // Options we add are
        bowtie2_options = bowtie2_args
        bowtie2_options +=  " --no-unal "  // We don't need unaligned reads in the BAM file

        // single-end / paired-end distinction. Might also be handled via params.single_end
        if (reads instanceof List) {
            readString = "-1 " + reads[0] + " -2 " + reads[1]
        }
        else {
            readString = "-U " + reads
        }

        index = params.genome["bowtie2"]
        bowtie2_name = name + "_" + params.genome["name"]

        """
        bowtie2 -x ${index} -p ${cores} ${bowtie2_options} ${readString}  2>${bowtie2_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie2_name}_bowtie2.bam
        """

}
