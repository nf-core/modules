nextflow.preview.dsl=2

process MULTIQC {

    // tag "FastQC - $sample_id"

    publishDir "${outdir}/multiqc",
        mode: "copy", overwrite: true

    input:
        path file
        val outdir
        val multiqc_args
        // multiqc_args are best passed into the workflow in the following manner:
        // --multiqc_args="--exlude STAR --title custom_report_title"
        val verbose

    output:
        path "*html",       emit: html

    script:

        if (verbose){
            println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
        }

        """
        multiqc $multiqc_args -x work .
        """

}
