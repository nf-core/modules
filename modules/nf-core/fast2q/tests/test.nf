include { FAST2Q } from './modules/local/fast2q'

workflow {
    Channel
        .fromPath('./tests')  // input is a folder
        .set { fastq_files }

    // Call the FAST2Q process
    FAST2Q(fastq_files)
}
