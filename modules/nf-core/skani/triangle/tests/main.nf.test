nextflow_process {

    name "Test Process SKANI_TRIANGLE"
    script "../main.nf"
    process "SKANI_TRIANGLE"

    tag "modules"
    tag "modules_nfcore"
    tag "skani"
    tag "skani/triangle"

    // Dependencies
    tag "skani/sketch"

    setup {
        run("SKANI_SKETCH", alias: "SKANI_SKETCH1") {
            script "../../sketch/main.nf"
            process {
                """
                input[0] = [
                    [ id: 'test' ],
                    file(params.modules_testdata_base_path + "genomics/sarscov2/genome/genome.fasta.gz", checkIfExists: true)
                    ]
                """
            }
        }

        run("SKANI_SKETCH", alias: "SKANI_SKETCH2") {
            script "../../sketch/main.nf"
            process {
                """
                input[0] = [
                    [ id: 'test' ],
                    file(params.modules_testdata_base_path + "genomics/sarscov2/illumina/fasta/contigs.fasta", checkIfExists: true)
                    ]
                """
            }
        }
    }



    test("fasta.gz fastq.gz") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test' ], // meta map
                    [
                        file(params.modules_testdata_base_path + "genomics/sarscov2/genome/genome.fasta.gz", checkIfExists: true),
                        file(params.modules_testdata_base_path + "genomics/sarscov2/illumina/fasta/contigs.fasta", checkIfExists: true)
                    ]
                ]
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }

    test("fasta.gz fastq.gz - stub") {

        options "-stub"

        when {
            process {
                """
                input[0] = [
                    [ id:'test' ], // meta map
                    [
                        file(params.modules_testdata_base_path + "genomics/sarscov2/genome/genome.fasta.gz", checkIfExists: true),
                        file(params.modules_testdata_base_path + "genomics/sarscov2/illumina/fasta/contigs.fasta", checkIfExists: true)
                    ]
                ]
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
