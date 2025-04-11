// TODO What if we took nf-test list?
params.module_dir = "./modules/nf-core/samtools/view"
params.repo = "../modules/"
params.nftest_config = "./tests/meta/nf-test.config"
params.nextflow_config = "./tests/meta/nf.config"


workflow {
    NFTEST_MODULE (
        file(params.module_dir + "/main.nf"),
        file(params.module_dir + "/tests"),
        file(params.nftest_config),
        file(params.nextflow_config),
        [], // TODO
        false
    )
}

process NFTEST_MODULE {
    tag "$module_file"

    input:
    path module_file
    path test_file
    path nftest_config
    path nextflow_config
    path full_repo // optional for subworkflow
    val update_snapshot

    output:
    path "*.snap", emit: snapshot, optional: true

    script:
    def snapshot = update_snapshot ? '--update-snapshot': ''
    """
    # TODO cd $full_repo
    nf-test test tests/*.nf.test \\
        --profile docker \\
        $snapshot \\
        --silent \\
        --verbose \\
        -c $nftest_config
    """
}
