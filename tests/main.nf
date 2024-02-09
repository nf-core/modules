params.module_dir = "./modules/nf-core/samtools/view"
params.repo = "../../modules/"

workflow {
    NFTEST_MODULE (
        file(params.module_dir + "/main.nf"),
        file(params.module_dir + "/tests"),
        [], // TODO
        false
    )
}

process NFTEST_MODULE {
    input:
    path module_file
    path test_file
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
        --verbose
    """
}
