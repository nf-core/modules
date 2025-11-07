process TRIDENT_FETCH {
    tag ""
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/poseidon-trident:1.6.7.1--hebebf5b_1'
        : 'biocontainers/poseidon-trident:1.6.7.1--hebebf5b_1'}"

    input:
    tuple path(archive_dir), val(fetch_s), path(fetch_fn)

    output:
    // All outputs are optional as fetch will check if the package already exists in the provided archive directories, and skip download if so.
    path "*/POSEIDON.yml", emit: poseidon_yml, optional: true
    path "*/*.{bed,geno,vcf}{,.gz}", emit: geno, optional: true
    path "*/*.{bim,snp,vcf}{,.gz}", emit: snp, optional: true
    path "*/*.{fam,ind,vcf}{,.gz}", emit: ind, optional: true
    path "*/*.janno", emit: janno, optional: true
    path "*/*.ssf", emit: ssf, optional: true
    path "*/*.bib", emit: bib, optional: true
    path "*/CHANGELOG.md", emit: changelog, optional: true
    path "*/README.md", emit: readme, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fetch_string = fetch_s ? "--fetchString ${fetch_s}" : ''
    def fetch_file = fetch_fn ? "--fetchFile ${fetch_fn}" : ''
    // fetch will always download to the first directory provided in `-d`, but check all provided dirs for already downloaded packages.
    // Handle multiple archive directories if provided
    def archives = archive_dir ? '-d ' + archive_dir.join(" -d ") : ''
    """
    trident fetch\\
        -d . \\
        ${archives} \\
        ${args} \\
        ${fetch_string} \\
        ${fetch_file}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trident: \$(trident --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def fetch_string = fetch_s ? "--fetchString ${fetch_s}" : ''
    def fetch_file = fetch_fn ? "--fetchFile ${fetch_fn}" : ''
    def archives = archive_dir ? archive_dir.join(" -d ") : ''
    """
    echo ${archives} ${fetch_string} ${fetch_file} ${args}
    
    mkdir dummy_package_dir
    touch dummy_package_dir/POSEIDON.yml
    touch dummy_package_dir/dummy_package.geno
    touch dummy_package_dir/dummy_package.snp
    touch dummy_package_dir/dummy_package.ind
    touch dummy_package_dir/dummy_package.janno

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trident: \$(trident --version)
    END_VERSIONS
    """
}
