process CUSTOM_RESOLVETAXONOMY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(sequences), val(taxonomy_required)

    output:
    // The fallback to 'fasta' matters when sequences has no extension at all (e.g. a
    // download URL with no filename suffix) -- Path.extension is '' there, which
    // would otherwise produce a malformed "*.resolved." glob matching nothing.
    tuple val(meta), path("*.resolved.tax"),                                emit: taxonomy
    tuple val(meta), path("*.resolved.${sequences.extension ?: 'fasta'}"),  emit: sequences
    tuple val(meta), path("*.warnings.txt"),                                emit: warnings
    tuple val("${task.process}"), val('biopython'), eval("python3 -c 'import Bio; print(Bio.__version__)'"), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Falls back to 'fasta' when sequences has no extension at all (e.g. a
    // download URL with no filename suffix) -- Path.extension is '' there, which
    // would otherwise produce a malformed "${prefix}.resolved." output name.
    def ext = sequences.extension ?: 'fasta'
    // Nextflow stages an absent optional path(taxonomy) as an empty list -- falsy in
    // Groovy -- rather than as a file, so this correctly distinguishes "no taxonomy
    // file given" from a real one, without ever interpolating the literal text "[]"
    // into the command line.
    def taxonomy_in = taxonomy ? "${taxonomy}" : ''
    """
    python3 - "${taxonomy_in}" "${sequences}" "${taxonomy_required}" "${prefix}.resolved.tax" "${prefix}.resolved.${ext}" "${prefix}.warnings.txt" << 'PYEOF'
import sys
from Bio import SeqIO

taxonomy_in, sequences_in, taxonomy_required, taxonomy_out, sequences_out, warnings_out = sys.argv[1:7]
taxonomy_required = taxonomy_required == 'true'

# Format sniffed from content (FASTA, Clustal or PHYLIP -- sequences can be any of
# these, depending on the caller).
with open(sequences_in) as fh:
    first_line = next((l.strip() for l in fh if l.strip()), '')
if first_line.startswith('>'):
    sequences_format = 'fasta'
elif first_line.upper().startswith('CLUSTAL'):
    sequences_format = 'clustal'
else:
    sequences_format = 'phylip-relaxed'

records = list(SeqIO.parse(sequences_in, sequences_format))
warnings = []

def embedded_taxonomy(record):
    # record.description is the *whole* header line (id + any trailing text);
    # record.id is just its first token -- GTDB's own single-file convention puts
    # the taxonomy string right after the id, space-separated.
    return record.description[len(record.id):].strip()

if taxonomy_in:
    # An explicit taxonomy file always wins. Warn (not fail) rather than silently
    # using the wrong source if the sequences also happen to carry embedded text.
    if any(embedded_taxonomy(record) for record in records):
        warnings.append(
            'An explicit taxonomy file was provided; ignoring embedded taxonomy '
            'text found in sequence record headers.'
        )
    with open(taxonomy_in) as fh_in, open(taxonomy_out, 'w') as fh_out:
        fh_out.write(fh_in.read())
else:
    missing = [record.id for record in records if not embedded_taxonomy(record)]
    if missing and taxonomy_required:
        sys.exit(
            'No taxonomy file was provided, and these records have no embedded '
            'taxonomy in their header either: ' + ', '.join(missing)
        )
    if missing:
        warnings.append(
            'No taxonomy file was provided; records with no embedded taxonomy in '
            'their header are omitted from the resolved taxonomy: ' + ', '.join(missing)
        )
    with open(taxonomy_out, 'w') as fh:
        for record in records:
            tax = embedded_taxonomy(record)
            if tax:
                print(f"{record.id}\\t{tax}", file=fh)

# Always strip headers down to a bare id -- some downstream tools keep the whole
# header line as the sequence/leaf name rather than just the first token.
for record in records:
    record.description = record.id
SeqIO.write(records, sequences_out, sequences_format)

with open(warnings_out, 'w') as fh:
    for warning in warnings:
        print(warning, file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = sequences.extension ?: 'fasta'
    """
    touch ${prefix}.resolved.tax ${prefix}.resolved.${ext} ${prefix}.warnings.txt
    """
}
