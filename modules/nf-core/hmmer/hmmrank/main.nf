process HMMER_HMMRANK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68261e24307fdf80b9988d6fd13cf3551735e6dc0e7e38003259f7d8efa84cdb/data' :
        'community.wave.seqera.io/library/duckdb-cli:1.5.5--9c6d18d9f687a45d' }"

    input:
    tuple val(meta), path(tblout), path(domtblout)    // Parquet tables from hmmer/formattsv + duckdb/table2parquet; domtblout is optional ([] when absent)

    output:
    tuple val(meta), path("*.hmmrank.tsv.gz"), emit: hmmrank
    tuple val("${task.process}"), val('duckdb'), eval("duckdb --version | sed 's/^v//; s/ .*//'"), topic: versions, emit: versions_duckdb

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A single quote inside a SQL string literal is escaped by doubling it, not by a backslash
    // (that's shell/awk's convention, not SQL's) -- meta.id/task.ext.prefix and file paths are
    // not under this module's control, so any embedded "'" (e.g. an apostrophe in a sample id)
    // would otherwise prematurely close the literal and corrupt the generated SQL.
    // .call(), not a bare sqlLit(...): Nextflow's `script:` block fails to resolve a
    // def-bound closure invoked with direct call syntax ("`sqlLit` is not defined"), even
    // though the exact same closure works fine passed by reference, e.g. to `.collect()`.
    def sqlLit = { s -> s.toString().replace("'", "''") }
    def tbloutSql = sqlLit.call(tblout)
    def domtbloutSql = domtblout ? sqlLit.call(domtblout) : null
    def prefixSql = sqlLit.call(prefix)
    // Reduces one coordinate set's (hmm/ali/env) domain rows to one row per (target, profile,
    // query): the outer bounds of the hit, plus the size and count of the union of its domains.
    // This is the SQL form of a classic "merge overlapping intervals" scan: sorted by each row's
    // own start column, a window MAX(...) tracks how far right everything up to (but not
    // including) the current row already reaches; a start past that opens a new island, and since
    // starts are sorted no later domain can close the gap either. Coordinates are inclusive, so a
    // domain starting exactly one position beyond the previous reach continues it rather than
    // opening a new island -- hence the +1. The two window steps can't be nested in one SQL
    // expression (DuckDB rejects a window function inside another window function's arguments),
    // hence the separate *_prev / *_isl statements below. A hit is keyed on the query name as well
    // as profile, since one HMM file may hold several models and the whole file is then searched
    // in one go, putting several models' domains in the same table.
    // Both windows order by (from, to), not from alone: two domains can share the same start (e.g.
    // repeat domains all starting at hmm position 1), and *_prev/*_isl are two separately
    // materialised temp tables, each free to break that tie differently (parallel threads give no
    // ordering guarantee) -- a total, identical order in both statements is what keeps a tied row's
    // prev_cummax and its running island count referring to the same relative position.
    def islands_sql = { set -> """
CREATE TEMP TABLE ${set}_prev AS
SELECT accno, profile, query, ${set}_from AS f, ${set}_to AS t,
    MAX(${set}_to) OVER (
        PARTITION BY accno, profile, query ORDER BY ${set}_from, ${set}_to
        ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
    ) AS prev_cummax
FROM dom_raw;

CREATE TEMP TABLE ${set}_isl AS
SELECT accno, profile, query, f, t,
    SUM(CASE WHEN f > COALESCE(prev_cummax, 0) + 1 THEN 1 ELSE 0 END) OVER (
        PARTITION BY accno, profile, query ORDER BY f, t
        ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
    ) AS island
FROM ${set}_prev;

CREATE TEMP TABLE ${set}_isl_agg AS
SELECT accno, profile, query, island, MIN(f) AS s, MAX(t) AS e
FROM ${set}_isl
GROUP BY accno, profile, query, island;

CREATE TEMP TABLE ${set}_coords AS
SELECT accno, profile, query,
    MIN(s) AS ${set}_from, MAX(e) AS ${set}_to,
    SUM(e - s + 1) AS ${set}_len, COUNT(*) AS ${set}_n_islands
FROM ${set}_isl_agg
GROUP BY accno, profile, query;
"""
    }
    // dom_raw is read once and reused by every *_prev table and dom_base below, rather than each
    // issuing its own read_parquet() of the same file (4 scans of one input down to 1).
    def domtbl_sql = domtblout ? """
CREATE TEMP TABLE dom_raw AS
SELECT target_name AS accno, profile, query_name AS query, target_length AS tlen, query_length AS qlen,
    hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
FROM read_parquet('${domtbloutSql}');
""" + ['hmm', 'ali', 'env'].collect(islands_sql).join('') + """
CREATE TEMP TABLE dom_base AS
SELECT DISTINCT accno, profile, query, tlen, qlen
FROM dom_raw;

CREATE TEMP TABLE domain_coords AS
SELECT dom_base.accno, dom_base.profile, dom_base.query, dom_base.tlen, dom_base.qlen,
    hmm_coords.hmm_from, hmm_coords.hmm_to, hmm_coords.hmm_len, hmm_coords.hmm_n_islands,
    ali_coords.ali_from, ali_coords.ali_to, ali_coords.ali_len, ali_coords.ali_n_islands,
    env_coords.env_from, env_coords.env_to, env_coords.env_len, env_coords.env_n_islands
FROM dom_base
LEFT JOIN hmm_coords USING (accno, profile, query)
LEFT JOIN ali_coords USING (accno, profile, query)
LEFT JOIN env_coords USING (accno, profile, query);
""" : ''
    // rank 1 is what downstream consumers select on; ties are resolved by profile and then model
    // name (alphabetical) so the order is deterministic even when one file holds several models
    // and score/e-value alone don't break the tie. The final ORDER BY matters for its own sake
    // too, separate from `rank`: DuckDB does not otherwise guarantee row order (query
    // parallelism can interleave rows arbitrarily), so without it the output's row order --
    // and therefore its checksum -- would vary run to run even when its content doesn't.
    def output_select = domtblout ? """
SELECT ranked.profile, ranked.accno, ranked.profile_desc, ranked.evalue, ranked.score, ranked.rank,
    domain_coords.tlen, domain_coords.qlen,
    domain_coords.hmm_from, domain_coords.hmm_to, domain_coords.hmm_len, domain_coords.hmm_n_islands,
    domain_coords.ali_from, domain_coords.ali_to, domain_coords.ali_len, domain_coords.ali_n_islands,
    domain_coords.env_from, domain_coords.env_to, domain_coords.env_len, domain_coords.env_n_islands
FROM ranked
LEFT JOIN domain_coords
    ON ranked.accno = domain_coords.accno AND ranked.profile = domain_coords.profile AND ranked.profile_desc = domain_coords.query
ORDER BY ranked.accno, ranked.rank
""" : 'SELECT * FROM ranked ORDER BY accno, rank\n'

    // Quoted heredoc: no shell expansion in the body -- sqlLit only defends the SQL layer, this
    // defends the shell layer. Delimiter is randomized per task so a real newline in a value
    // can't forge a line that closes it early.
    def heredocTag = "SQL_${java.util.UUID.randomUUID().toString().replace('-', '')}"
    """
    duckdb <<'${heredocTag}'
    SET threads=${task.cpus};
    SET memory_limit='${task.memory.toGiga()}GB';
    SET temp_directory='.';

    CREATE TEMP TABLE ranked AS
    SELECT
        profile, target_name AS accno, query_name AS profile_desc,
        full_evalue AS evalue, full_score AS score,
        ROW_NUMBER() OVER (
            PARTITION BY target_name
            ORDER BY full_score DESC, full_evalue ASC, profile ASC, query_name ASC
        ) AS rank
    FROM read_parquet('${tbloutSql}');
    ${domtbl_sql}
    COPY (${output_select}) TO '${prefixSql}.hmmrank.tsv.gz' (FORMAT CSV, DELIMITER '\\t', HEADER, COMPRESSION 'gzip', NULLSTR 'NA');
${heredocTag}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def coord_columns = domtblout ? '\ttlen\tqlen\thmm_from\thmm_to\thmm_len\thmm_n_islands\tali_from\tali_to\tali_len\tali_n_islands\tenv_from\tenv_to\tenv_len\tenv_n_islands' : ''
    """
    echo 'profile\taccno\tprofile_desc\tevalue\tscore\trank${coord_columns}'  > ${prefix}.hmmrank.tsv
    gzip ${prefix}.hmmrank.tsv
    """
}
