---
name: nf-core-seqera-containers
description: Move a batch of nf-core/modules modules from quay.io biocontainers to Seqera Containers and open a draft PR. Use when asked to migrate, convert or port modules to Seqera Containers or Wave, or to submit the next batch of quay.io modules (nf-core/modules#13143).
---

# Seqera Containers batch PR

Moves a small set of modules off quay.io onto Seqera Containers using nf-core tools, tests them, and opens one draft PR on nf-core/modules. Run it from a clone of nf-core/modules. Follow the repo's `AGENTS.md` throughout.

`convert.sh` sits next to this file.

## 1. Pick the batch

- Keep it under 12 modules. Bigger PRs overload GitHub Actions and take longer to fix.
- Prefer single-tool modules with few md5 checksums in their snapshots.
- Check open PRs by the files they change, not by title. A title search misses PRs that touch a module without naming it:

  ```bash
  gh pr list -R nf-core/modules --state open --limit 1000 --json number,title,author,isDraft,updatedAt,files \
    | jq -r --arg m "modules/nf-core/<module>/" '.[] | select(any(.files[]; .path | startswith($m))) | "#\(.number) \(.title) [\([.files[].path | select(startswith($m))] | join(", "))]"'
  ```

- Drop any module that another PR already migrates or changes in `main.nf`. If a PR touches only other files in the module, such as `meta.yml`, keep the module but mention the overlap, because the migration also edits `meta.yml`.
- Tell the user which modules you dropped and why before converting.

## 2. Set up

```bash
git fetch upstream master && git switch -c <branch> upstream/master
uv tool install --force "nf-core @ git+https://github.com/nf-core/tools@dev"
nf-core --version   # expect a .dev0 version
```

`TOWER_ACCESS_TOKEN` must be set, or Wave rate-limits the builds.

## 3. Convert

```bash
.claude/skills/nf-core-seqera-containers/convert.sh <module> <module> ...
```

For each module the script runs `nf-core modules bump-versions`, then `nf-core modules containers create` if the module is still on quay.io. It runs them one after another, because parallel runs fail on nf-core's shared git cache. Logs go to `$TMPDIR/seqera-containers/`.

The script leaves you with:

- the new container lines in `main.nf`
- a `containers:` block in `meta.yml`
- `.conda-lock/` files
- updated `environment.yml` and snapshot version strings for any module that was bumped

Check every version bump before keeping it:

- Compare the old and new version in `environment.yml`. `bump-versions` can pick an older release (foldseek 10.941cd33 went to 8.ef4e960). Revert any bump that isn't a newer version.
- Check the version the tool reports in the snapshot. Some binaries report a different version from their package (graphmap 0.6.4 reports 0.6.5).

If Wave can't build a module, the automatic bumps may be the cause. Before dropping it:

1. Restore the original pins in `environment.yml` and build again with no bumps.
2. If that fails, run a local conda dry-run solve to find the dependency that conflicts, and bump only that one to the lowest version that solves. Record what you tried.

A container that builds can still fail when it runs (summarizedexperiment 1.32.0 failed with `no package called 'GenomeInfoDbData'`). If the only fix is a tool version bump, ask the user before keeping it.

If nothing works, drop the module:

```bash
git checkout -- modules/nf-core/<module>
rm -rf modules/nf-core/<module>/.conda-lock
```

An arm64 build can fail because bioconda has no arm64 package. That is fine if the old quay.io image was amd64-only too.

## 4. Lint

```bash
for m in <modules>; do nf-core modules lint $m; done
nextflow lint modules/nf-core/{<modules>}
git add -A modules/ && uvx pre-commit run --files $(git diff --cached --name-only)
```

Fix any new failure. Leave warnings that already existed on master, such as the deprecated `when` section or non-standard labels. `AGENTS.md` forbids editing `when`.

## 5. Test

Run the module tests, then the tests of every subworkflow that includes a converted module:

```bash
nf-test test modules/nf-core/<module>/tests --profile=+docker
grep -rlE "nf-core/(<m1>|<m2>)/main" subworkflows --include='main.nf'
nf-test test subworkflows/nf-core/<subworkflow>/tests --profile=+docker
```

If a subworkflow fails before any process runs, run the same test on a clean copy of master (`git archive upstream/master | tar -x -C <scratch>`). If master fails the same way, the failure isn't caused by the batch.

### Snapshot changes

Local runs on an arm64 Mac use the amd64 image under Rosetta, which is not the same as CI. Tools that write their own command line into the output, like featureCounts and gffread, record the full binary path under Rosetta but the bare command name on CI. A snapshot regenerated locally can then pass here and fail on every CI profile.

When a module's md5 changes:

1. Look at the output file. If only a header line differs (a version or path), the data is unchanged.
2. Ask the user whether to keep the module or drop it. The default rule for these batches is to drop any module whose snapshot changes.
3. If they keep it, prefer making the assertion independent of the path over copying md5s from CI. For a BAM, snapshot the header with `@PG` lines removed, plus the records. `getReadsMD5()` alone failed lint on ivar/trim. Changing what a test checks needs the user's OK, and the PR should say why.
4. Otherwise, push and let CI fail once. Then take the md5s from the CI nf-test log, check that docker, singularity and conda all agree, and put those values in the snapshot. Say in the PR that the md5s came from CI.

R tools in Wave images run with no `LANG` set, so R prints plain-ASCII stand-ins for accented names and quotes. The conda CI runner uses UTF-8, so a log md5 differs on conda only. Snapshot only the lines that are plain ASCII.

Read the md5s from the CI logs:

```bash
gh api --allow-escape-sequences repos/nf-core/modules/actions/jobs/<job_id>/logs \
  | sed 's/\x1b\[[0-9;]*m//g' | grep -E '^\S+\s+[-+]\s+".*:md5,'
```

## 6. Open the PR

Ask before pushing, unless the user or a coordinating session has already approved the push. End the commit message with one "Generated by <model name>" line and no Co-Authored-By line. Then:

```bash
git push -u origin <branch>
gh pr create -R nf-core/modules --base master --head <fork-owner>:<branch> --draft \
  --title "Move <n> modules from quay.io to Seqera Containers" --body-file <body.md>
```

Your agent's permission checks may block `gh pr create`, `gh pr ready` or `gh pr merge --auto`. Don't retry or work around it. Keep the body in a file, stop, and give the user the exact command to run.

`AGENTS.md` only allows draft PRs. Use `.github/pull_request_template.md` for the body. Include:

- the module list, with any version bumps
- why each changed snapshot changed
- the modules you dropped and why
- `Refs #13143`
- "Generated by <model name>" above the checklist

Tick only the test profiles you actually ran. Then bind the PR and offer Auto-fix.

If a CI job fails on a network timeout pulling an image for a module outside the batch, re-run it instead of migrating that module too. Once the run has finished:

```bash
gh run rerun <run_id> -R nf-core/modules --failed
```

Adding a module to the batch needs the user's OK.
