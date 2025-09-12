/**
 * Update nf-test snapshots from CI artifacts
 */

/**
 * Find snapshot artifacts from failed nf-test workflow runs for the current PR
 */
async function findSnapshotArtifacts({ github, context, core }) {
  // React to comment
  await github.rest.reactions.createForIssueComment({
    owner: context.repo.owner,
    repo: context.repo.repo,
    comment_id: context.payload.comment.id,
    content: "eyes",
  });

  // Get the PR head commit SHA
  const prNumber = context.payload.issue.number;
  const pr = await github.rest.pulls.get({
    owner: context.repo.owner,
    repo: context.repo.repo,
    pull_number: prNumber,
  });
  const headSha = pr.data.head.sha;
  console.log("PR head commit:", headSha);

  // Find workflow runs for this specific commit
  const allRuns = await github.rest.actions.listWorkflowRunsForRepo({
    owner: context.repo.owner,
    repo: context.repo.repo,
    head_sha: headSha,
    status: "completed",
    per_page: 100,
  });

  // Filter for failed nf-test workflows
  const targetWorkflows = ["Run nf-tests", "Run GPU nf-tests"];
  const failedNfTestRuns = allRuns.data.workflow_runs.filter(
    (run) => targetWorkflows.includes(run.name) && run.conclusion === "failure",
  );

  if (failedNfTestRuns.length === 0) {
    await github.rest.issues.createComment({
      owner: context.repo.owner,
      repo: context.repo.repo,
      issue_number: context.payload.issue.number,
      body: "❌ **No failed nf-test runs found for this PR**\n\nI couldn't find any failed 'Run nf-tests' or 'Run GPU nf-tests' workflow runs for this PR's latest commit. Make sure the CI tests have run and failed due to snapshot mismatches first.",
    });
    return;
  }

  console.log(
    `Found ${failedNfTestRuns.length} failed nf-test runs:`,
    failedNfTestRuns.map((r) => `${r.name} (${r.id})`),
  );

  // Collect artifacts from all failed runs
  let allSnapshotArtifacts = [];
  let allRunIds = [];

  for (const run of failedNfTestRuns) {
    console.log(`Checking artifacts for run ${run.id} (${run.name})`);
    const artifacts = await github.rest.actions.listWorkflowRunArtifacts({
      owner: context.repo.owner,
      repo: context.repo.repo,
      run_id: run.id,
    });

    const snapshotArtifacts = artifacts.data.artifacts.filter((artifact) =>
      artifact.name.startsWith("updated-snapshots-"),
    );

    console.log(`Found ${snapshotArtifacts.length} snapshot artifacts in run ${run.id}`);
    allSnapshotArtifacts = allSnapshotArtifacts.concat(snapshotArtifacts);
    if (snapshotArtifacts.length > 0) {
      allRunIds.push(run.id);
    }
  }

  if (allSnapshotArtifacts.length === 0) {
    await github.rest.issues.createComment({
      owner: context.repo.owner,
      repo: context.repo.repo,
      issue_number: context.payload.issue.number,
      body: "🤔 **No snapshot artifacts found**\n\nThe failed CI runs didn't upload any snapshot artifacts. This might mean:\n- Tests failed for reasons other than snapshot mismatches\n- The test failures occurred before snapshots could be generated\n- The failed jobs didn't run long enough to generate updated snapshots",
    });
    return;
  }

  console.log(`Total snapshot artifacts found: ${allSnapshotArtifacts.length} across ${allRunIds.length} runs`);
  core.setOutput("run_ids", allRunIds.join(","));
  core.setOutput("has_artifacts", "true");
  core.setOutput("artifact_count", allSnapshotArtifacts.length.toString());
}

/**
 * Commit updated snapshots and provide feedback
 */
async function commitSnapshots({ github, context, require }) {
  const { execSync } = require("child_process");

  // Check if there are changes
  try {
    execSync("git diff --quiet", { stdio: "inherit" });
    // No changes
    await github.rest.issues.createComment({
      owner: context.repo.owner,
      repo: context.repo.repo,
      issue_number: context.payload.issue.number,
      body: "🤔 **Artifacts found but no changes**\n\nI found snapshot artifacts but they didn't result in any changes. The snapshots may already be up to date.",
    });
  } catch (error) {
    // There are changes, commit them
    execSync('git config user.email "core@nf-co.re"');
    execSync('git config user.name "nf-core-bot"');
    execSync('git add "**/*.nf.test.snap"');
    execSync('git commit -m "[automated] Update nf-test snapshots from CI artifacts"');
    execSync("git push");

    await github.rest.issues.createComment({
      owner: context.repo.owner,
      repo: context.repo.repo,
      issue_number: context.payload.issue.number,
      body: "✅ **Snapshots updated from CI artifacts!**\n\nI found updated snapshots from the failed CI run and committed them to this PR. Much faster than re-running tests! 🚀",
    });
  }
}

module.exports = {
  findSnapshotArtifacts,
  commitSnapshots,
};
