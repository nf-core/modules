/**
 * Update nf-test snapshots from CI artifacts
 */

/**
 * Find snapshot artifacts from failed nf-test workflow runs
 */
async function findSnapshotArtifacts({ github, context, core }) {
    // React to comment
    await github.rest.reactions.createForIssueComment({
        owner: context.repo.owner,
        repo: context.repo.repo,
        comment_id: context.payload.comment.id,
        content: "eyes",
    });

    // Find latest failed nf-test workflow runs (both regular and GPU)
    const workflows = ["nf-test.yml", "nf-test-gpu.yml"];
    let allFailedRuns = [];

    for (const workflowId of workflows) {
        const runs = await github.rest.actions.listWorkflowRuns({
            owner: context.repo.owner,
            repo: context.repo.repo,
            workflow_id: workflowId,
            status: "completed",
            conclusion: "failure",
            per_page: 10,
        });

        allFailedRuns = allFailedRuns.concat(runs.data.workflow_runs);
    }

    if (allFailedRuns.length === 0) {
        await github.rest.issues.createComment({
            owner: context.repo.owner,
            repo: context.repo.repo,
            issue_number: context.payload.issue.number,
            body: "❌ **No failed nf-test runs found**\n\nI couldn't find any recent failed nf-test or nf-test-gpu workflow runs. Make sure the CI tests have run and failed due to snapshot mismatches first.",
        });
        return;
    }

    // Sort by created_at to get the most recent failed run
    allFailedRuns.sort((a, b) => new Date(b.created_at) - new Date(a.created_at));
    const latestFailedRun = allFailedRuns[0];
    console.log("Found failed run:", latestFailedRun.id);

    // Get artifacts from the failed run
    const artifacts = await github.rest.actions.listWorkflowRunArtifacts({
        owner: context.repo.owner,
        repo: context.repo.repo,
        run_id: latestFailedRun.id,
    });

    const snapshotArtifacts = artifacts.data.artifacts.filter((artifact) =>
        artifact.name.startsWith("updated-snapshots-"),
    );

    if (snapshotArtifacts.length === 0) {
        await github.rest.issues.createComment({
            owner: context.repo.owner,
            repo: context.repo.repo,
            issue_number: context.payload.issue.number,
            body: "🤔 **No snapshot artifacts found**\n\nThe failed CI run didn't upload any snapshot artifacts. This might mean:\n- Tests failed for reasons other than snapshot mismatches\n- The test failures occurred before snapshots could be generated",
        });
        return;
    }

    console.log("Found snapshot artifacts:", snapshotArtifacts.length);
    core.setOutput("run_id", latestFailedRun.id);
    core.setOutput("has_artifacts", "true");
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
