/**
 * Check for nf-test snapshot failures and post helpful comment
 */
async function checkSnapshotFailures({ github, context, core }) {
    // Only run on pull requests
    if (!context.payload.pull_request) {
        console.log("Not a pull request, skipping snapshot failure check");
        return;
    }

    const prNumber = context.payload.pull_request.number;

    // Look for existing bot comments about snapshot failures
    const comments = await github.rest.issues.listComments({
        owner: context.repo.owner,
        repo: context.repo.repo,
        issue_number: prNumber,
    });

    const botComments = comments.data.filter(
        (comment) => comment.user.login === "nf-core-bot" && comment.body.includes("nf-test failures detected"),
    );

    const commentBody = `⚠️ **nf-test failures detected**

The nf-test workflow has failed. If the failures are due to snapshot mismatches (expected outputs don't match actual outputs), you can automatically update them by commenting:

\`\`\`
@nf-core-bot update snapshots
\`\`\`

This will:
- Download any updated snapshots from this failed CI run
- Commit them directly to your PR
- Much faster than re-running tests locally! 🚀

**Note:** This only works if the failures are due to snapshot mismatches. For other test failures, you'll need to fix the underlying issues first.`;

    if (botComments.length > 0) {
        // Update the existing comment
        await github.rest.issues.updateComment({
            owner: context.repo.owner,
            repo: context.repo.repo,
            comment_id: botComments[0].id,
            body: commentBody,
        });
        console.log("Updated existing snapshot failure comment");
    } else {
        // Create a new comment
        await github.rest.issues.createComment({
            owner: context.repo.owner,
            repo: context.repo.repo,
            issue_number: prNumber,
            body: commentBody,
        });
        console.log("Created new snapshot failure comment");
    }
}

module.exports = { checkSnapshotFailures };
