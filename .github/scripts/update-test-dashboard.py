#!/usr/bin/env python3
"""
Update the monthly test dashboard issue with batch results.
"""

import json
import os
import sys
from datetime import datetime, timezone
from typing import Optional
from urllib.error import URLError
from urllib.request import Request, urlopen

TABLE_HEADER = "| Batch | Status |"
TABLE_SEPARATOR = "|-------|--------|"
FAILED_SECTION_HEADER = "## Failed Tests"


def github_api_request(method: str, endpoint: str, data: Optional[dict] = None) -> dict:
    token = os.environ.get("GITHUB_TOKEN")
    repo = os.environ.get("GITHUB_REPOSITORY")

    url = f"https://api.github.com/repos/{repo}/{endpoint}"
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github.v3+json",
        "Content-Type": "application/json",
    }

    req = Request(url, headers=headers, method=method)
    if data:
        req.data = json.dumps(data).encode()

    try:
        with urlopen(req) as response:
            return json.loads(response.read().decode())
    except URLError as e:
        print(f"API request failed: {e}", file=sys.stderr)
        return {}


def get_existing_issue() -> Optional[tuple[int, str]]:
    """Return (number, body) of the existing dashboard issue, or None."""
    issues = github_api_request("GET", "issues?labels=test-dashboard&state=open&per_page=1")
    if not isinstance(issues, list) or not issues:
        return None
    issue = issues[0]
    return issue["number"], issue.get("body", "")


def create_batch_row(
    batch_num: str, batch_label: str, status: str, passed: int, total: int, run_number: str, run_url: str
) -> str:
    status_emoji = {"passed": "✅", "failed": "❌", "cancelled": "⚠️", "skipped": "⏭️"}.get(status, "❓")
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return f"| Batch {batch_num} ({batch_label}) | {status_emoji} {status} | {passed}/{total} | {timestamp} | [Run #{run_number}]({run_url}) |"


def create_failed_tests_section(batch_num: str, failed_count: int, failed_tests: list[str]) -> str:
    if not failed_tests or failed_count == 0:
        return ""

    failed_list = "\n".join(failed_tests)
    return f"""
<details>
<summary>❌ {failed_count} failing test(s) in batch {batch_num}</summary>

```
{failed_list}
```
</details>
"""


def update_body(
    current_body: str,
    batch_num: str,
    batch_row: str,
    failed_section: str,
) -> str:
    lines = current_body.split("\n")
    updated_lines = []
    in_table = False
    batch_updated = False
    in_failed_section = False
    skip_until_detail_close = False

    for line in lines:
        if TABLE_HEADER in line:
            in_table = True
            updated_lines.append(line)
            continue

        if in_table and line.strip() and not line.startswith("|"):
            in_table = False

        if in_table and line.strip().startswith(f"| Batch {batch_num}"):
            updated_lines.append(batch_row)
            batch_updated = True
            continue

        if FAILED_SECTION_HEADER in line:
            in_failed_section = True
            updated_lines.append(line)
            continue

        if in_failed_section and f"batch {batch_num}" in line:
            skip_until_detail_close = True
            continue

        if skip_until_detail_close:
            if "</details>" in line:
                skip_until_detail_close = False
            continue

        updated_lines.append(line)

    if not batch_updated:
        for i, line in enumerate(updated_lines):
            if TABLE_SEPARATOR in line:
                updated_lines.insert(i + 1, batch_row)
                break

    if failed_section:
        for i, line in enumerate(updated_lines):
            if FAILED_SECTION_HEADER in line:
                updated_lines.insert(i + 2, failed_section)
                break

    return "\n".join(updated_lines)


def create_issue_body(batch_row: str, failed_section: str, repo: str, workflow_url: str) -> str:
    failed_text = failed_section if failed_section else "_No failures currently tracked_"

    return f"""# 🧪 Monthly Module Test Dashboard

This issue tracks the status of monthly full module testing across all batches.

## Current Month Status

| Batch | Status | Tests Passed | Last Run | Workflow |
|-------|--------|--------------|----------|----------|
{batch_row}

## Failed Tests

{failed_text}

---
*This dashboard is automatically updated by the [Monthly Full Test workflow]({workflow_url})*
"""


def main():
    batch_num = os.environ["BATCH_NUM"]
    batch_label = os.environ["BATCH_LABEL"]
    status = os.environ["STATUS"]
    run_url = os.environ["RUN_URL"]
    run_number = os.environ["RUN_NUMBER"]
    failed_count = int(os.environ["FAILED_COUNT"])
    passed_count = int(os.environ["PASSED_COUNT"])
    total_count = int(os.environ["TOTAL_COUNT"])
    repo = os.environ["GITHUB_REPOSITORY"]
    workflow_url = f"https://github.com/{repo}/actions/workflows/monthly-full-test.yml"

    failed_tests = []
    try:
        with open("failed_tests.txt") as f:
            failed_tests = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        pass

    batch_row = create_batch_row(batch_num, batch_label, status, passed_count, total_count, run_number, run_url)
    failed_section = create_failed_tests_section(batch_num, failed_count, failed_tests)

    existing = get_existing_issue()

    if existing is None:
        body = create_issue_body(batch_row, failed_section, repo, workflow_url)
        result = github_api_request(
            "POST",
            "issues",
            {"title": "🧪 Monthly Module Test Dashboard", "body": body, "labels": ["test-dashboard"]},
        )
        if result:
            print("✅ Created new dashboard issue")
        else:
            print("❌ Failed to create issue", file=sys.stderr)
            sys.exit(1)
    else:
        issue_number, current_body = existing
        updated_body = update_body(current_body, batch_num, batch_row, failed_section)
        result = github_api_request("PATCH", f"issues/{issue_number}", {"body": updated_body})
        if result:
            print(f"✅ Updated dashboard issue #{issue_number}")
        else:
            print(f"❌ Failed to update issue #{issue_number}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
