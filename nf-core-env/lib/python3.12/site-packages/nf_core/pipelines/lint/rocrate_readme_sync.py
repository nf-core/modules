import json
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def rocrate_readme_sync(self):
    """
    Check if the RO-Crate description in ro-crate-metadata.json matches the README.md content.
    If not, the RO-Crate description will be automatically updated to match the README.md content during linting.
    """

    passed = []
    ignored = []
    fixed = []

    # Check if the file exists before trying to load it
    metadata_file = Path(self.wf_path, "ro-crate-metadata.json")
    readme_file = Path(self.wf_path, "README.md")

    # Only proceed if both files exist
    if not (metadata_file.exists() and readme_file.exists()):
        if not metadata_file.exists():
            ignored.append("`ro-crate-metadata.json` not found")
        if not readme_file.exists():
            ignored.append("`README.md` not found")
        return {"passed": passed, "fixed": fixed, "ignored": ignored}

    try:
        metadata_content = metadata_file.read_text(encoding="utf-8")
        metadata_dict = json.loads(metadata_content)
    except json.JSONDecodeError as e:
        log.error("Failed to decode JSON from `ro-crate-metadata.json`: %s", e)
        ignored.append("Invalid JSON in `ro-crate-metadata.json`")
        return {"passed": passed, "fixed": fixed, "ignored": ignored}
    readme_content = readme_file.read_text(encoding="utf-8")
    graph = metadata_dict.get("@graph")

    if not graph or not isinstance(graph, list) or not graph[0] or not isinstance(graph[0], dict):
        ignored.append("Invalid RO-Crate metadata structure.")
    else:
        # Check if the 'description' key is present
        if "description" not in graph[0]:
            metadata_dict.get("@graph")[0]["description"] = readme_content
            fixed.append("Fixed: add the same description from `README.md` to the RO-Crate metadata.")

    rc_description_graph = metadata_dict.get("@graph", [{}])[0].get("description")

    # Compare the two strings and add a linting error if they don't match
    if readme_content != rc_description_graph:
        metadata_dict.get("@graph")[0]["description"] = readme_content
        with metadata_file.open("w", encoding="utf-8") as f:
            json.dump(metadata_dict, f, indent=4)
        passed.append("RO-Crate description matches the `README.md`.")
        fixed.append("Mismatch fixed: RO-Crate description updated from `README.md`.")
    else:
        passed.append("RO-Crate descriptions are in sync with `README.md`.")
    return {"passed": passed, "fixed": fixed, "ignored": ignored}
