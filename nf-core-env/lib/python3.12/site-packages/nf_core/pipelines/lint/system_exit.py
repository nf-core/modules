import logging
from pathlib import Path

log = logging.getLogger(__name__)


def system_exit(self):
    """Check for System.exit calls in groovy/nextflow code

    Calls to System.exit(1) should be replaced by throwing errors

    This lint test looks for all calls to `System.exit`
    in any file with the `.nf` or `.groovy` extension
    """
    passed = []
    warned = []

    root_dir = Path(self.wf_path)

    # Get all groovy and nf files
    groovy_files = [f for f in root_dir.rglob("*.groovy")]
    nf_files = [f for f in root_dir.rglob("*.nf")]
    to_check = nf_files + groovy_files

    for file in to_check:
        try:
            with file.open() as fh:
                for i, line in enumerate(fh.readlines(), start=1):
                    if "System.exit" in line and "System.exit(0)" not in line:
                        warned.append(f"`System.exit` in {file.name}: _{line.strip()}_  [line {i}]")
        except FileNotFoundError:
            log.debug(f"Could not open file {file.name} in system_exit lint test")

    if len(warned) == 0:
        passed.append("No `System.exit` calls found")

    return {"passed": passed, "warned": warned}
