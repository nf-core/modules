import mimetypes
import re


def template_strings(self):
    """Check for template placeholders.

    The ``nf-core pipelines create`` pipeline template uses
    `Jinja <https://jinja.palletsprojects.com/en/2.11.x/>`_ behind the scenes.

    This lint test fails if any Jinja template variables such as
    ``{{ pipeline_name }}`` are found in your pipeline code.

    Finding a placeholder like this means that something was probably copied and pasted
    from the template without being properly rendered for your pipeline.

    This test ignores any double-brackets prefixed with a dollar sign, such as
    ``${{ secrets.AWS_ACCESS_KEY_ID }}`` as these placeholders are used in GitHub Actions workflows.

    .. tip:: You can choose to ignore lint test tests by editing the file called
        ``.nf-core.yml`` in the root of your pipeline and setting the test to false:

        .. code-block:: yaml

            lint:
                template_strings: False

        To disable this test only for specific files, you can specify a list of file paths to ignore.
        For example, to ignore a pdf you added to the docs:

        .. code-block:: yaml

            lint:
                template_strings:
                    - docs/my_pdf.pdf
    """
    passed = []
    failed = []
    ignored = []
    # Files that should be ignored according to the linting config
    ignore_files = self.lint_config.get("template_strings", []) if self.lint_config is not None else []

    files = self.list_files()
    # Loop through files, searching for string
    num_matches = 0
    for fn in files:
        if str(fn.relative_to(self.wf_path)) in ignore_files:
            ignored.append(f"Ignoring Jinja template strings in file `{fn}`")
            continue
        # Skip binary files
        binary_ftypes = ["image", "application/java-archive"]
        (ftype, encoding) = mimetypes.guess_type(fn)
        if encoding is not None or (ftype is not None and any([ftype.startswith(ft) for ft in binary_ftypes])):
            continue

        with open(fn, encoding="latin1") as fh:
            lnum = 0
            for line in fh:
                lnum += 1
                cc_matches = re.findall(r"[^$]{{[^:}]*}}", line)
                if len(cc_matches) > 0:
                    for cc_match in cc_matches:
                        failed.append(f"Found a Jinja template string in `{fn}` L{lnum}: {cc_match}")
                        num_matches += 1
    if num_matches == 0:
        passed.append(f"Did not find any Jinja template strings ({len(self.files)} files)")

    return {"passed": passed, "failed": failed, "ignored": ignored}
