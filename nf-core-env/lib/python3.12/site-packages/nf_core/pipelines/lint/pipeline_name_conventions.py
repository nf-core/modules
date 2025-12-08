def pipeline_name_conventions(self):
    """Checks that the pipeline name adheres to nf-core conventions.

    In order to ensure consistent naming, pipeline names should contain only lower case, alphanumeric characters.
    Otherwise a warning is displayed.

    .. warning::
        DockerHub is very picky about image names and doesn't even allow hyphens (we are ``nfcore``).
        This is a large part of why we set this rule.
    """
    passed = []
    warned = []
    failed = []

    if self.pipeline_name.islower() and self.pipeline_name.isalnum():
        passed.append("Name adheres to nf-core convention")
    if not self.pipeline_name.islower():
        warned.append("Naming does not adhere to nf-core conventions: Contains uppercase letters")
    if not self.pipeline_name.isalnum():
        warned.append("Naming does not adhere to nf-core conventions: Contains non alphanumeric characters")

    return {"passed": passed, "warned": warned, "failed": failed}
