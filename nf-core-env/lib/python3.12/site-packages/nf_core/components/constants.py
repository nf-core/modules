import os

# Constants for the nf-core/modules repo used throughout the module files
NF_CORE_MODULES_NAME = os.environ.get("NF_CORE_MODULES_NAME", "nf-core")
NF_CORE_MODULES_REMOTE = os.environ.get("NF_CORE_MODULES_REMOTE", "https://github.com/nf-core/modules.git")
NF_CORE_MODULES_DEFAULT_BRANCH = os.environ.get("NF_CORE_MODULES_DEFAULT_BRANCH", "master")
