"""
Constant variables for yacman package
"""

IK = "__internal"

FILEPATH_KEY = "file_path"
RO_KEY = "ro"
USE_LOCKS_KEY = "locks"
ORI_STATE_KEY = "ori_state"
WAIT_MAX_KEY = "wait_time"
ALIASES_KEY = "aliases"
ALIASES_KEY_RAW = "aliases_raw"
WRITE_VALIDATE_KEY = "write_validate"
SCHEMA_KEY = "schema"

ATTR_KEYS = (
    IK,
    USE_LOCKS_KEY,
    FILEPATH_KEY,
    RO_KEY,
    ORI_STATE_KEY,
    WAIT_MAX_KEY,
    ALIASES_KEY,
    ALIASES_KEY_RAW,
    WRITE_VALIDATE_KEY,
    SCHEMA_KEY,
)

LOCK_PREFIX = "lock."
DEFAULT_RO = False
DEFAULT_WAIT_TIME = 60
