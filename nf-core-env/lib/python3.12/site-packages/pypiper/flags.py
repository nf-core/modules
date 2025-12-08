"""Status flags"""

# TODO: ultimately, these should migrate to pep.
RUN_FLAG = "running"
COMPLETE_FLAG = "completed"
FAIL_FLAG = "failed"
WAIT_FLAG = "waiting"
PAUSE_FLAG = "partial"
FLAGS = [RUN_FLAG, COMPLETE_FLAG, FAIL_FLAG, WAIT_FLAG, PAUSE_FLAG]

__all__ = ["COMPLETE_FLAG", "FAIL_FLAG", "FLAGS", "PAUSE_FLAG", "RUN_FLAG", "WAIT_FLAG"]
