"""Custom pypiper exceptions"""

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = [
    "PipelineError",
    "PipelineHalt",
    "IllegalPipelineDefinitionError",
    "IllegalPipelineExecutionError",
    "MissingCheckpointError",
    "UnknownPipelineStageError",
    "UnsupportedFiletypeException",
    "SubprocessError",
]


class PipelineError(Exception):
    """General pipeline error."""

    pass


class SubprocessError(Exception):
    pass


class IllegalPipelineDefinitionError(PipelineError):
    pass


class IllegalPipelineExecutionError(PipelineError):
    """Represent cases of illogical start/stop run() declarations."""

    pass


class MissingCheckpointError(Exception):
    """Represent case of expected but absent checkpoint file."""

    def __init__(self, checkpoint, filepath):
        msg = "{}: '{}'".format(checkpoint, filepath)
        super(MissingCheckpointError, self).__init__(msg)


class UnknownPipelineStageError(Exception):
    """
    Triggered by use of unknown/undefined name for a pipeline stage.

    :param str stage_name: Name of the stage triggering the exception.
    :param pypiper.Pipeline pipeline: Pipeline for which the stage is unknown/undefined.
    """

    def __init__(self, stage_name, pipeline=None):
        message = stage_name
        if pipeline is not None:
            try:
                stages = pipeline.stages()
            except AttributeError:
                # Just don't contextualize the error with known stages.
                pass
            else:
                message = "{}; defined stages: {}".format(
                    message, ", ".join(map(str, stages))
                )
        super(UnknownPipelineStageError, self).__init__(message)


class PipelineHalt(Exception):
    """
    Execution-stopping exception for halting a pipeline.

    This is useful for stopping execution of a truly script-like pipeline.
    That is, a pipeline that doesn't bundle/define stages or wrap run() calls
    in functions. In this case, we want to be able to stop the Python process
    as it chugs through a pipeline script, and we can do that by having a
    PipelineManager's halt method raise this exception.

    """

    def __init__(self, checkpoint=None, finished=None):
        if checkpoint is None:
            super(PipelineHalt, self).__init__()
        else:
            if isinstance(checkpoint, str):
                last_stage_done = checkpoint
            else:
                last_stage_done = getattr(checkpoint, "name", None) or getattr(
                    checkpoint, "__name__", None
                )
            if not last_stage_done:
                super(PipelineHalt, self).__init__()
            else:
                if finished is None:
                    msg = last_stage_done
                elif finished:
                    msg = "Finished '{}'".format(last_stage_done)
                else:
                    msg = "Stopped at '{}'".format(last_stage_done)
                super(PipelineHalt, self).__init__(msg)


class UnsupportedFiletypeException(Exception):
    """Restrict filetype domain."""

    # Use superclass ctor to allow file name/path or extension to pass
    # through as the message for why this error is occurring.
    pass
