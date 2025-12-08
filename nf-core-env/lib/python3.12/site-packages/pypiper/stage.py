"""Conceptualize a pipeline processing phase/stage."""

import copy

from .utils import translate_stage_name

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["Stage"]


class Stage(object):
    """
    Single stage/phase of a pipeline; a logical processing "unit". A stage is a
    collection of commands that is checkpointed.
    """

    def __init__(
        self,
        func,
        f_args=None,
        f_kwargs=None,
        name=None,
        checkpoint=True,
        *,
        nofail=False
    ):
        """
        A function, perhaps with arguments, defines the stage.

        :param callable func: The processing logic that defines the stage
        :param tuple f_args: Positional arguments for func
        :param dict f_kwargs: Keyword arguments for func
        :param str name: name for the phase/stage
        :param callable func: Object that defines how the stage will execute.
        :param bool nofail: Allow a failure of this stage to not fail the pipeline
            in which it's running
        """
        if isinstance(func, Stage):
            raise TypeError("Cannot create Stage from Stage")
        super(Stage, self).__init__()
        self.f = func
        self.f_args = f_args or tuple()
        self.f_kwargs = f_kwargs or dict()
        self.name = name or func.__name__
        self.checkpoint = checkpoint
        self.nofail = nofail

    @property
    def checkpoint_name(self):
        """
        Determine the checkpoint name for this Stage.

        :return str | NoneType: Checkpoint name for this stage; null if this
            Stage is designated as a non-checkpoint.
        """
        return translate_stage_name(self.name) if self.checkpoint else None

    def run(self, *args, **kwargs):
        """Alternate form for direct call; execute stage."""
        self(*args, **kwargs)

    def __call__(self, *args, **update_kwargs):
        """Execute the stage, allowing updates to args/kwargs."""
        kwargs = copy.deepcopy(self.f_kwargs)
        kwargs.update(update_kwargs)
        args = args or self.f_args
        self.f(*args, **kwargs)

    def __eq__(self, other):
        return (
            isinstance(other, Stage)
            and self.f.__name__ == other.f.__name__
            and (
                {k: v for k, v in self.__dict__.items() if k != "f"}
                == {k: v for k, v in other.__dict__.items() if k != "f"}
            )
        )

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return (
            "{klass} '{n}': f={f}, args={pos}, kwargs={kwd}, "
            "checkpoint={check}".format(
                klass=self.__class__.__name__,
                f=self.f,
                n=self.name,
                pos=self.f_args,
                kwd=self.f_kwargs,
                check=self.checkpoint,
            )
        )

    def __str__(self):
        return "{}: '{}'".format(self.__class__.__name__, self.name)
