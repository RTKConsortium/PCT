from itk import pctConfig
from itk.rtkargumentparser import RTKArgumentParser

"""PCT application argument parser.

``PCTArgumentParser`` reuses :class:`RTKArgumentParser` and only overrides
the version string so that ``--version`` reports PCT's version.
"""

__all__ = ["PCTArgumentParser"]


class PCTArgumentParser(RTKArgumentParser):
    """Argument parser for PCT Python applications.

    Reuses ``RTKArgumentParser`` and only overrides the version string so that
    ``--version`` reports PCT's version.
    """

    def __init__(self, description=None, version=None, **kwargs):
        super().__init__(
            description=description,
            version=version or pctConfig.PCT_GLOBAL_VERSION_STRING,
            **kwargs,
        )
