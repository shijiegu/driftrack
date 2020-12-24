"""affinewarp - Time warping under affine warping functions"""

__version__ = '0.1.0'
__author__ = 'Alex Williams <ahwillia@stanford.edu>, Shijie Gu <shijiegu@berkeley.edu>'

import sys
import warnings
if sys.version_info[0] != 3:
    warnings.warn(
        "This package is supported for Python 3.+, "
        "you are currently using {}.{}.{}".format(
            sys.version_info[0],
            sys.version_info[1],
            sys.version_info[2]
        )
    )

from .piecewisewarp import PiecewiseWarping
from .rank1_piecewisewarp import Rank1_PiecewiseWarping
from .shiftwarp import ShiftWarping
from .shiftwarp_PCA import ShiftWarping_PCA
from .spikedata import SpikeData
