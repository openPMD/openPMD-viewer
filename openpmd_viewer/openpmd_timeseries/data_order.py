"""
This file is part of the openPMD-viewer.

It defines an ordering Enum to allow for flexible ordering of components

Copyright 2022, openPMD-viewer contributors
Authors: Axel Huebl, Ryan Sandberg
License: 3-Clause-BSD-LBNL
"""

import enum
 
# creating enumerations using class
class RZorder(enum.Enum):
    """The index names of RZ axes in C order"""
    mrz = 1  # z is the fastest varying index in memory
    mzr = 2  # r is the fastest varying index in memory


order_error_msg = 'Data order is unupported. '
order_error_msg += 'Allowed orderings: '
order_error_msg += ' '.join([val.name for val in RZorder])
