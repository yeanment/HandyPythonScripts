#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This code is used to study paraview for postprocessing of OpenFOAM data.
"""

import paraview
paraview.compatibility.major = 5
# paraview.compatibility.minor = 0
from paraview.simple import *

# Reading OpenFOAM case