import numpy as np
import time
from functools import reduce
import taichi as ti


@ti.data_oriented
class WAVESolver:
    method_wave_equation=1
    method_wave_mesh_free=2
    methods={
        'wave_equation':method_wave_equation,
        'wave_mesh_free':method_wave_mesh_free
    }

    def __init__(self,
                 c0,
                 fz,
                 during_t,
                 source_position
                 ):
        self.c0=c0

