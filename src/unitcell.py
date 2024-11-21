import numpy as np
import os
import copy
import json
import molmath as mm
import sm_new
import molecule as mol


class Unitcell:
    def __init__(self,**kwargs):
        # self.a = 1
        # self.b = 1
        # self.c = 1
        # self.alpha = np.pi/2
        # self.beta = np.pi/2
        # self.gamma = np.pi/2
        #vertices should be generated from cell params.
        self.tv1 = np.array([])
        self.tv2 = np.array([])
        self.tv2 = np.array([])

