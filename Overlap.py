#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from prody import *
from numpy import *
from matplotlib.pyplot import *


a_aa = parsePDB('a_TgCDPK1.pdb')
a_ca = a_aa.select('calpha')

i_aa = parsePDB('i_TgCDPK1.pdb')
i_ca = i_aa.select('calpha')

anm_a = ANM('a_TgCDPK1')
anm_a.buildHessian(a_ca)
anm_a.calcModes()

anm_i = ANM('i_TgCDPK1')
anm_i.buildHessian(i_ca)
anm_i.calcModes()

calcRMSD(i_ca, a_ca)
aligned_i_ca, T = superpose(i_ca, a_ca)
calcRMSD(aligned_i_ca, a_ca)
defvec = calcDeformVector(a_ca, aligned_i_ca)
showOverlap(defvec.getNormed(), anm_a)
showCumulOverlap(defvec.getNormed(), anm_a, 'r')

calcRMSD(a_ca, i_ca)
aligned_a_ca, T = superpose(a_ca, i_ca)
calcRMSD(aligned_a_ca, i_ca)
defvec = calcDeformVector(i_ca, aligned_a_ca)
showOverlap(defvec.getNormed(), anm_i)
showCumulOverlap(defvec.getNormed(), anm_i, 'r')
