# -*- coding: utf-8 -*-
"""
Created on Thu May 11 01:44:09 2023

@author: ppynr2
"""

import pymeshlab
import os

path = r'R:\DRS-KidsOPM\Paediatric_OPM_Notts_AdultData_individual\Data\BIDS\derivatives\coregistration\sub-101'
os.chdir(path)
ms = pymeshlab.MeshSet()

ms.load_new_mesh(path + '\sub-101_scalp_points.ply')

ms.compute_normal_for_point_clouds(k=8, smoothiter=100)

ms.generate_surface_reconstruction_screened_poisson(pointweight=40)

ms.save_current_mesh(path + '\sub-101_scalp_mesh_test.ply', save_face_color=True)

