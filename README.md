# pedestrian-segmentation
Pedestrian Pose Estimation via Shape Segmentation using LB Eigenfunctions for USC EE 575 Final Project

PYTHON CODE:

**segment_mesh.py**: Main Code... update file path at bottom to target .off mesh and add a name to save all matlab data
**folder_to_video.py**: Take output of animate_flooding.m and create .avi video

MATLAB PLOTTING CODE:

**Plot Cuts**: Plots segmentation 

**Plot Crits**: Plots all critical points (saddles,mins,maxes) as hollow circles and remaining critical points after noise removal as filled circles over the first eigenfunction on the mesh.

**Plot Boundary**: Plots on the level set triangles

**Plot Efuncs**: Plots the first four non-zero eigenvalue LB eigenfunctions of the mesh

**Read OFF**: Load .off files to faces and vertices

**Animate Flooding**: Produces a folder of different steps (once every 100 steps) of the flooding algorithm to create segmentation from selected extremum
