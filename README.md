# pedestrian-segmentation
Pedestrian Pose Estimation via Shape Segmentation using LB Eigenfunctions for USC EE 575 Final Project

HOW TO USE:
python estimate_pose.py --name <Add meaningful run name> --mesh <path/to/off/> will start process and prompt you to manually switch eigenfunctions. A .m file at ./runs/$name/initial_eigenfunctions.mat is saved with the initial eigenfunctions. To switch eigenfunctions, run eigenfunction_visualizer.m with the updated filename in line 1. Go back to the python script and flip the eigenfunctions if needed. The program will output the best matching pose as well as a ./runs/$name/final_results.mat file that allows you to run a bunch of different visualizations!
  
  
