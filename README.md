# pedestrian-segmentation
Pedestrian Pose Estimation via Shape Segmentation using LB Eigenfunctions for USC EE 575 Final Project <br>
Equal Work contributed by Theodore Lewitt (tedlewitt@gmail.com) and Oscar Bashaw <br>
Any improvements (especially porting the plotting code to Python) are welcomed!

## Installation

1. Download [https://github.com/Deep-MI/LaPy](LaPy) and place inside of the repo
2. Install required python libraries
``` pip install -r requirements.txt```

## Comparing Poses
```python estimate_pose.py --name <Add meaningful run name> --mesh <path/to/off/>``` <br>
This command will start segmentation process and prompt you to manually switch eigenfunctions to ensure correct pose comparisons. A .mat file at ./runs/$name/initial_eigenfunctions.mat is saved with the initial eigenfunctions. To switch eigenfunctions, run compare_eigenfunctions.m with the updated filename in line 1. Go back to the python script and flip the eigenfunctions if needed. The program will output the best matching pose as well as a ./runs/$name/all_data.mat file that allows you to run a bunch of different visualizations!
  
## Visualizing results
- **plot_eigenfunctions.m** plots the first 4 eigenfunctions of the mesh.
- **plot_segmented_mesh.m** plots the segmented mesh, giving each segement a unique color.
- **plot_critical_points.m** plots both the original critical points (unfilled circles) and the selected critical points (filled circles) after topological noise removal.
- **plot_landmarks.m** plots the ordered landmark points on the mesh. 
- **compare_eigenfunctions.m** plots the first 4 eigenfunctions of the base mesh and the supplied other mesh to easily visualize if eigenfunctions need to be flipped.
  
  
