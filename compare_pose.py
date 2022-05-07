from argparse import ArgumentParser
from pathlib import Path
from segmenter import Segmenter

RUNS_ROOT = Path("./runs/")

parser = ArgumentParser()
parser.add_argument("--name",required=True, type=str)
parser.add_argument("--mesh", required=True,type = str)
args = parser.parse_args()
try:
    results_dir = RUNS_ROOT / args.name
    results_dir.mkdir()
except FileExistsError:
    print("Please choose a different name or delete the previous run with the same name!")

segmenter = Segmenter(args.mesh, results_dir)
segmenter.flip_eigenfunctions()
segmenter.get_min_max_saddles()
segmenter.construct_MS_complex()
segmenter.topological_noise_removal()
#This step is a hyperparameter that controls the number of saddle points
step = max(segmenter.noise_removal_steps.keys()) - 3
segmenter.get_selected_extremum(step)
segmenter.get_landmark_matrix()
segmenter.color_mesh()
segmenter.save_to_matlab()


