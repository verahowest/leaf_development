from default_setup import *
import visualization as vis
from growth import *
import matplotlib.pyplot as plt

STEPS = 30 #simulation steps
INTERPOLATION = 2 #number of points to interpolate between margin points (min 1)
GR = 0.2 #growth rate
CP_TH = 7 #threshold distance for cp creation
KV = 0.1 #vasculatory auxin movement rate
LEAF_PATH = "../img/plot_data/" #where to save plot data
BASE_NAME = "leaf_" #base name of plots

def main():
    """Adapted from Runions et al. 2017 approach for leaf development"""
    leaf = initialize_default_leaf()
    vis.plot_leaf_segments(leaf, 15, LEAF_PATH, BASE_NAME, 0)
    for i in range(STEPS):
        print(f"----------STEP {i}----------")
        # growstep driven by expansion of veins
        expand_veins(leaf, GR, INTERPOLATION, CP_TH)
        # modification of morphogen distribution

        # introducing new cp's after margin growth
        introduce_new_cp(leaf, CP_TH, INTERPOLATION)
        print(f"len of margin: {len(leaf.margin.points)}")
        print(f"all cp len: {len(leaf.margin.all_cp)}")

        # new vein addition
        vein_addition(leaf, KV)

        # plot leaf
        vis.plot_leaf_segments(leaf, 15, LEAF_PATH, BASE_NAME, i+1)


if __name__ == "__main__":
    main()

