from default_setup import *
import visualization as vis
from growth import *
import matplotlib.pyplot as plt

STEPS = 15 #simulation steps
INTERPOLATION = 3 #number of points to interpolate between margin points (min 1)
GR = 0.3
CP_TH = 6.5
KV = 1 #vasculatory auxin movement rate
LEAF_PATH = "../img/plot_data/" #where to save plot data
BASE_NAME = "leaf_" #base name of plots

def main():
    """Adapted from Runions et al. 2017 approach for leaf development"""
    leaf = initialize_default_leaf()
    vis.plot_leaf(leaf, 15, LEAF_PATH, BASE_NAME, 0)
    for i in range(STEPS):
        # growstep driven by expansion of veins
        expand_veins(leaf, GR, INTERPOLATION, CP_TH)
        # modification of morphogen distribution

        # new convergence points & possible new morphogen distribution
        if i < 0:
            hard_coded_cp_addition(leaf)
        else:
            introduce_new_cp(leaf, CP_TH, INTERPOLATION)
            print(f"len of margin: {len(leaf.margin.points)}")
            print(f"all cp len: {len(leaf.margin.all_cp)}")

        # new vein addition
        vein_addition(leaf, KV)

        # plot leaf
        vis.plot_leaf_segments(leaf, 15, LEAF_PATH, BASE_NAME, i+1)


if __name__ == "__main__":
    main()

