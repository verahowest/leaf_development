from default_setup import *
import visualization as vis
from growth import *
import matplotlib.pyplot as plt


# adapted from Runions et al. 2017 approach for leaf development
STEPS = 2
LEAF_PATH = "../img/plot_data/"
BASE_NAME = "leaf_"
def main():

    leaf = initialize_default_leaf()
    for i in range(STEPS):
        # growstep driven by expansion of veins
        expand_veins(leaf, 0.5)
        # modification of morphogen distribution

        # new convergence points & possible new morphogen distribution
        hard_coded_cp_addition(leaf)

        # new vein addition
        vein_addition(leaf, leaf.primordium_vein)

        # plot leaf
        vis.plot_leaf(leaf, 15, LEAF_PATH, BASE_NAME, i)


if __name__ == "__main__":
    main()
