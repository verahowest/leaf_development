from default_setup import *
import visualization as vis
from growth import *
import matplotlib.pyplot as plt


# adapted from Runions et al. 2017 approach for leaf development
STEPS = 10
LEAF_PATH = "../img/plot_data/"
BASE_NAME = "leaf_"
def main():

    leaf = initialize_default_leaf()
    vis.plot_leaf(leaf, 15, LEAF_PATH, BASE_NAME, 0)
    for i in range(STEPS):
        # growstep driven by expansion of veins
        leaf = expand_veins(leaf, 0.2)
        # modification of morphogen distribution

        # new convergence points & possible new morphogen distribution
        if i == 0:
            leaf = hard_coded_cp_addition(leaf)
        else:
            leaf = introduce_new_cp(leaf, 3.5)
            print(f"len of margin: {len(leaf.margin.points)}")
            print(f"all cp len: {len(leaf.margin.all_cp)}")

        # new vein addition
        leaf = vein_addition(leaf)

        # plot leaf
        vis.plot_leaf(leaf, 15, LEAF_PATH, BASE_NAME, i+1)


if __name__ == "__main__":
    main()

