from default_setup import *
import visualization as vis
from growth import *


# adapted from Runions et al. 2017 approach for leaf development


def main():

    leaf = initialize_default_leaf()

    # growstep driven by expansion of veins

    # modification of morphogen distribution

    # new convergence points & possible new morphogen distribution
    hard_coded_cp_addition(leaf)

    # new vein addition
    vein_addition(leaf, leaf.primordium_vein)

    # expand_veins(leaf, 3)

    # plot leaf
    vis.plot_leaf(leaf)


if __name__ == "__main__":
    main()
