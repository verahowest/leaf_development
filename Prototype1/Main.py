from default_setup import *
import visualization as vis
from growth import *



def main():

    leaf = initialize_default_leaf()

    # growstep
    hard_coded_cp_addition(leaf)
    connect_new_cp(leaf, leaf.primordium_vein)




    # print_points(margin)
    # vis.print_points(leaf.margin.all_cp)

    # plot leaf
    vis.plot_leaf(leaf)


if __name__ == "__main__":
    main()
