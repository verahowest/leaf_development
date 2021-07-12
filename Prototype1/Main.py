from default_setup import *
import visualization as vis
from growth import *
import export as export

STEPS = 30  # simulation steps
PLT_DIM = 7
INTERPOLATION = 1  # number of points to interpolate between margin points (min 1)
GR = 0.2  # growth rate
GD = 0.2  # directional growth
CP_TH = 4  # threshold distance for cp creation
KV = 0.4  # vasculatory auxin movement rate
LEAF_PATH = "../img/plot_data/"  # where to save plot data
BASE_NAME = "leaf_"  # base name of plots
CSV_PATH = "../data/leaf_export"


def main():
    """Adapted from Runions et al. 2017 approach for leaf development"""
    # leaf_param = [STEPS, INTERPOLATION, GR, CP_TH, KV]
    leaf = initialize_default_leaf()
    vis.plot_leaf_segments(leaf, PLT_DIM, LEAF_PATH, BASE_NAME, 0)
    for i in range(1, (STEPS + 1)):
        print(f"----------------------------")
        print(f"----------STEP {i}----------")
        print(f"----------------------------")

        # growstep driven by expansion of veins
        expand_veins(leaf, GR, GD, INTERPOLATION, CP_TH)
        # modification of morphogen distribution
        # introducing new cp's after margin growth
        introduce_new_cp(leaf, CP_TH, KV)
        print(f"--------------------->>>after cp addition.")
        s = 0
        for segment in leaf.segments:
            print(f"segments {s}: {segment.vein_segment}")
            s += 1
        for vein in leaf.all_veins:
            print(f"len of anchor_pts: {vein.anchor_pts} {vein.start_point} {vein.end_point}")

        # plot leaf
        vis.plot_leaf_segments(leaf, PLT_DIM, LEAF_PATH, BASE_NAME, i)
    export.export_leaf(leaf, CSV_PATH)


if __name__ == "__main__":
    main()
