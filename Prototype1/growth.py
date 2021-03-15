import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from leafstructure import Point, Vein, Margin, Leaf



def hard_coded_cp_addition(leaf):
    # create new cp for vein insertion
    # this is temporarily hard coded for poc, later done dynamically based on threshold during growth

    leaf.margin.points[4].is_cp = 1
    leaf.margin.points[7].is_cp = 1

    leaf.margin.check_conv_points()


def find_shortest_path(leaf, vein):
    _, _, all_cp_pos = leaf.margin.get_cp_pos()
    _, _, vein_pos = vein.get_points_pos()
    euc_res = euclidean_distances(all_cp_pos, vein_pos)
    min_index = np.argmin(euc_res, axis=1)
    return min_index, euc_res


def connect_new_cp(leaf, vein):
    min_index, euc_res = find_shortest_path(leaf, vein)
    for i in range(0, len(leaf.margin.all_cp)):
        distance = euc_res[i, min_index[i]]
        if distance == 0:  # skip if already connected to a vein
            continue
        else:
            start_point = vein.points[min_index[i]]
            end_point = leaf.margin.all_cp[i]
            new_vein = Vein([start_point, end_point], start_point, end_point)
            end_point.connect_to_new_vein(new_vein)


