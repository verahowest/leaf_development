import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from leafstructure import Point, Vein, Margin, Leaf

def normalize_vec(vec):
    """Normalizes a directional vector"""
    norm_vec = vec/np.sqrt((vec**2).sum())
    return norm_vec

def vector_projection(a, b):
    """Calculates the vector projection of a onto b."""
    proj = (np.dot(a, b) / np.dot(b, b)) * b
    return proj
    # k = ((y2 - y1) * (x3 - x1) - (x2 - x1) * (y3 - y1)) / ((y2 - y1) ^ 2 + (x2 - x1) ^ 2)
    # x4 = x3 - k * (y2 - y1)
    # y4 = y3 + k * (x2 - x1)

def hard_coded_cp_addition(leaf):
    # create new cp for vein insertion
    # this is temporarily hard coded for poc, later done dynamically based on threshold during growth

    leaf.margin.points[4].is_cp = 1
    leaf.margin.points[7].is_cp = 1

    leaf.margin.check_conv_points()


# def find_shortest_path(leaf, vein):
#     """Calculated euclidean distances between the CP and the given Vein."""
#     _, _, all_cp_pos = leaf.margin.get_cp_pos()
#     _, _, vein_pos = vein.get_points_pos()
#     euc_res = euclidean_distances(all_cp_pos, vein_pos)
#     min_index = np.argmin(euc_res, axis=1)
#     return min_index, euc_res

# def vein_addition(leaf, vein):
#     """Connects unconnected Cp's with vein and adds new vein to leaf."""
#     min_index, euc_res = find_shortest_path(leaf, vein)
#     for i in range(0, len(leaf.margin.all_cp)):
#         distance = euc_res[i, min_index[i]]
#         if distance == 0:  # skip if already connected to a vein
#             continue
#         else:
#             print(vein.points[min_index[i]].print_point())
#             start_point = vein.points[min_index[i]]
#             end_point = leaf.margin.all_cp[i]
#             new_vein = Vein([start_point, end_point], start_point, end_point)
#             # print(new_vein.print_points())
#             end_point.connect_to_new_vein(new_vein)

def create_anchor_point(cp, vein_assoc):
    """draws a perpendicular line to the given vein.
    Later this can be done with a theta adjustment."""
    # find anchor point
    vein = vein_assoc[-1]
    anchor = vector_projection(cp.pos, vein.get_vector())
    anchor_pt = Point(anchor.tolist(), 0, vein_assoc, 0, 0)

    # add anchor point to corresponding vein
    vein.insert_point(anchor_pt)

    return anchor_pt


def vein_addition(leaf, vein):
    """Connects unconnected Cp's with vein and adds new vein to leaf."""
    for cp in leaf.margin.all_cp:
        if cp.has_vein == 1:  # skip if already connected to a vein
            continue
        else:
            # create new vein
            anchor_pt = create_anchor_point(cp, cp.vein_assoc)
            end_point = cp
            new_vein = Vein([anchor_pt, end_point], anchor_pt, end_point)
            leaf.add_vein(new_vein)

            # add vein to vein assoc
            anchor_pt.connect_to_new_vein(new_vein)
            anchor_pt.has_vein = 1
            cp.connect_to_new_vein(new_vein)
            cp.has_vein = 1


# def normalize_vein(vein_start, vein_end):
#     """Normalized a 2D vector by its given vein start point and end point."""
#     vein_vec = np.array(vein_end) - np.array(vein_start)
#     norm_vein = normalize_vec(vein_vec)
#     return norm_vein



# TO DO
def expand_veins(leaf, gr):
    """Expand the cp on margin in the direction of their veins
    by a growth rate (gr)"""
    veins = []
    for cp in leaf.margin.all_cp:
        exp_vein = cp.vein_assoc[-1]
        # start = cp.vein_assoc[-1].start_point.pos
        # end = cp.vein_assoc.end_point.pos
        # vein_vec = np.array(end) - np.array(start)
        dir_vein = normalize_vec(exp_vein.get_vector())
        cp.pos = cp.pos + (gr * dir_vein)

