import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from leafstructure import Point, Vein, Margin, Leaf


# HELPER FUNCTIONS

def normalize_vec(vec):
    """Normalizes a directional vector"""
    norm_vec = vec/np.sqrt((vec**2).sum())
    return norm_vec

def vector_projection(a, b):
    """Calculates the vector projection of a onto b."""
    proj = (np.dot(a, b) / np.dot(b, b)) * b
    return proj


def distance_matrix(margin, cp):
    """Calculates euclidean distance between each point and a given convergence point"""
    _, _, all_cp_pos = margin.get_points_pos()
    euc_res = euclidean_distances(all_cp_pos, [cp.pos])
    return euc_res

def hard_coded_cp_addition(leaf):
    """create new cp for vein insertion this is temporarily hard coded for poc,
    later done dynamically based on threshold during growth"""
    leaf.margin.points[4].is_cp = 1
    leaf.margin.points[10].is_cp = 1
    leaf.margin.check_conv_points()


# LEAF DEVELOPMENT FUNCTIONS

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
            new_vein = Vein([anchor_pt, cp], anchor_pt, cp)

            # add vein to vein assoc
            anchor_pt.connect_to_new_vein(new_vein)
            anchor_pt.has_vein = 1
            cp.connect_to_new_vein(new_vein)
            cp.has_vein = 1

            # add vein to leaf
            leaf.add_vein(new_vein)

# def calculate_gr(distance_matrix, gr, th):
#


def expand_veins(leaf, gr):
    """Expand the cp on margin in the direction of their veins
    by a growth rate (gr)"""
    veins = []
    for i in range(1, len(leaf.margin.all_cp)):
        cp = leaf.margin.all_cp[i]
        # only feed positions left and right of cp to the distance function
        # multiply these by growth rate in direction
        # add to existing gr values
        # at the end of the loop expand all in 1 go

        # print(distance_matrix(leaf.margin, cp))

        exp_vein = cp.vein_assoc[-1]
        dir_vein = normalize_vec(exp_vein.get_vector())
        cp.pos = cp.pos + (gr * dir_vein)

