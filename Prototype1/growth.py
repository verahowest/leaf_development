import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from scipy.interpolate import interp1d
from leafstructure import Point, Vein, Margin, Leaf


# HELPER FUNCTIONS

def normalize_vec(vec):
    """Normalizes a directional vector"""
    norm_vec = vec / np.sqrt((vec ** 2).sum())
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


def bounding_distances(list_of_pts, left_pt, right_pt):
    """Calculates the distance to the left and right vein for each point in a margin section.
    For this purpose the base point is treated as a convergence point."""
    pos = []

    for point in list_of_pts:
        pos.append(point.pos)
    # calculate distances to left
    euc_res_left = euclidean_distances(pos, [left_pt.pos])
    # calculate distances to right
    euc_res_right = euclidean_distances(pos, [right_pt.pos])
    # apply thresholding (for later)

    return euc_res_left, euc_res_right


# LEAF DEVELOPMENT FUNCTIONS

def interpolate_pts(pt_A, pt_B, leaf):
    """Given two positions a new point is created between these points by linear interpolation"""
    print(pt_A)
    print(pt_B)
    x = [pt_A.pos[0], pt_B.pos[0]]
    y = [pt_A.pos[1], pt_B.pos[1]]
    print(f"x: {x} y: {y}")

    xnew = np.linspace(x[0], x[1], num=3, endpoint=True)
    f = interp1d(x, y, kind='linear')
    ynew = f(xnew)
    print(f"xnew: {xnew} ynew: {ynew}")
    new_pt = Point([xnew[1], ynew[1]], 0, leaf.primordium_vein, 0, 0)
    leaf.margin.insert_point(new_pt)


def initialize_growth(leaf):
    cp_indicators, cp_index = leaf.margin.get_cp_indicators()
    print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
    # check if points need to be interpolated
    if len(cp_index) > 1:
        prev_i = 0
        for i in range(1, len(cp_index)):
            if (cp_index[i] - cp_index[prev_i]) <= 1:
                print(f"previ{cp_index[prev_i]} i{cp_index[i]}")
                interpolate_pts(leaf.margin.points[cp_index[prev_i]], leaf.margin.points[cp_index[i]], leaf)
                # recalculate cp_indicators
                cp_indicators, cp_index = leaf.margin.get_cp_indicators()
                print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
            prev_i = i
    return cp_indicators, cp_index


def hard_coded_cp_addition(leaf):
    """create new cp for vein insertion this is temporarily hard coded for poc,
    later done dynamically based on threshold during growth"""
    leaf.margin.points[2].is_cp = 1
    leaf.margin.points[3].is_cp = 1
    leaf.margin.points[7].is_cp = 1
    leaf.margin.points[8].is_cp = 1
    leaf.margin.check_conv_points()


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


def calculate_gr(prev_data, next_data, gr):
    # gr_pts = np.zeros(prev.size())
    # for prev and next
    # for i in [prev_data, next_data]:
    # temp_dist = i[0]
    # temp_dir = i[1]

    # TODO ! multiply direction * gr * normalized distance in correct math
    # next_growth = temp_dist * gr * temp_dir

    # calculate total
    # gr_pts = gr_pts + next_growth

    return 1


def expand_veins(leaf, gr):
    """Expand the cp on margin in the direction of their veins
    by a growth rate (gr)"""
    cp_indicators, cp_index = initialize_growth(leaf)
    prev_cp_i = 0
    gr_total = np.zeros(len(cp_indicators))
    # maybe change to 0 later or just to a minimal value instead?
    prev_vein_dir = -1 * normalize_vec(leaf.primordium_vein.get_vector())
    for i in range(1, len(cp_indicators)):
        if cp_indicators[i] == 1:

            # when next cp is found define next margin part
            if i == cp_index[-1]:  # if the last Cp is find connect basepoint as next cp
                next_cp_i = 0
                print(f"leaf.margin.points[{prev_cp_i + 1}:{next_cp_i}]")
                margin_part = leaf.margin.points[(prev_cp_i + 1):]
            else:
                next_cp_i = cp_index[0]
                print(f"leaf.margin.points[{prev_cp_i + 1}:{next_cp_i}]")
                if (prev_cp_i + 1) == next_cp_i:
                    margin_part = [leaf.margin.points[next_cp_i]]
                else:
                    print(f"leaf.margin.points[{prev_cp_i + 1}:{next_cp_i}]")
                    margin_part = leaf.margin.points[(prev_cp_i + 1):next_cp_i]
                # print(f"margin_part: {margin_part}")
            # print(f"i: {i}, next_cp_i = {cp_index[0]}")

            prev_cp = leaf.margin.points[prev_cp_i]
            next_cp = leaf.margin.points[next_cp_i]

            # only feed positions left and right of cp to the distance function
            prev_dist, next_dist = bounding_distances(margin_part, prev_cp, next_cp)

            # get growth direction of new cp
            exp_vein = next_cp.vein_assoc[-1]
            next_vein_dir = normalize_vec(exp_vein.get_vector())

            # add to array of gr values
            prev_data = [prev_dist, prev_vein_dir]
            next_data = [next_dist, next_vein_dir]
            gr_total[(prev_cp_i + 1):next_cp_i] = calculate_gr(prev_data, next_data, gr)
            print(f" gr_total[{prev_cp_i + 1}:{next_cp_i}]: {gr_total}")
            # move margin for next iteration
            cp_index = cp_index[1:]
            # print(f"cp_index: {cp_index}")
            prev_cp_i = next_cp_i
            prev_vein_dir = next_vein_dir

    # leaf.margin.grow(gr_total)
