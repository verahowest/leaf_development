import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from scipy.interpolate import interp1d
from leafstructure import Point, Vein, Margin, Leaf


# HELPER FUNCTIONS

def normalize_vec(vec):
    """Normalizes a directional vector"""
    norm_vec = vec / np.sqrt((vec ** 2).sum())
    return norm_vec

def normalize_to_range(vec, old_range, fit_range):
    """Normalizes a vector to a given range. as adapted from sklearn:
    X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    X_scaled = X_std * (max - min) + min"""
    vec_std = (vec - old_range[0]) / (old_range[1] - old_range[0])
    vec_scaled = vec_std * (fit_range[1] - fit_range[0]) + fit_range[0]
    return vec_scaled

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
    # print(f"pos: {pos}, left_pt.pos {[left_pt.pos]}")
    euc_res_left = euclidean_distances(pos, [left_pt.pos])
    # calculate distances to right
    euc_res_right = euclidean_distances(pos, [right_pt.pos])
    # apply thresholding (for later)

    return euc_res_left, euc_res_right


def interpolate_pts(pt_A, pt_B, leaf, mustAdd, interpolation):
    """Given two positions a new point is created between these points by linear interpolation"""
    x = [pt_A.pos[0], pt_B.pos[0]]
    y = [pt_A.pos[1], pt_B.pos[1]]
    num_interp = 2 + interpolation

    xnew = np.linspace(x[0], x[1], num=num_interp, endpoint=True)
    f = interp1d(x, y, kind='linear')
    ynew = f(xnew)

    # for each interpolation point
    if mustAdd:
        for i in range(1, interpolation + 1):
            new_pt = Point([xnew[i], ynew[i]], 0, leaf.primordium_vein, 0, 0)
            leaf.margin.insert_point(new_pt)
        return 0
    return [xnew[1], ynew[1]]


# LEAF DEVELOPMENT FUNCTIONS

def init_cp_indicators(leaf, interpolation):
    cp_indicators, cp_index = leaf.margin.get_cp_indicators()
    print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
    # check if points need to be interpolated
    if len(cp_index) > 1:
        prev_i = 0
        # interpolate for the first and last index (separate due to special case index slicing)
        if cp_index[0] == 1:
            interpolate_pts(leaf.margin.points[0], leaf.margin.points[1], leaf, True, interpolation)
            cp_indicators, cp_index = leaf.margin.get_cp_indicators()
            print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
        if cp_index[-1] == (len(leaf.margin.points) - 1):
            interpolate_pts(leaf.margin.points[cp_index[-1]], leaf.margin.points[0], leaf, True, interpolation)
            cp_indicators, cp_index = leaf.margin.get_cp_indicators()
            print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")

        # interpolate if there are no points between cp
        for i in range(1, len(cp_index)):
            if (cp_index[i] - cp_index[prev_i]) <= 1:
                interpolate_pts(leaf.margin.points[cp_index[prev_i]], leaf.margin.points[cp_index[i]], leaf, True, interpolation)
                # recalculate cp_indicators
                cp_indicators, cp_index = leaf.margin.get_cp_indicators()
                print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
            prev_i = i
    return cp_indicators, cp_index


def hard_coded_cp_addition(leaf):
    """create new cp for vein insertion this is temporarily hard coded for poc,
    later done dynamically based on threshold during growth"""
    leaf.margin.points[2].is_cp = 1
    # leaf.margin.points[3].is_cp = 1
    leaf.margin.points[7].is_cp = 1
    # leaf.margin.points[8].is_cp = 1
    leaf.margin.check_conv_points()
    return leaf


def calculate_gr(dist, dir, gr, cp_th):
    """multiply direction * gr * distance and
    normalize this to half the max distance."""

    # TODO ! improve fit, division by 0 not possible
    old_range = [-cp_th, cp_th]
    fit_range = [-cp_th/2, cp_th/2]

    temp_growth = (gr * normalize_to_range(dir, old_range, fit_range)) / dist
    # print(f"dist { dist} temp_growth: {temp_growth}")
    return temp_growth


def expand_veins(leaf, gr, interpolation, cp_th):
    """Expand the cp on margin in the direction of their veins by a growth rate (gr)"""
    # initialize growth and interpolate margin where necessary
    init_cp_indicators(leaf, interpolation)
    segments = leaf.define_segments()
    gr_total = np.zeros((len(leaf.margin.points), 2))
    prev_vein_dir = -1 * normalize_vec(leaf.primordium_vein.get_vector())

    # calculate growth rate for each segment
    for s in segments:
        # only feed positions left and right of cp to the distance function
        prev_cp = s.pts_segment[0]
        next_cp = s.pts_segment[-1]
        prev_dist, next_dist = bounding_distances(s.pts_segment[1:-1], prev_cp, next_cp)

        # get new growth dir
        if s.margin_slices[1] != 0:
            next_vein = next_cp.vein_assoc[-1]
            next_vein_dir = normalize_vec(next_vein.get_vector())
        else:
            next_vein_dir = -1 * normalize_vec(leaf.primordium_vein.get_vector())

        # add to array of gr values
        temp_prev = calculate_gr(prev_dist, prev_vein_dir, gr, cp_th)
        temp_next = calculate_gr(next_dist, next_vein_dir, gr, cp_th)
        if s.margin_slices[1] != 0:
            # add cp growth rates
            gr_total[(s.margin_slices[1])-1] = np.multiply(next_vein_dir, gr)
            # add non cp growth rates
            gr_total[(s.margin_slices[0]+1):(s.margin_slices[1]-1)] = temp_next + temp_prev
        else:
            gr_total[(s.margin_slices[0]+1):] = temp_next + temp_prev
            leaf.margin.grow(gr_total)
            return leaf
        prev_vein_dir = next_vein_dir
    return 1


def calculate_margin_distance(margin_part):
    """calculates the length of a margin section, from cp to cp"""
    dist_array = np.zeros(len(margin_part))

    prev_pt = margin_part[0]
    for i in range(len(margin_part)):
        pt = margin_part[i]
        dist_array[i] = euclidean_distances([pt.pos], [prev_pt.pos])
        prev_pt = pt

    dist_sum = dist_array.sum()

    return dist_array, dist_sum


def insert_cp(vein_assoc, margin_part, dist_array, dist_sum):
    """Inserts a new cp in the middle of two cps."""

    middle = dist_sum / 2
    temp_dist = 0
    pos = np.zeros(2)
    for i in range(1, len(dist_array)):
        temp_dist += dist_array[i]
        if temp_dist >= middle:
            pos = interpolate_pts(margin_part[i - 1], margin_part[i], None, False, 1)
            print(f"interpolated position: {pos}")
            break

    new_cp = Point(pos, 1, vein_assoc, 0, 0)

    return new_cp


def introduce_new_cp(leaf, cp_th, interpolation):
    """Introduces new cp where the boundary distance exceeds a certain distance threshold"""
    cp_indicators, cp_index = init_cp_indicators(leaf, interpolation)
    prev_cp_i = 0
    new_cps = []
    for i in range(1, len(cp_indicators)):
        # define next margin part
        if cp_indicators[i] == 1 or cp_index == []:
            if cp_index == []:
                margin_part = leaf.margin.points[prev_cp_i:(len(cp_indicators) - 1)]
                margin_part = margin_part + [leaf.margin.end_point]
            else:
                next_cp_i = cp_index[0]
                margin_part = leaf.margin.points[prev_cp_i:(next_cp_i + 1)]
            print(f"cp_index = {cp_index}")

            # insert new cp if it exceeds the margin
            dist_array, dist_sum = calculate_margin_distance(margin_part)
            if cp_th < dist_sum:
                new_cp = insert_cp([leaf.primordium_vein], margin_part, dist_array, dist_sum)
                new_cps.append(new_cp)
            print(f"number of cp's to add: {len(new_cps)}")

            if cp_index == []:
                for cp in new_cps:
                    print(f"new_cp: {cp}")
                    leaf.margin.insert_point(cp)
                    leaf.margin.check_conv_points()
                return leaf
            # prepare for next iteration
            else:
                cp_index = cp_index[1:]
                prev_cp_i = next_cp_i

    return 1


def find_optimal_angle(segment, kv):
    """Finds the optimal vein angle by applying cos-1(kv/km)"""

    theta = np.arccos(kv/1)

    return theta


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


def vein_addition(leaf, kv):
    """Connects unconnected Cp's with vein and adds new vein to leaf."""
    segments = leaf.define_segments()
    # theta_anchor_pt = find_optimal_angle(segment, kv)

    # print(f"segments: {segments}")

    for cp_i in range(len(leaf.margin.all_cp)):
        cp = leaf.margin.all_cp[cp_i]
        print(f"creating new vein for: {cp}")
        if cp.has_vein == 1:  # skip if already connected to a vein
            continue
        else:
            # create new vein
            # TODO! replace create_anchor_point by find_optimal_angle etc

            anchor_pt = create_anchor_point(cp, cp.vein_assoc)
            new_vein = Vein([anchor_pt, cp], anchor_pt, cp)

            # add vein to vein assoc
            anchor_pt.connect_to_new_vein(new_vein)
            anchor_pt.has_vein = 1
            cp.connect_to_new_vein(new_vein)
            cp.has_vein = 1

            # add vein to leaf
            leaf.add_vein(new_vein)
    return leaf
