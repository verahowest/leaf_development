import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from math import tan
from scipy.interpolate import interp1d
from leafstructure import Point, Vein, Margin, Leaf


# HELPER FUNCTIONS

def normalize_vec(vec):
    """Normalizes a directional vector"""
    norm_vec = vec / np.sqrt((vec ** 2).sum())
    return norm_vec

def coord_to_vec(start_pos, end_pos):
    "takes two coordinate lists and converts them into a 2D vector"
    a = start_pos[0]
    b = start_pos[1]
    vec = np.array(start_pos) - np.array(end_pos)
    return vec

def get_magnitude(vec):
    """Get the magnitude of a vector"""
    mag = np.sqrt(vec.dot(vec))
    return mag

def normalize_to_range(vec, old_range, fit_range):
    """Normalizes a vector to a given range. as adapted from sklearn:
    X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    X_scaled = X_std * (max - min) + min"""
    vec_std = (vec - old_range[0]) / (old_range[1] - old_range[0])
    vec_scaled = vec_std * (fit_range[1] - fit_range[0]) + fit_range[0]
    return vec_scaled

# def vector_projection(a, b):
#     """Calculates the vector projection of a onto b."""
#     print(f"vector_projection of a: {a}, b: {b}")
#     proj = (np.dot(a, b) / np.dot(b, b)) * b
#     return proj

def vector_projection(start_point, end_point, point):
    """Calculates the vector projection of a point onto a line."""
    point_vec = point - start_point
    line_vec = end_point - start_point
    proj = start_point + np.dot(point_vec, line_vec)/np.dot(line_vec, line_vec) * line_vec
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
    # print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
    # check if points need to be interpolated
    if len(cp_index) > 1:
        prev_i = 0
        # interpolate for the first and last index (separate due to special case index slicing)
        if cp_index[0] == 1:
            interpolate_pts(leaf.margin.points[0], leaf.margin.points[1], leaf, True, interpolation)
            cp_indicators, cp_index = leaf.margin.get_cp_indicators()
            # print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
        if cp_index[-1] == (len(leaf.margin.points) - 1):
            interpolate_pts(leaf.margin.points[cp_index[-1]], leaf.margin.points[0], leaf, True, interpolation)
            cp_indicators, cp_index = leaf.margin.get_cp_indicators()
            # print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")

        # interpolate if there are no points between cp
        for i in range(1, len(cp_index)):
            if (cp_index[i] - cp_index[prev_i]) <= 1:
                interpolate_pts(leaf.margin.points[cp_index[prev_i]], leaf.margin.points[cp_index[i]], leaf, True, interpolation)
                # recalculate cp_indicators
                cp_indicators, cp_index = leaf.margin.get_cp_indicators()
                # print(f"cp_indicators:  {cp_indicators}, cp_index: {cp_index}")
            prev_i = i
    return cp_indicators, cp_index

def calculate_gr(dist, dir, gr, dir_growth, cp_th):
    """multiply direction * gr * distance and
    normalize this to half the max distance."""

    # TODO ! improve fit, division by 0 not possible
    old_range = [-cp_th , cp_th]
    fit_range = [-cp_th/2, cp_th/2]

    temp_growth = (gr * normalize_to_range(dir, old_range, fit_range)) / dist
    # temp_growth = (gr * dir) / dist
    # print(f"dist { dist} temp_growth: {temp_growth}")
    temp_growth = temp_growth + dir_growth
    # print(f"dir_growth { dir_growth} temp_growth: {temp_growth}")
    return temp_growth


def expand_veins(leaf, gr, gd, interpolation, cp_th):
    """Expand the cp on margin in the direction of their veins by a growth rate (gr)"""
    # initialize growth and interpolate margin where necessary
    init_cp_indicators(leaf, interpolation)
    segments = leaf.define_segments()
    gr_total = np.zeros((len(leaf.margin.points), 2))
    prim_vein_dir = normalize_vec(leaf.primordium_vein.get_vector())
    prev_vein_dir = 0 * prim_vein_dir
    dir_growth = gd * prim_vein_dir

    # calculate growth rate for each segment
    for s in segments:
        # only feed positions left and right of cp to the distance function
        prev_cp = s.margin_pts_segment[0]
        next_cp = s.margin_pts_segment[-1]
        prev_dist, next_dist = bounding_distances(s.margin_pts_segment[1:-1], prev_cp, next_cp)

        # get new growth dir
        if s.margin_slices[1] != 0:
            next_vein = next_cp.vein_assoc[-1]
            next_vein_dir = normalize_vec(next_vein.get_vector())
        else:
            next_vein_dir = 0 * normalize_vec(leaf.primordium_vein.get_vector())

        # add to array of gr values
        temp_prev = calculate_gr(prev_dist, prev_vein_dir, gr, dir_growth, cp_th)
        temp_next = calculate_gr(next_dist, next_vein_dir, gr, dir_growth, cp_th)
        if s.margin_slices[1] != 0:
            # add cp growth rates
            gr_total[(s.margin_slices[1])-1] = np.multiply(next_vein_dir, gr) + (dir_growth * 2)
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

def insert_cp_pos(segment, dist_array, dist_sum):
    """Inserts a new cp in the middle of two cps."""
    margin_part = segment.margin_pts_segment
    middle = dist_sum / 2
    temp_dist = 0
    pos = np.zeros(2)
    for i in range(1, len(dist_array)):
        temp_dist += dist_array[i]
        if temp_dist >= middle:
            pos = interpolate_pts(margin_part[i - 1], margin_part[i], None, False, 1)
            # print(f"interpolated position: {pos}")
            break
    return pos

def introduce_new_cp(leaf, cp_th, interpolation, kv):
    """Introduces new cp where the boundary distance exceeds a certain distance threshold"""
    cp_indicators, cp_index = init_cp_indicators(leaf, interpolation)
    segments = leaf.define_segments()
    print(f"segments {len(segments)}")
    for segment in segments:
        print(f"segment {segment}")
        # insert new cp if it exceeds the margin
        dist_array, dist_sum = calculate_margin_distance(segment.margin_pts_segment)
        if cp_th < dist_sum:
            print(f"new cp to add at: ")
            new_cp_pos = insert_cp_pos(segment, dist_array, dist_sum)
            print(f"pos :{new_cp_pos}")
            new_vein, new_cp = create_vein(leaf, segment, kv, new_cp_pos, [leaf.primordium_vein])
            print(f"with new vein :{new_vein}")
            leaf.margin.insert_point(new_cp)
            leaf.margin.check_conv_points()
    leaf.define_segments()
    return leaf

def find_closest_vein(cp, vein_segment, theta):
    """Given a vein segment this functions find the shortest anchor point with angle theta to create a vein"""

    def find_closest_point(cp, projection, vein, proj_min, mag_min, vein_min):
        """Compares the magnitude of the cp projection on the vein segments. Returning the min projection."""
        orthogonal_side = get_magnitude(projection - np.array(cp.pos))
        if mag_min == 0 or mag_min >= orthogonal_side:
            mag_min = orthogonal_side
            proj_min = projection
            vein_min = vein
        return proj_min, mag_min, vein_min

    if theta == 0:
        anchor_pos = vein_segment[0][0]
    else:
        proj_min = []
        mag_min = 0
        vein_min = []
        # find the shortest connection between all the vein projections
        print(f"vein_segment: {vein_segment}")
        for vein_part in vein_segment:
            print(f"vein_part: {vein_part}")
            vein = coord_to_vec(vein_part[0], vein_part[1])
            projection = vector_projection(vein_part[0], vein_part[1], cp.pos)
            print(f"projection: {projection}")
            proj_min, mag_min, vein_min = find_closest_point(cp, projection, vein, proj_min, mag_min, vein_min)
        offset_dir = normalize_vec(vein_min)
        offset = mag_min / tan(theta)
        anchor_pos = proj_min + offset_dir * offset
    anchor_pt = Point(anchor_pos, 0, cp.vein_assoc, 0, 0)
    # add anchor point to corresponding vein
    # vein.insert_point(anchor_pt)

    return anchor_pt


def create_vein(leaf, segment, kv, pos, vein_assoc):
    """Connects unconnected Cp's with vein and adds new vein to leaf."""
    # segments ok but vein addition problem cp has vein = 0
    km = 1
    theta = np.arccos(kv / km)
    # create new vein
    print(f"creating new vein for: {pos}")
    new_cp = Point(pos, 1, vein_assoc, 0, 0)
    anchor_pt = find_closest_vein(new_cp, segment.vein_segment, theta)
    new_vein = Vein([anchor_pt, new_cp], anchor_pt, new_cp)
    print(f"new_vein:, {new_vein}")
    # add vein to vein assoc
    anchor_pt.has_vein = 1
    new_cp.connect_to_new_vein(new_vein)
    new_cp.has_vein = 1
    # add vein to leaf
    leaf.add_vein(new_vein)
    return new_vein, new_cp

