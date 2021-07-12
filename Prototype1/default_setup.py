from leafstructure import Point, Vein, Margin, Leaf
from growth import interpolate_pts

def initialize_default_leaf():
    # initialize single leaf instance
    base_point = Point([0, 0], 0, [], 0, 1)
    final_point = Point([0, 5], 1, [], 0, 1)

    # create primordium vein
    primordium_vein = Vein([base_point, final_point], base_point, final_point)
    final_point.vein_assoc = [primordium_vein]
    primordium_vein.insert_point(final_point)

    base_point.vein_assoc = [primordium_vein]

    # create margin
    margin = Margin([base_point, final_point], base_point, base_point)
    width_points = [0, 1.5, 2, 1.75, 1, 0]
    for i in [-1, 1]:
        for j in range(1, 5):
            x = i * width_points[j]
            is_cp = 0
            has_morphogen = 0
            new_point = Point([x, j], is_cp, [primordium_vein], has_morphogen, 0)
            margin.insert_point(new_point)
    margin.check_conv_points()

    default_leaf = Leaf(base_point, primordium_vein, margin, [primordium_vein])

    # interpolate to create more points on margin
    ptA = default_leaf.margin.start_point
    iter_pts = default_leaf.margin.points[1::] + [default_leaf.margin.end_point]
    for ptB in iter_pts:
        interpolate_pts(ptA, ptB, default_leaf, True, interpolation=1)
        ptA = ptB

    return default_leaf