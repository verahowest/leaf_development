from leafstructure import Point, Vein, Margin, Leaf



def initialize_default_leaf():
    # initialize single leaf instance
    base_point = Point([0, 0], 0, [], 0, 1)

    # create primordium vein
    primordium_vein = Vein([base_point], base_point, base_point)

    base_point.vein_assoc = [primordium_vein]

    has_morphogen = 0
    for i in range(6):
        if i == 0:  # skip base_point
            continue
        elif i == 5:
            is_cp = 1
            new_point = Point([0, i], is_cp, [primordium_vein], has_morphogen, 1)
            primordium_vein.end_point = new_point
        else:
            is_cp = 0
            new_point = Point([0, i], is_cp, [primordium_vein], has_morphogen, 0)
        primordium_vein.add_to_end(new_point)
        primordium_vein.start_point = base_point

    # create margin
    margin = Margin([], base_point, base_point)
    width_points = [0, 1.5, 2, 1.5, 1, 0]
    for i in [-1, 1]:
        for j in range(6):
            if j == 0:
                new_point = base_point
            elif j == 5:
                new_point = primordium_vein.end_point
            else:
                x = i * width_points[j]
                is_cp = 0
                has_morphogen = 0
                new_point = Point([x, j], is_cp, [primordium_vein], has_morphogen, 0)
            margin.add_to_end(new_point)
    margin.points[6:11] = reversed(margin.points[6:11])  # Reverse, So margin is one continuous line

    # --------------------------------------------------------
    default_leaf = Leaf(base_point, primordium_vein, margin, [primordium_vein])
    return default_leaf