import csv


def get_point_data(pt, leaf_part):
    pt_data = [pt.pos[0], leaf_part, pt.pos[1]]
    return pt_data


def make_leaf_data(leaf):
    # initialize
    attributes = ["P", "leaf_part", "P"]
    data_margin = []
    data_prim = []
    data_veins = []

    # import margin points
    for pt in leaf.margin.points:
        data_margin.append(get_point_data(pt, 0))
    # import primordium vein separately
    data_veins.append(get_point_data(leaf.primordium_vein.start_point, 1))
    data_veins.append(get_point_data(leaf.primordium_vein.end_point, 1))
    # import veins
    v = 1
    for vein in leaf.all_veins:
        data_veins.append(get_point_data(vein.start_point, v))
        data_veins.append(get_point_data(vein.end_point, v))
        v += 1

    data = {"margin": data_margin, "primordium_vein": data_prim, "veins": data_veins}
    return attributes, data


def write_to_csv(attributes, data, file):
    for key in data:
        with open(file + key + ".csv", 'w') as f:
            write = csv.writer(f)
            write.writerow(attributes)
            write.writerows(data[key])


def export_leaf(leaf, file):
    attributes, data = make_leaf_data(leaf)
    write_to_csv(attributes, data, file)
