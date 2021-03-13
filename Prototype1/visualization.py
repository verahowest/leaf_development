import matplotlib.pyplot as plt
from matplotlib import collections  as mc


def print_points(points):
    for i in range(0, len(points)):
        print(str(i) + "   " + ', '.join("%s: %s" % item for item in vars(points[i]).items()))  # print class data


def plot_leaf(leaf):  # base_point, all_veins, margin, all_cp):
    """Plots an image of the leaf's veins, convergence points and margin."""

    # plot margin
    x_margin, y_margin, margin_pos = leaf.margin.get_points_pos()
    plt.plot(x_margin, y_margin, '-ko')

    # plot cp points
    x_cp, y_cp, pos_cp = leaf.margin.get_cp_pos()
    plt.plot(x_cp, y_cp, 'yo')

    # plot veins
    for vein in leaf.all_veins:
        print("next vein of length: " + str(len(vein.points)))
        x = vein.end_point.pos[0]
        y = vein.end_point.pos[1]
        plt.plot(x, y, 'bo')

    # plot format
    plt.axis([-7, 7, 0, 7])
    plt.show()
