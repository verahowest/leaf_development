import matplotlib.pyplot as plt
from matplotlib import collections  as mc


# def print_points(points):
#     for i in range(0, len(points)):
#         print(str(i) + "   " + ', '.join("%s: %s" % item for item in vars(points[i]).items()))  # print class data


def plot_leaf(leaf, dim, path, name, step):  # base_point, all_veins, margin, all_cp):
    """Plots an image of the leaf's veins, convergence points and margin."""

    # plot margin
    x_margin, y_margin, margin_pos = leaf.margin.get_points_pos()
    x_margin = x_margin + [x_margin[0]]
    y_margin = y_margin + [y_margin[0]]
    plt.plot(x_margin, y_margin, '-ko')

    # plot veins
    for i in range(len(leaf.all_veins)):
        vein = leaf.all_veins[i]
        x = [vein.end_point.pos[0], vein.start_point.pos[0]]
        y = [vein.end_point.pos[1], vein.start_point.pos[1]]
        plt.plot(x, y, '-bo')

    # plot cp points
    x_cp, y_cp, pos_cp = leaf.margin.get_cp_pos()
    plt.plot(x_cp, y_cp, 'yo')

    plt.axis([-dim, dim, 0, dim])
    plt.savefig(path + name + str(step) + ".png", format="PNG")
    plt.clf()

    # plt.show()

