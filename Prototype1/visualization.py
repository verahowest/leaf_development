import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections  as mc


# def print_points(points):
#     for i in range(0, len(points)):
#         print(str(i) + "   " + ', '.join("%s: %s" % item for item in vars(points[i]).items()))  # print class data

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

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

    plt.axis([-dim, dim, -dim/4, dim])
    plt.savefig(path + name + str(step) + ".png", format="PNG")
    plt.clf()

    # plt.show()

def plot_leaf_segments(leaf, dim, path, name, step):  # base_point, all_veins, margin, all_cp):
    """Plots an image of the leaf's veins, convergence points and margin."""

    # plot margin line
    x_margin, y_margin, margin_pos = leaf.margin.get_points_pos()
    x_margin = x_margin + [x_margin[0]]
    y_margin = y_margin + [y_margin[0]]
    plt.plot(x_margin, y_margin, '-ko')

    # plot margin segments
    segments_x, segments_y = leaf.define_segments_pos()
    # print(f"total len: {len(leaf.margin.points)}, nr of seg: {len(segments_x)}")
    for i in range(0, len(segments_x)):
        # print(f"part len: {len(segments_x[i])}")
        plt.plot(segments_x[i], segments_y[i], color=np.random.rand(3,))

    # plot veins
    for i in range(len(leaf.all_veins)):
        vein = leaf.all_veins[i]
        x = [vein.end_point.pos[0], vein.start_point.pos[0]]
        y = [vein.end_point.pos[1], vein.start_point.pos[1]]
        plt.plot(x, y, '-bo')

    # plot cp points
    x_cp, y_cp, pos_cp = leaf.margin.get_cp_pos()
    plt.plot(x_cp, y_cp, 'yo')

    plt.axis([-dim, dim, -dim/4, dim])
    plt.savefig(path + name + "_segments_"+ str(step) + ".png", format="PNG")
    plt.clf()

    # plt.show()
