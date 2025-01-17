from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from default_setup import initialize_default_leaf, initialize_default_leaf2
import visualization as vis
from leafstructure import Point, Vein, Margin, Leaf
from scipy.spatial import ConvexHull


# pts1 = [[0.0, 10.0], [0.0, 20.0], [10.0, 0.0], [20.0, 0.0]]
# pts = np.array(pts1)
#
# hull = ConvexHull(pts)
#
# plt.fill(pts[hull.vertices,0], pts[hull.vertices,1],'red',alpha=0.5)
# plt.show()

def plot_segments_hull(leaf):
    # plot margin segments
    _, _, segments_pos = leaf.define_segments_pos()
    # print(f"total len: {len(leaf.margin.points)}, nr of seg: {len(segments_x)}")
    for segment in segments_pos:
        pts = np.array(segment)
        hull = ConvexHull(pts)

        plt.fill(pts[hull.vertices, 0], pts[hull.vertices, 1],  color=np.random.rand(3, ), alpha=0.5)
    plt.show()
#
# x = np.linspace(0, 10, num=11, endpoint=True)
# y = np.cos(-x**2/9.0)
# # f = interp1d(x, y)
# # f2 = interp1d(x, y, kind='cubic')
#
# # xnew = np.linspace(0, 10, num=41, endpoint=True)
# #
# # plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
# # plt.legend(['data', 'linear', 'cubic'], loc='best')
# # plt.show()
#
#
# leaf = initialize_default_leaf()
# print(leaf.margin)
# points_x, points_y, _ = leaf.margin.get_points_pos()
#
# left_point = 3
# right_point = 4
# print(f"points_x{points_x}")
# print(f"points_y{points_y}")
# print(points_x[left_point:right_point+1], points_x[left_point], points_x[right_point])
#
# xnew = np.linspace(points_x[left_point], points_x[right_point], num=3, endpoint=True)
# print(xnew)
#
# f_leaf = interp1d(points_x[left_point:right_point+1], points_y[left_point:right_point+1], kind='linear')
# ynew = f_leaf(xnew)
# print(ynew)
# plt.plot(xnew, ynew, 'o', xnew, ynew, '-')
# plt.legend(['data', 'cubic'], loc='best')
# plt.show()
# vis.plot_leaf(leaf, 15, ".", "test", 0)

