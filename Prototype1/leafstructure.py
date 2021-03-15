from operator import itemgetter, attrgetter
from functools import cmp_to_key
import numpy as np

class Point:
    def __init__(self, pos, is_cp, vein_assoc, has_morphogen, has_vein):
        self.pos = pos
        self.is_cp = is_cp
        self.vein_assoc = vein_assoc
        self.has_morphogen = has_morphogen
        self.has_vein = has_vein

    def __repr__(self):
        return repr((self.pos, self.is_cp, self.vein_assoc, self.has_morphogen, self.has_vein))

    def comparator(ptA, ptB):
        """Sort by increasing y pos, and if equal increasing x pos"""
        if ptA.pos[1] < ptB.pos[1]:
            return -1
        elif ptA.pos[1] > ptB.pos[1]:
            return 1
        else:
            if ptA.pos[0] < ptB.pos[0]:
                return -1
            elif ptA.pos[0] > ptB.pos[0]:
                return 1
        return 0

    def connect_to_new_vein(self, new_vein):
        self.vein_assoc.append(new_vein)

    def print_point(self):
        """prints out the class attributes into a string"""
        print(f"pos: {self.pos} is_cp: {self.is_cp} vein_assoc: {self.vein_assoc} has_morphogen: {self.has_morphogen}")

class PointCollection:
    def __init__(self, points, start_point, end_point):
        self.points = points
        self.start = start_point
        self.end = end_point

    def add_to_end(self, point_to_add):
        self.points.append(point_to_add)
        self.end_point = point_to_add

    def print_points(self):
        """Prints out the class data of each point."""
        for i in range(0, len(self.points)):
            # print(str(i) + "   ")
            print(self.points[i].print_point())

    def get_points_pos(self):
        """Returns lists with x_pos, y_pos, pos respectively."""
        x_pos = []
        y_pos = []
        pos = []
        for point in self.points:
            pos.append(point.pos)
            x_pos.append(point.pos[0])
            y_pos.append(point.pos[1])
        return x_pos, y_pos, pos

    def get_endpoints_pos(self):
        """Returns a list containing the positions for [start,end] points."""
        return [self.start_point.pos, self.end_point.pos]


class Vein(PointCollection):
    def __init__(self, points, start_point, end_point):
        super().__init__(points, start_point, end_point)

    def get_vector(self):
        vein_vec = np.array(self.end_point.pos) - np.array(self.start_point.pos)
        return vein_vec

    def insert_point(self, point_to_add):
        """Adds a new point to the vein, sorting it and resetting the endpoints."""
        self.points.append(point_to_add)
        self.points = sorted(self.points, key=cmp_to_key(Point.comparator))
        self.start_point = self.points[0]
        self.end_point = self.points[-1]


class Margin(PointCollection):
    def __init__(self, points, start_point, end_point):
        super().__init__(points, start_point, end_point)
        self.all_cp = []

    def check_conv_points(self):
        self.all_cp = []
        for i in self.points:
            if i.is_cp == 1:
                self.all_cp.append(i)

    def get_cp_pos(self):
        x_pos = []
        y_pos = []
        pos = []
        for point in self.all_cp:
            pos.append(point.pos)
            x_pos.append(point.pos[0])
            y_pos.append(point.pos[1])
        return x_pos, y_pos, pos

class Leaf:
    def __init__(self, base_point, primordium_vein, margin, all_veins):
        self.base_point = base_point
        self.primordium_vein = primordium_vein
        self.margin = margin
        self.all_veins = all_veins

    def add_vein(self, new_vein):
        self.all_veins.append(new_vein)





