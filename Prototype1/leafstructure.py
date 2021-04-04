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
<<<<<<< HEAD
        return repr((self.pos, self.is_cp, self.vein_assoc, self.has_morphogen, self.has_vein))
=======
        return repr((self.pos, self.is_cp, len(self.vein_assoc), self.has_morphogen, self.has_vein))
>>>>>>> origin/angled_vein_addition

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
        self.start_point = start_point
        self.end_point = end_point

    def __repr__(self):
        return repr((self.points, self.start_point, self.end_point))

    def add_to_start(self, point_to_add):
        self.points = self.points + [point_to_add]
        self.start_point = point_to_add

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

    def insert_point(self, point_to_add):
        """Adds a new point sorting it in ascending y."""
        self.points.append(point_to_add)
        self.points = sorted(self.points, key=cmp_to_key(Point.comparator))
        self.start_point = self.points[0]
        self.end_point = self.points[-1]

    # TODO
    # def grow(self, translation):
    #     """Adds a given translation defined by [x,y] to a point position"""
    #     for i in range(self.points):
    #         self.points[i].pos += translation[i]


class Vein(PointCollection):
    def __init__(self, points, start_point, end_point):
        super().__init__(points, start_point, end_point)

    def get_vector(self):
        vein_vec = np.array(self.end_point.pos) - np.array(self.start_point.pos)
        return vein_vec

<<<<<<< HEAD
    def insert_point(self, point_to_add):
        """Adds a new point to the vein, sorting it and resetting the endpoints."""
        self.points.append(point_to_add)
        self.points = sorted(self.points, key=cmp_to_key(Point.comparator))
        self.start_point = self.points[0]
        self.end_point = self.points[-1]


=======
>>>>>>> origin/angled_vein_addition
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

    def insert_point(self, point_to_add):
        """Adds a new point to the margin, sorting it and resetting the endpoints."""
        m_neg = []
        m_pos = []
        m_0 = []
        if point_to_add.pos[0] < 0:
            m_neg.append(point_to_add)
        if point_to_add.pos[0] > 0:
            m_pos.append(point_to_add)
        if point_to_add.pos[0] == 0:
            m_0.append(point_to_add)

        # split points into left and right margin based on x values
        for pt in self.points:
            if pt.pos[0] < 0:
                m_neg.append(pt)
            if pt.pos[0] > 0:
                m_pos.append(pt)
            if pt.pos[0] == 0:
                m_0.append(pt)

        # sort sides individually
        m_neg = sorted(m_neg, key=cmp_to_key(Point.comparator))
        m_pos = sorted(m_pos, key=cmp_to_key(Point.comparator), reverse=True) #sort top to bottom
        m_0 = sorted(m_0, key=cmp_to_key(Point.comparator))
        # print(f"m_neg {m_neg}, m_pos {m_pos}, m_0 {m_0}")

        # add list parts and set to attributes
        self.points = [m_0[0]] + m_neg + m_0[1:] + m_pos
        self.start_point = self.points[0]
        self.end_point = self.points[0]
        return len(self.points)

    def get_cp_indicators(self):
        """Returns a list of 1 or 0 values depending on the existence of a cp
        in the respective margin position."""
        cp_indicators = []
        cp_index = []
        for i in range(len(self.points)):
            if self.points[i].is_cp:
                cp_indicators.append(1)
                cp_index.append(i)
            else:
                cp_indicators.append(0)

        return cp_indicators, cp_index

    def grow(self, gr_total):
        """Adds the growth array to the current margin positions"""
        for i in range(len(gr_total)):
            self.points[i].pos += gr_total[i]


class Leaf:
    def __init__(self, base_point, primordium_vein, margin, all_veins):
        self.base_point = base_point
        self.primordium_vein = primordium_vein
        self.margin = margin
        self.all_veins = all_veins

    def add_vein(self, new_vein):
        self.all_veins.append(new_vein)





