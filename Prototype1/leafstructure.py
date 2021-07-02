from operator import itemgetter, attrgetter
from functools import cmp_to_key
import numpy as np
import math

class Point:
    def __init__(self, pos, is_cp, vein_assoc, has_morphogen, has_vein):
        self.pos = pos
        self.is_cp = is_cp
        self.vein_assoc = vein_assoc
        self.has_morphogen = has_morphogen
        self.has_vein = has_vein

    def __repr__(self):
        return repr((self.pos))

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

class Vein(PointCollection):
    def __init__(self, points, start_point, end_point):
        super().__init__(points, start_point, end_point)

    def get_vector(self):
        vein_vec = np.array(self.end_point.pos) - np.array(self.start_point.pos)
        return vein_vec


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
        self.segments = self.define_segments()

    def add_vein(self, new_vein):
        self.all_veins.append(new_vein)
        self.segments = self.define_segments()

    def define_segments(self):
        """Defines segments as a collection of margin segments
        with their bounding veins."""
        segments = []
        cp_indicators, cp_index = self.margin.get_cp_indicators()
        cp_amount = len(self.margin.all_cp)
        # print(f", cp_amount: {cp_amount}")
        prev_i = 0
        for i in range(0, len(cp_index)+1):
            # print(f"cp_i : {i}")
            if i >= len(cp_index):
                end_segment_pts = self.margin.points[prev_i:]
                end_segment_pts.append(self.margin.points[0])
                end_segment = self.Segment(end_segment_pts, [prev_i, 0], cp_amount)
                segments.append(end_segment)
            else:
                segment_pts = self.margin.points[prev_i:cp_index[i]+1]
                segment = self.Segment(segment_pts, [prev_i, cp_index[i]+1], cp_amount)
                segments.append(segment)
                prev_i = cp_index[i]

        if len(segments) == (len(cp_index) + 1):
            self.segments = segments
            return segments
        else:
            return 1

    def define_segments_pos(self):
        """Returns x and y pos of segments as a collection of margin segments
        with their bounding veins."""
        # self.segments = self.define_segments()
        segments_x = []
        segments_y = []
        segments_pos = []
        for segment in self.segments:
            pts_x = []
            pts_y = []
            pts_pos = []
            for pt in segment.margin_pts_segment:
                pts_x.append(pt.pos[0])
                pts_y.append(pt.pos[1])
                pts_pos.append(pt.pos)
            segments_x.append(pts_x)
            segments_y.append(pts_y)
            segments_pos.append(pts_pos)

        return segments_x, segments_y, segments_pos

    class Segment:
        def __init__(self, margin_pts_segment, margin_slices, cp_amount):
            self.margin_pts_segment = margin_pts_segment
            self.margin_slices = margin_slices
            self.cp_amount = cp_amount
            self.vein_segment = self.find_surrounding_veins()
            self.all_pts_pos = self.get_all_pts_pos()

        def __repr__(self):
            return repr((self.margin_pts_segment, self.cp_amount, self.vein_segment, self.vein_segment, self.all_pts_pos))

        # helper function
        def find_intersection(self, left_vein, right_vein):
            """Define the segment points of two intersecting veins,
            defined by a list with their start and end points.
            - segment is enclosed by 2 veins -> returns [left_seg, right_seg]
            - segment is enclosed by 3+ veins -> returns []"""

            # left on right
            a = right_vein[1].pos
            b = right_vein[0].pos
            c = left_vein[0].pos
            left_on_right = math.dist(a, c) + math.dist(b, c) == math.dist(a, b);

            # right on left
            a = left_vein[1].pos
            b = left_vein[0].pos
            c = right_vein[0].pos
            right_on_left = math.dist(a, c) + math.dist(b, c) == math.dist(a, b);

            if right_on_left or left_on_right:
                if left_on_right:
                    # print("left_on_right")
                    left_seg = [left_vein[0].pos, left_vein[1].pos]
                    # check whether x of endpoint vein is > or < 0
                    if left_vein[1].pos[0] < 0:
                        right_seg = [left_vein[0].pos, right_vein[1].pos]
                    else:
                        right_seg = [right_vein[0].pos, left_vein[0].pos]
                if right_on_left:
                    # print("right_on_left")
                    right_seg = [right_vein[0].pos, right_vein[1].pos]
                    if right_vein[1].pos[0] < 0:
                        left_seg = [left_vein[0].pos, right_vein[0].pos]
                    else:
                        left_seg = [right_vein[0].pos, left_vein[1].pos]
                return [left_seg, right_seg]
            # TODO! check if this breaks everything
            else:
                # print("defining middle segment")
                left_seg = [left_vein[0].pos, left_vein[1].pos]
                right_seg = [right_vein[0].pos, right_vein[1].pos]
                middle_seg = [left_vein[0].pos, right_vein[0].pos]
                return [left_seg, middle_seg, right_seg]

        def find_surrounding_veins(self):
            """Given a margin segment of points with their two surrounding cp's,
            this function finds the veins that surround that segment."""

            # helper function
            def evaluate_cases(left_vein, right_vein):
                """Evaluates three cases for an iterative call:
                - only the primordium vein exists -> returns [left_vein]
                - segment is enclosed by 2 veins -> returns [left_seg, right_seg]
                - segment is enclosed by 3+ veins -> returns []"""
                # case: only primordium vein
                if (left_vein == right_vein) & (self.cp_amount == 1):
                    return [[left_vein[0].pos, left_vein[1].pos]]
                if (left_vein == right_vein) & (self.cp_amount > 1):
                    return [[]]
                else:
                    return self.find_intersection(left_vein, right_vein)

            # find surrounding veins of cp
            assoc_left = self.margin_pts_segment[0].vein_assoc
            assoc_right = self.margin_pts_segment[-1].vein_assoc
            # right_side = (self.margin_pts_segment[0].pos[0]) >= 0 and (self.margin_pts_segment[-1].pos[0]) >= 0

            # initialize base case of intersecting veins
            left_cp_vein = [assoc_left[-1].start_point, assoc_left[-1].end_point]
            right_cp_vein = [assoc_right[-1].start_point, assoc_right[-1].end_point]
            # print(f"left_cp_vein {left_cp_vein}")
            # print(f"right_cp_vein{right_cp_vein}")
            # print("evaluating cases....")
            total_seg = evaluate_cases(left_cp_vein, right_cp_vein)
            # print(f"len{len(assoc_left), len(assoc_right)}, total_seg: {total_seg}")
            # if ((len(assoc_left) == 1) & (len(assoc_right) == 1)) & (self.cp_amount >= 2):
            #     total_seg = []
            if total_seg:
                # print("gonna return total_seg")
                return total_seg
            # # recursion is necessery
            # else:
            #     total_seg = [left_cp_vein, right_cp_vein]
            #     # find intersection of the next associated veins
            #     for i in (range(2, max(len(assoc_left), len(assoc_right)))):
            #         i_left = min(len(assoc_left), i)
            #         i_right = min(len(assoc_right), i)
            #         # only consider parts within segment, so cut of the vein to the required part
            #         if right_side:
            #             # print(f"len{len(assoc_left), len(assoc_right)}, i: {i}, i-1: {-(i-1)}")
            #             left_vein = [assoc_left[-i_left].start_point, assoc_left[-(i-1)].start_point]
            #             right_vein = [assoc_right[-(i-1)].start_point, assoc_right[- i_right].end_point]
            #         else:
            #             # print(f"len{len(assoc_left), len(assoc_right)}, i: {i}, i-1: {-(i-1)}")
            #             left_vein = [assoc_left[-(i-1)].start_point, assoc_left[-i_left].end_point]
            #             right_vein = [assoc_right[- i_right].start_point, assoc_right[-(i-1)].start_point]
            #         temp_seg = evaluate_cases(left_vein, right_vein)
            #         # add found segment intersection to segment list
            #         if temp_seg:
            #             total_seg[-1:-1] = temp_seg
            #             return total_seg

        def get_all_pts_pos(self):
            """returns a string of all the positions for each point that define the boundaries of a segment.
            This is defined as [margin_seg_pos] + [vein_segment]"""

            pos = []
            for pt in self.margin_pts_segment:
                pos.append(pt.pos)
            for vein_seg in self.vein_segment:
                for vein_pt in vein_seg:
                    pos.append(vein_pt)

            return pos

