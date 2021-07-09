from operator import itemgetter, attrgetter
from functools import cmp_to_key
import numpy as np
import math
import itertools
from collections import OrderedDict


class Point:
    def __init__(self, pos, is_cp, vein_assoc, has_morphogen, has_vein, neighbours=[None, None]):
        self.pos = pos
        self.is_cp = is_cp
        self.vein_assoc = vein_assoc
        self.has_morphogen = has_morphogen
        self.has_vein = has_vein
        self.neighbours = neighbours

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
    def __init__(self, points, start_point, end_point, anchor_pts=[]):
        super().__init__(points, start_point, end_point)
        self.anchor_pts = anchor_pts

    def get_vector(self):
        vein_vec = np.array(self.end_point.pos) - np.array(self.start_point.pos)
        return vein_vec

    def define_parts(self):
        parts = []
        if len(self.anchor_pts) == 0:
            parts.append([self.start_point, self.end_point])
        else:
            start = self.start_point
            for anchor in self.anchor_pts:
                parts.append([start, anchor])
                start = anchor
            parts.append([start, self.end_point])
        return parts


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

    def add_neighbours(self, new_pts, pt_left=None, pt_right=None, is_interpolation=False):
        # add multiple interpolated points or single pt
        if is_interpolation:
            prev_pt = pt_left
            margin_index = self.points.index(prev_pt)
            for i in range(len(new_pts)):
                if i < len(new_pts) - 1:
                    new_pts[i].neighbours = [prev_pt, new_pts[i + 1]]
                    prev_pt = new_pts[i]
                else:
                    new_pts[i].neighbours = [prev_pt, pt_right]
                # update neighbours in margin
                self.points[margin_index].neighbours[1] = [new_pts[i]]
                # print(f"margin_index: {margin_index} len(new_pts): {len(new_pts)}, len(points): {len(self.points)}")
                # print(f"self.points[margin_index+1]: {self.points[margin_index+1]}")
                if not margin_index + 1 == len(self.points):
                    self.points[margin_index + 1].neighbours[0] = [new_pts[i]]
                # insert new point
                self.points.insert(margin_index + 1, new_pts[i])
        if not is_interpolation:
            for new_pt in new_pts:
                margin_index = self.points.index(new_pt.neighbours[0])
                self.points[margin_index].neighbours[1] = [new_pt]
                self.points[margin_index + 1].neighbours[0] = [new_pt]
                self.points.insert(margin_index + 1, new_pt)

    def insert_point(self, point_to_add):
        """ DEPRECATED : DO NOT USE except for initialization!
        Adds a new point to the margin, sorting it and resetting the endpoints."""
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
        m_pos = sorted(m_pos, key=cmp_to_key(Point.comparator), reverse=True)  # sort top to bottom
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


class Leaf:
    def __init__(self, base_point, primordium_vein, margin, all_veins):
        self.base_point = base_point
        self.primordium_vein = primordium_vein
        self.margin = margin
        self.all_veins = all_veins
        self.segments = self.define_segments()

    def add_vein(self, new_vein, do_define=True):
        self.all_veins.append(new_vein)
        if do_define:
            self.segments = self.define_segments()

    def remove_vein(self, vein_to_remove):
        self.all_veins.remove(vein_to_remove)

    def grow(self, gr_total, anchors=None):
        """First adds the growth array to the current margin positions,
        Then grows the anchor points in their respective direction"""
        for i in range(len(gr_total)):
            self.margin.points[i].pos += gr_total[i]
        # if anchors:
        #     for anchor_pts, growth in anchors:
        #         for pt in anchor_pts:
        #             pt.pos += growth

    def print_all_anchor_pts(self):
        for vein in self.all_veins:
            for anchor_pt in vein.anchor_pts:
                print(anchor_pt)

    def find_anchor_pt(self, new_pos, new_vein_assoc):
        """checks if an anchor point already exists and otherwise creates a new one.
        If it creates a new one it also returns the True boolean."""
        for vein in self.all_veins:
            for anchor_pt in vein.anchor_pts:
                if np.allclose(anchor_pt.pos, new_pos):
                    return anchor_pt, new_vein_assoc.append(vein), False
        # no anchor point was found at that position
        anchor_pt = Point(new_pos, 0, new_vein_assoc, 0, 0)
        return anchor_pt, new_vein_assoc, True

    def find_vein_for_part(self, vein_part):
        """When given a part between two anchor points  it finds which vein it is"""
        pt_a = vein_part[0]
        pt_b = vein_part[1]
        for vein in self.all_veins:
            search_points = vein.anchor_pts.copy()
            search_points.append(vein.start_point)
            search_points.append(vein.end_point)
            print(f"search_points: {search_points}")
            if pt_a in search_points and pt_b in search_points:
                print("found vein_part in vein")
                return vein

    def define_segments(self):
        """Defines segments as a collection of margin segments
        with their bounding veins."""
        segments = []
        cp_indicators, cp_index = self.margin.get_cp_indicators()
        cp_amount = len(self.margin.all_cp)
        # print(f", cp_amount: {cp_amount}")
        prev_i = 0
        for i in range(0, len(cp_index) + 1):
            # print(f"cp_i : {i}")
            if i >= len(cp_index):
                end_segment_pts = self.margin.points[prev_i:]
                end_segment_pts.append(self.margin.points[0])
                end_segment = self.Segment(end_segment_pts, [prev_i, 0], cp_amount)
                segments.append(end_segment)
            else:
                segment_pts = self.margin.points[prev_i:cp_index[i] + 1]
                segment = self.Segment(segment_pts, [prev_i, cp_index[i] + 1], cp_amount)
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
            # self.vein_segment = self.find_surrounding_veins()
            self.vein_segment = self.find_surrounding_veins_2()
            self.all_pts_pos = self.get_all_pts_pos()

        def __repr__(self):
            return repr(
                (self.margin_pts_segment, self.cp_amount, self.vein_segment, self.vein_segment, self.all_pts_pos))

        # @staticmethod
        # def find_intersection(left_vein, right_vein):
        #     """Define the segment points of two intersecting veins,
        #     defined by a list with their start and end points.
        #     - segment is enclosed by 2 veins -> returns [left_seg, right_seg]
        #     - segment is enclosed by 3+ veins -> returns []"""
        #
        #     # left on right
        #     a = right_vein[1].pos
        #     b = right_vein[0].pos
        #     c = left_vein[0].pos
        #     left_on_right = math.dist(a, c) + math.dist(b, c) == math.dist(a, b);
        #
        #     # right on left
        #     a = left_vein[1].pos
        #     b = left_vein[0].pos
        #     c = right_vein[0].pos
        #     right_on_left = math.dist(a, c) + math.dist(b, c) == math.dist(a, b);
        #
        #     if right_on_left or left_on_right:
        #         if left_on_right:
        #             # print("left_on_right")
        #             left_seg = [left_vein[0].pos, left_vein[1].pos]
        #             # check whether x of endpoint vein is > or < 0
        #             if left_vein[1].pos[0] < 0:
        #                 right_seg = [left_vein[0].pos, right_vein[1].pos]
        #             else:
        #                 right_seg = [right_vein[0].pos, left_vein[0].pos]
        #         if right_on_left:
        #             # print("right_on_left")
        #             right_seg = [right_vein[0].pos, right_vein[1].pos]
        #             if right_vein[1].pos[0] < 0:
        #                 left_seg = [left_vein[0].pos, right_vein[0].pos]
        #             else:
        #                 left_seg = [right_vein[0].pos, left_vein[1].pos]
        #         return [left_seg, right_seg]
        #     # TODO! check if this breaks everything
        #     else:
        #         # print("defining middle segment")
        #         left_seg = [left_vein[0].pos, left_vein[1].pos]
        #         right_seg = [right_vein[0].pos, right_vein[1].pos]
        #         middle_seg = [left_vein[0].pos, right_vein[0].pos]
        #         return [left_seg, middle_seg, right_seg]

        @staticmethod
        def check_if_intersecting(l_vein_in, r_vein_in, ass_index):
            if ass_index <= len(l_vein_in):
                l_vein = l_vein_in[-ass_index]
            else:
                l_vein = l_vein_in[0]
            if ass_index <= len(r_vein_in):
                r_vein = r_vein_in[-ass_index]
            else:
                r_vein = r_vein_in[0]

            is_inters = False
            inters = None
            # check if start left vein is somewhere on the right vein
            if (l_vein.start_point in r_vein.anchor_pts) or (l_vein.start_point in r_vein.points):
                is_inters = True
                inters = l_vein.start_point
            # check if the left anchor points are somewhere on the right vein
            else:
                if len(l_vein.anchor_pts) > 0:
                    for anchor_pt in l_vein.anchor_pts:
                        is_inters = (anchor_pt in r_vein.anchor_pts) or (anchor_pt in r_vein.points)
                        if is_inters:
                            inters = anchor_pt
                            break
            return is_inters, inters

        @staticmethod
        def add_to_vein_seg(l_vein, r_vein, vein_seg, inters_pt):
            """if these parts are intersecting, add all the vein segments up to the intersection"""
            l_connected = False
            i = 1
            while not l_connected:
                # exception for the primordium vein
                if np.allclose(np.array(l_vein[0][0].pos), np.array([0, 0])) and r_vein[-1][1].pos[0] < 0:
                    # print("pos is 0,0")
                    l_part = l_vein[0]
                else:
                    l_part = l_vein[-i]
                # print(f"--?l_part {l_part} inters_pt: {inters_pt}")
                if l_part not in vein_seg:
                    vein_seg.append(l_part)
                # if it contains the intersection stop
                if inters_pt in l_part:
                    # print("l_connected")
                    l_connected = True
                i = i + 1
            r_connected = False
            i = 1
            while not r_connected:
                # exception for the primordium vein
                if np.allclose(np.array(r_vein[0][0].pos), np.array([0, 0])) and l_vein[-1][1].pos[0] > 0:
                    # print("pos is 0,0")
                    r_part = r_vein[0]
                else:
                    r_part = r_vein[-i]
                # print(f"--?r_part {r_part} inters_pt: {inters_pt}")
                if r_part not in vein_seg:
                    vein_seg.append(r_part)
                # if it starts in the intersection stop
                if inters_pt in r_part:
                    # print("r_connected")
                    r_connected = True
                i = i + 1
            if l_connected and r_connected:
                return vein_seg

        @staticmethod
        def cut_vein_assoc(left_part, left_cut, right_part, right_cut):
            # TODO! primordium exception
            i = 0
            for part in left_part:
                i += 1
                if left_cut == part[0]:
                    break
            left_part = left_part[i-1::]
            j = 0
            for part in right_part:
                j += 1
                if right_cut==part[1]:
                    break
            right_part = right_part[:j]
            print(f"i: {i} j: {j}")

            return left_part, right_part

        def find_surrounding_veins_2(self):
            """Given a margin segment of points with their two surrounding cp's,
            this function finds the veins that surround that segment."""

            # find surrounding veins of cp
            assoc_left = self.margin_pts_segment[0].vein_assoc
            assoc_right = self.margin_pts_segment[-1].vein_assoc
            print(f"len(assoc_left): {len(assoc_left)}, len(assoc_right): {len(assoc_right)}")
            vein_segment = []
            assoc = 1
            is_intersecting, intersection = self.check_if_intersecting(assoc_left, assoc_right, assoc)
            print(f"is_intersecting: {is_intersecting} intersection: {intersection}")
            while not is_intersecting:
                temp_left = assoc_left[-assoc].define_parts()
                left_cut = temp_left[0][0]
                temp_right = assoc_right[-assoc].define_parts()
                right_cut = temp_right[0][0]
                # print(f"temp_left {temp_left}, left_cut: {left_cut}")
                # print(f"temp_right {temp_right}, right_cut: {right_cut}")

                # this is new
                if assoc < len(assoc_left):
                    vein_segment.append(temp_left)
                if assoc < len(assoc_right):
                    vein_segment.append(temp_right)
                # vein_segment.append(temp_left)
                # vein_segment.append(temp_right)
                print(f"vein_segment while not intersecting: {vein_segment}")
                assoc += 1
                is_intersecting, intersection = self.check_if_intersecting(assoc_left, assoc_right, assoc)
                print(f"is_intersecting: {is_intersecting} intersection: {intersection}")

            if is_intersecting:
                if assoc > 1:
                    if assoc <= len(assoc_left):
                        left_part = assoc_left[-assoc]
                    else:
                        left_part = assoc_left[0]
                    if assoc <= len(assoc_right):
                        right_part = assoc_right[-assoc]
                    else:
                        right_part = assoc_right[0]

                    left_part = left_part.define_parts()
                    right_part = right_part.define_parts()
                    print(f"<<<<<<intersection for assoc > 1")
                    print(f"left_part: {left_part}")
                    print(f"right_part: {right_part}")

                    if left_part == right_part:
                        # TODO! update alll ditt
                        i = 0
                        cut = []
                        # print(f"left_cut: {left_cut}, right_cut: {right_cut}")
                        for part in left_part:
                            i += 1
                            # print(f"left_cut: {left_cut}, part[1]: {np.array(part[1])}")
                            if left_cut == part[1] or right_cut == part[1]:
                                cut.append(i)
                        middle_part = left_part[cut[0]:cut[1]]
                        # print(f"middle_part: {middle_part}")
                        vein_segment.append(middle_part)

                    else:
                        print(f"left_cut: {left_cut}, right_cut: {right_cut}")
                        left_part, right_part = self.cut_vein_assoc(left_part, left_cut, right_part, right_cut)
                        print(f"-> the cut left_part: {left_part}")
                        print(f"-> the cut right_part: {right_part}")
                        # intersection = left_cut

                        vein_segment = self.add_to_vein_seg(left_part, right_part, vein_segment, intersection)

                else:
                    left_part = assoc_left[-assoc].define_parts()
                    right_part = assoc_right[-assoc].define_parts()
                    vein_segment = self.add_to_vein_seg(left_part, right_part, vein_segment, intersection)
                # print(f"left_part: {left_part}")
                # print(f"right_part: {right_part}")

                print(f"????FINAL VEIN SEGMENT: {vein_segment}")
                return vein_segment
            # case left vein and right vein don't connect -> recursion (maybe this isn't necessary yet)


        # def find_surrounding_veins(self):
        #     """Given a margin segment of points with their two surrounding cp's,
        #     this function finds the veins that surround that segment."""
        #
        #     # helper function
        #     def evaluate_cases(left_vein, right_vein):
        #         # case: only primordium vein
        #         if (left_vein == right_vein) & (self.cp_amount == 1):
        #             return [[left_vein[0].pos, left_vein[1].pos]]
        #         if (left_vein == right_vein) & (self.cp_amount > 1):
        #             return [[]]
        #         else:
        #             return self.find_intersection(left_vein, right_vein)
        #
        #     # find surrounding veins of cp
        #     assoc_left = self.margin_pts_segment[0].vein_assoc
        #     assoc_right = self.margin_pts_segment[-1].vein_assoc
        #
        #     # initialize base case of intersecting veins
        #     left_cp_vein = [assoc_left[-1].start_point, assoc_left[-1].end_point]
        #     right_cp_vein = [assoc_right[-1].start_point, assoc_right[-1].end_point]
        #     total_seg = evaluate_cases(left_cp_vein, right_cp_vein)
        #     if total_seg:
        #         return total_seg

        def get_all_pts_pos(self):
            """returns a string of all the positions for each point that define the boundaries of a segment.
            This is defined as [margin_seg_pos] + [vein_segment]"""

            pos = []
            for pt in self.margin_pts_segment:
                pos.append(np.array(pt.pos).astype(float))
            for vein_seg in self.vein_segment:
                for vein_pt in vein_seg:
                    if isinstance(vein_pt, list):
                        for pt in vein_pt:
                            pos.append(np.array(pt.pos).astype(float))
                    else:
                        pos.append(np.array(vein_pt.pos).astype(float))


            return pos

