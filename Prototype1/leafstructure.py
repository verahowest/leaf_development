class Point:
    def __init__(self, pos, is_cp, vein_assoc, has_morphogen):
        self.pos = pos
        self.is_cp = is_cp
        self.vein_assoc = vein_assoc
        self.has_morphogen = has_morphogen

    def connect_to_new_vein(self, new_vein):
        self.vein_assoc.append(new_vein)

    def print_point(self):
        """prints out the class attributes into a string"""
        print("%s: %s" % item for item in vars(self).items())


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
            print(str(i) + "   ")
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





