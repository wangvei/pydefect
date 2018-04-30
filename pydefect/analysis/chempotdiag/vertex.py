from __future__ import print_function
from collections import OrderedDict, Iterable, UserList
import json
import string

import numpy as np


class Vertex:
    """
        Object for vertex in chemical potential diagram.
    """

    def __init__(self, label, elements, coords):
        """
        Create a Vertex object.
        Args:
            label (None/str): Label of vertex. This is used to draw a diagram.
            elements (List of str/List of Element(pmg) ):
                Elemental elements of coords, like ["Sr", "Ti", "O"]
            coords (numpy.array/list of float):
                Coordinates of vertex in chemical potential diagram.
        """
        # for 1 unary system, coords may be float
        if not isinstance(coords, Iterable):
            coords = [coords]
        self.label = label
        self._elem_coords = OrderedDict()
        for elem, coord in zip(elements, coords):
            self._elem_coords[elem] = coord

    @property
    def elements(self):
        """
            (list of str) Elemental elements of the composition vector.
        """
        return np.array([ec[0] for ec in self._elem_coords.items()])

    @property
    def coords(self):
        """
            (numpy.array) Coordinates of vertex.
        """
        return np.array([ec[1] for ec in self._elem_coords.items()])

    def as_dict(self):
        d = {"label": self.label,
             "elem_coords": self._elem_coords}
        return d

    @classmethod
    def from_dict(cls, d):
        label = d["label"]
        elements = list(d["elem_coords"].keys())
        coords = list(d["elem_coords"].values())
        return cls(label, elements, coords)

    def set_elements(self, elements):
        """
        Change elements of the object.
        This also changes order of composition vector.
        If self.order is subset of input order,
        then new species of atoms are added.
        If self.order is superset of input order and
        compositions of extra elements are nearly zero,
        then those atoms are removed.
        Otherwise, raise ValueError.
        Args:
            elements (list of str): New elements list.
        """
        for e, c in self._elem_coords.items():
            # case [c_a, c_b, c_c(=nonzero)] -> [c_b, c_a]
            # should be raise error
            if e not in elements and self._elem_coords[e] > 1e-5:
                raise ValueError("Invalid elements list. "
                                 "Composition of removed element {}"
                                 "must be nearly zero ( <=1e-5 ), "
                                 "but actually {} ."
                                 .format(e, self._elem_coords[e]))
        new_elem_coords = OrderedDict()
        for e in elements:
            if e in self.elements:
                new_elem_coords[e] = self._elem_coords[e]
            else:
                new_elem_coords[e] = 0
        self._elem_coords = new_elem_coords

    def __repr__(self):
        pretty_coords = np.round(self.coords, 5)
        return ("Vertex("
                "Label: {0} , "
                "Elements: {1} , "
                "Coordinates: {2})"
                .format(self.label, "-".join(self.elements), pretty_coords))

    def __eq__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self == other.
        """
        if not isinstance(other, Vertex):
            raise TypeError("Vertex class can not be"
                            " compared with other class.")
        if any([e1 != e2 for e1, e2 in zip(self.elements, other.elements)]):
            return False
        elif any([c1 != c2 for c1, c2 in zip(self.coords, other.coords)]):
            return False
        return True

    def __ne__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self != other.
        """
        return not self == other

    def __lt__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self < other.
        """
        if not isinstance(other, Vertex):
            raise TypeError("Vertex class can not be"
                            " compared with other class.")
        if any([e1 != e2 for e1, e2 in zip(self.elements, other.elements)]):
            for e1, e2 in zip(self.elements, other.elements):
                if e1 != e2:
                    return e1 < e2
        elif any([c1 != c2 for c1, c2 in zip(self.coords, other.coords)]):
            for c1, c2 in zip(self.coords, other.coords):
                if c1 != c2:
                    return c1 < c2
        else:  # all elements are strictly the same.
            return False

    def __gt__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self > other.
        """
        return (not self < other) and self != other

    def __ge__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self >= other.
        """
        return not self < other

    def __le__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self > other.
        """
        return not self > other

    def almost_equal(self, other, tol=1e-5):
        """
        Check if two equilibrium_points is almost same.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
            tol (float): Tolerance for numeric error of coordinates.
        Returns (bool): If self almost equals to other.
        """
        if not isinstance(other, Vertex):
            raise TypeError("Vertex class can not be"
                            " compared with {0} type.".format(type(other)))
        if any([e1 != e2 for e1, e2 in zip(self.elements, other.elements)]):
            return False
        elif any([abs(c1 - c2) > tol
                  for c1, c2 in zip(self.coords, other.coords)]):
            return False
        return True


class VertexOnBoundary(Vertex):

    def __init__(self, label, elements, coords, boundary_range,
                 boundary_range_limit):
        super().__init__(label, elements, coords)
        self._boundary_range = boundary_range
        self._boundary_range_limit = boundary_range_limit
        self._element_of_boundary_coords = {}
        for element, coord in zip(elements, coords):
            if np.abs(coord - boundary_range) < 1e-6:
                self._element_of_boundary_coords[element] = True
            else:
                self._element_of_boundary_coords[element] = False

    def set_boundary_range(self, boundary_range):
        if boundary_range > self._boundary_range_limit:
            raise ValueError("Boundary range is shallower than deepest vertex.")
        for element in self._element_of_boundary_coords:
            if self._element_of_boundary_coords[element]:
                self._elem_coords[element] = boundary_range

    @property
    def boundary_range_limit(self):
        return self._boundary_range_limit


class VerticesList(list):

    _TYPE_ERROR_MESSAGE = "VerticesArray must be contains only vertex."

    def __init__(self, *args, **kw):
        list.__init__(self, *args, **kw)
        for v in self:
            if not isinstance(v, Vertex):
                raise TypeError(VerticesList._TYPE_ERROR_MESSAGE)

    def append(self, vertex):
        if not isinstance(vertex, Vertex):
            raise TypeError(VerticesList._TYPE_ERROR_MESSAGE)
        list.append(self, vertex)

    def extend(self, vertices):
        if isinstance(vertices, list) and \
                all([isinstance(v, Vertex) for v in vertices]):
            list.extend(self, vertices)
        else:
            raise TypeError(VerticesList._TYPE_ERROR_MESSAGE)

    def __add__(self, vertices):
        return VerticesList(list.__add__(self, vertices))

    def __setitem__(self, key, vertex):
        if not isinstance(vertex, Vertex):
            raise TypeError(VerticesList._TYPE_ERROR_MESSAGE)
        list.__setitem__(self, key, vertex)

    def set_elements(self, order):
        for i in range(len(self)):
            self[i].set_elements(order)

    def sorted_to_loop_in_3d(self):
        """
        Args:
            sort indices of vertex to loop (for drawing diagram)
        """
        if any([len(v.coords) != 3 for v in self]):
            raise ValueError("This function can be applied to 3D vector,"
                             "But input vectors contain non-3D vector.")
        vertices_coords = [v.coords for v in self]
        mean = np.zeros(len(self[0].coords))
        for v in vertices_coords:
            mean += v / len(self)
        from_mean = [v - mean for v in vertices_coords]
        n = np.cross(from_mean[0], from_mean[1])
        n = n / np.linalg.norm(n)

        def angle_between_v0(index):
            """
            Args:
                index (int): index of vertices_coords.
            Returns (float):
                angle between from_mean[index] and from_mean[0]
            """
            v0 = from_mean[0]
            v = from_mean[index]
            v0 = v0 / np.linalg.norm(v0)
            v = v / np.linalg.norm(v)
            if np.linalg.norm(v - v0) < 1e-10:
                return 0
            dot = np.dot(v0, v)
            det = v0[0]*v[1]*n[2] + v[0]*n[1]*v0[2] + n[0]*v0[1]*v[2] - \
                v0[2]*v[1]*n[0] - v[2]*n[1]*v0[0] - n[2]*v0[1]*v[0]
            angle = np.arctan2(det, dot)
            # For easy debug
            angle = np.rad2deg(angle)
            return angle
        indices = [i for i in range(len(vertices_coords))]
        indices.sort(key=angle_between_v0)
        return VerticesList([self[i] for i in indices])

    def clear_label(self):
        """
            Clear labels of all equilibrium_points.
        """
        for i, _ in enumerate(self):
            self[i].label = None

    @property
    def is_labeled(self):
        """
            (bool) Is equilibrium_points are labeled.
        """
        return any([v.label for v in self])

    def set_alphabetical_label(self):
        """
            Label equilibrium_points alphabetically.
        """
        name_list = list(string.ascii_uppercase)
        count = 0
        for i in range(len(self)):
            if isinstance(self[i], VertexOnBoundary):
                continue
            self[i].label = name_list[count]
            count += 1

    def get_indices_and_vertices(self, label):
        """
        Find object of Compound from self by name(str) of compound.
        Args:
            label (str):
        Returns (None/VerticesList):
            Matched vertex data from label. If no compounds match, return None.
        """
        matched_list = [(i, v) for i, v in enumerate(self)
                        if v.name == label]
        if len(matched_list) == 0:
            return None
        return VerticesList(matched_list)

    def set_boundary_range(self, boundary_range):
        """
            Change coordinates of equilibrium_points on boundary and range of drawing.
        """
        for v in self:
            if isinstance(v, VertexOnBoundary):
                v.set_boundary_range(boundary_range)

    @property
    def boundary_range_limit(self):
        limit = None
        for v in self:
            if isinstance(v, VertexOnBoundary):
                limit = v.boundary_range_limit
        if limit is None:
            raise("All equilibrium_points are not on boundary."
                  "It is unexpected error."
                  "Please contact author.")
        return limit

    # TODO: document, unittest
    def almost_equal(self, other, tol=1e-5):
        """
        Check if this object and other object almost equal.
        If only order of elements like ([Ca, O] and [O, Ca]),
        judges they are same.
        Args:
            other (list of Compound or CompoundsList) : Compared object.
            tol (float) : Numerical tolerance. (Of energy, composition)
        Returns (bool) : If this object and other object almost equal.

        """
        if len(self) != len(other):
            return False
        else:
            for v1, v2 in zip(sorted(self), sorted(other)):
                if not v1.almost_equal(v2, tol=tol):
                    return False
        return True


