from __future__ import print_function
import argparse
import copy
from collections import OrderedDict

import ruamel.yaml
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.spatial import HalfspaceIntersection

from chempotdiag.compound import Compound, DummyCompoundForDiagram,\
    CompoundsList
from chempotdiag.vertex import Vertex, VertexOnBoundary, VerticesList


class ChemPotDiag:
    """
        Object for chemical potentials diagram.
    """

    def __init__(self, element_energy, stable_compounds, unstable_compounds,
                 vertices,
                 compounds_to_vertex_list, vertex_to_compounds_list):
        self._element_energy = element_energy
        self._stable_compounds = stable_compounds
        self._unstable_compounds = unstable_compounds
        self._vertices = vertices
        self._compounds_to_vertex_list = compounds_to_vertex_list
        self._vertex_to_compounds_list = vertex_to_compounds_list

    @classmethod
    def from_calculation(cls, input_compounds_array):
        """
        Create a object of ChemPot.
        Args:
            input_compounds_array (CompoundsList): Considered compounds array
        """
        compounds_array = copy.deepcopy(input_compounds_array)
        # Temporary
        compounds_array.set_elements(input_compounds_array.elements)
        compounds_array.standardize_energies()
        element_energy = compounds_array.element_energies
        dim = len(compounds_array.elements)
        if dim == 1:
            min_energy = min(compounds_array, key=lambda x: x.energy).energy
            stable_compounds \
                = CompoundsList([comp for comp in compounds_array
                                 if comp.energy == min_energy])
            unstable_compounds \
                = CompoundsList([comp for comp in compounds_array
                                 if comp.energy != min_energy])
            vertices = VerticesList([])
            for comp in compounds_array:
                vertex = Vertex(None, compounds_array.elements, comp.energy)
                vertices.append(vertex)
            compounds_to_vertex_list = list(range(len(compounds_array)))
            vertex_to_compounds_list = list(range(len(compounds_array)))
            return cls(element_energy, stable_compounds, unstable_compounds,
                       vertices,
                       compounds_to_vertex_list, vertex_to_compounds_list)
        # set boundary range
        BOUNDARY_RATE = 1.1  # can be arbitrary number > 1.0
        intersections = np.zeros((len(compounds_array), dim))
        for i, comp_dat in enumerate(compounds_array):
            e = comp_dat.energy
            c = comp_dat.composition
            for j in range(dim):
                if c[j] != 0:
                    intersections[i][j] = e / c[j]
                else:  # This case does not related to search_vertices_range.
                    intersections[i][j] = float("inf")
        search_vertices_range = np.min(intersections) * BOUNDARY_RATE
        elements = compounds_array.elements
        for e in compounds_array.elements:
            energy = search_vertices_range
            boundary \
                = DummyCompoundForDiagram.construct_boundary(elements,
                                                             e,
                                                             -energy)
            compounds_array.append(boundary)
        # Calculate HalfspaceIntersection
        eps = 1e-2
        interior_point \
            = np.array([search_vertices_range + eps] * dim)
        halfspaces = []
        for comp_dat in compounds_array:
            h = np.append(comp_dat.composition, -comp_dat.energy)
            halfspaces.append(h)
        halfspaces = np.array(halfspaces)
        # TODO: Probably there is a more simple way to make 2D np.array.
        hi = HalfspaceIntersection(halfspaces, interior_point)
        # facets_by_halfspace
        n = len(hi.halfspaces)
        facets_by_halfspace = []
        for i in range(n):
            indices = []
            for j, _ in enumerate(hi.dual_facets):
                if i in hi.dual_facets[j]:
                    indices.append(j)
            facets_by_halfspace.append(indices)
        # classify all compounds into dummy, stable, and unstable.
        stable_compounds = CompoundsList([])
        unstable_compounds = CompoundsList([])
        stable_original_index_list = []
        unstable_original_index_list = []
        # (un)stable_compounds_list[i] =
        # compounds_array[(un)stable_original_index_list[i]]
        which_vertex_on_boundary = set()
        for i, compound in enumerate(compounds_array):
            if isinstance(compound, DummyCompoundForDiagram):
                which_vertex_on_boundary |= set(facets_by_halfspace[i])
            elif facets_by_halfspace[i]:
                stable_compounds.append(compounds_array[i])
                stable_original_index_list.append(i)
            else:
                unstable_compounds.append(compounds_array[i])
                unstable_original_index_list.append(i)
        stable_compounds.set_elements(compounds_array.elements)
        unstable_compounds.set_elements(compounds_array.elements)
        # make vertices_array
        draw_criterion = min([item for item in hi.intersections.flatten()
                              if abs(item - search_vertices_range) > 1e-6])
        draw_vertices = copy.deepcopy(hi.intersections)
        draw_range = draw_criterion * 1.1
        # sometimes boundary range too small like -100
        # due to an unstable substance, then the value is changed to
        # draw_range
        vertices = VerticesList([])
        for i in range(len(draw_vertices)):
            if i in which_vertex_on_boundary:
                is_element_boundary \
                    = [abs(v - search_vertices_range) < 1e-6
                       for v in draw_vertices[i]]
                for j, flag in enumerate(is_element_boundary):
                    if flag:
                        draw_vertices[i][j] = draw_range
                        vertex = VertexOnBoundary(None,
                                                  elements,
                                                  draw_vertices[i],
                                                  draw_range,
                                                  draw_criterion)
            else:
                vertex = Vertex(None, elements, draw_vertices[i])
            vertices.append(vertex)
        # make compounds_to_vertex_list, vertex_to_compounds_list
        compounds_to_vertex_list = [l for l in facets_by_halfspace if l]
        vertex_to_compounds_list = []
        for i, l in enumerate(hi.dual_facets):
            stable_index_list \
                = [stable_original_index_list.index(j)
                   for j in l if j in stable_original_index_list]
            vertex_to_compounds_list.append(stable_index_list)
        return cls(element_energy, stable_compounds, unstable_compounds,
                   vertices,
                   compounds_to_vertex_list, vertex_to_compounds_list)

    @classmethod
    def from_file(cls, filename):
        compounds_array = CompoundsList.from_file(filename)
        return cls.from_calculation(compounds_array)

    @classmethod
    def from_vasp_calculations_files(cls,
                                     poscar_paths,
                                     output_paths,
                                     fmt="outcar"):
        """
        Args:
            poscar_paths (list of str):
            output_paths (list of str):
            fmt (str): "outcar" or "vasprun".
        Returns:
            (ChemPotDiag) ChemPotDiag object from vasp files.
        """
        compounds_list = \
            CompoundsList.from_vasp_calculations_files(poscar_paths,
                                                       output_paths,
                                                       fmt=fmt)
        return cls.from_calculation(compounds_list)

    @property
    def dim(self):
        """
            (int) Number of considered atoms.
        """
        return self.stable_compounds.dim

    @property
    def elements(self):
        """
            (list of str) Considered elements.
            Order of list is used as order of coordinates.
        """
        return self.stable_compounds.elements

    @property
    def element_energy(self):
        """
            (list of float) Element energy.
        """
        return self._element_energy

    @property
    def stable_compounds(self):
        """
            (CompoundsList) CompoundList
        """
        return copy.deepcopy(self._stable_compounds)

    @property
    def unstable_compounds(self):
        """
            (CompoundsList) CompoundList
        """
        return copy.deepcopy(self._unstable_compounds)

    @property
    def all_compounds(self):
        """
            (CompoundsList) CompoundList
        """
        return self.stable_compounds + self.unstable_compounds

    def set_elements(self, elements):
        """
        Change elements of elements. Internal composition data will be rearranged.
        Args:
            elements (list of str) :
                Elemental elements, like ["Sr", "O", "Ti"]
                Element of Pymatgen can be also accepted,
                like [Element("Sr"), ...]
        """
        self._stable_compounds.set_elements(elements)
        self._unstable_compounds.set_elements(elements)
        self._vertices.set_elements(elements)

    @property
    def vertices(self):
        """
            (VerticesList) Vertices, including vertices on drawing boundary
        """
        return self._vertices

    @property
    def equilibrium_points(self):
        """
            (VerticesList) Vertices, excluding vertices on drawing boundary
                           and only physically meaningful points are included.
        """
        return VerticesList([v for v in self._vertices
                             if not isinstance(v, VertexOnBoundary)])

    def get_neighbor_vertices(self, compound):
        """
        Search equilibrium_points on remarked compound.
        Args:
            compound (int/str/Compound):
                Index of compound or compound name or Compound object.
        Returns (VerticesList): Vertices on input compound.
        """
        if isinstance(compound, int):
            index = compound
        elif isinstance(compound, str):
            matched = self.stable_compounds.get_indices_and_compounds(compound)
            if len(matched) >= 2:
                raise ValueError\
                    ("More than one compounds matched the name {0}."
                     .format(compound))
            if matched:
                index = matched[0][0]
        elif isinstance(compound, Compound):
            index = self.stable_compounds.index(compound)
        else:
            raise TypeError("Input compound must be int or str or Compound.")
        if index is None:
            raise ValueError\
                ("{0} did not found in stable_compounds.".format(compound))
        return VerticesList([self._vertices[i]
                             for i in self._compounds_to_vertex_list[index]])

    def get_neighbor_compounds(self, vertex):
        """
        Search equilibrium_points on remarked vertex.
        Args:
            vertex (int/str/Vertex):
                Input compound (index or label or vertex).
        Returns (CompoundsList): Compounds on the vertex.
        """
        if isinstance(vertex, int):
            index = vertex
        elif isinstance(vertex, str):
            matched = self._vertices.get_indices_and_vertices(vertex)
            if len(matched) >= 2:
                raise ValueError(
                    "More than one equilibrium_points matched the name {0}.".
                    format(vertex))
            if matched:
                index = matched[0][0]
        elif isinstance(vertex, Vertex):
            index = self._vertices.index(vertex)
        else:
            raise TypeError("Input vertex must be int or str or Vertex.")

        if index is None:
            raise ValueError(
                "{0} did not found in equilibrium_points.".format(vertex))
        return CompoundsList([self.stable_compounds[i]
                              for i in self._vertex_to_compounds_list[index]])

    def get_neighbor_vertices_as_dict(self, remarked_compound, **kwargs):
        d = {"compound": remarked_compound}
        standard_energy = {str(el): float(en) for el, en
                           in zip(self.elements, self.element_energy)}
        d["standard_energy"] = standard_energy
        vertices_list =\
            VerticesList(self.get_neighbor_vertices(remarked_compound))
        if self.dim == 3:
            vertices_list.sorted_to_loop_in_3d()
        vertices_list.set_alphabetical_label()
        for v in vertices_list:
            if not isinstance(v, VertexOnBoundary):
                d[v.label] = {str(el): float(en) for el, en
                              in zip(self.elements, v.coords)}
        for k, v in kwargs.items():
            d[k] = v
        return d

    def dump_yaml(self, file_path, remarked_compound, **kwargs):
        """

        Args:
            file_path(str):
            remarked_compound(str):
            **kwargs(dict): other property, like {"comment" : "foo,var"}

        Returns:

        """
        d = self.get_neighbor_vertices_as_dict(remarked_compound, **kwargs)
        filename = file_path + "/vertices_" + remarked_compound + ".yaml"
        with open(filename, 'w') as f:
            ruamel.yaml.dump(d, f)

    def draw_diagram(self,
                     title=None,
                     save_file_name=None,
                     remarked_compound=None,
                     with_label=True,
                     draw_range=None):
        """
        Draw chemical potential diagram.
        Args:
            title (str): Title of diagram.
            save_file_name (None/str): If you will save diagram as image file,
                                       specify name of the file by this arg.
            remarked_compound (None/str): Vertices on the specified compound
                                          will be labeled.
            with_label (bool): Whether if you will display names of compounds.
            draw_range (None/float): Range to draw diagram.
                                     If none, range will be
                                     determined automatically.
        """
        if self.dim >= 4:
            # TODO: fix some parameter and plot in 3D or 2D?
            raise NotImplementedError("4D or more data can not be drawn.")

        if self.dim != 1:
            if not draw_range:
                draw_range = self._vertices.boundary_range_limit * 1.1
            self._vertices.set_boundary_range(draw_range)

        #  1D, 2D, and 3D dimension. More than 4D has not yet implemented.
        if self.dim == 1:
            ax = plt.figure().add_subplot(111)
            x = np.zeros(len(self.stable_compounds)
                         + len(self._unstable_compounds))
            y = [cd.energy
                 for cd in self.stable_compounds + self._unstable_compounds]
            ax.scatter(x, y)
            ax.set_ylabel('Chemical potential of ' + self.elements[0])
            ax.set_xlim(-1, 1)
            y_max = np.max(y)
            ax.set_ylim(-y_max * 0.1, y_max * 1.2)
            for cd in self.all_compounds:
                x_shift = -0.2
                ax.text(x_shift,
                        cd.energy,
                        cd.name,
                        size='smaller')
            if title:
                ax.set_title(title)
            else:
                ax.set_title \
                    ('Chemical potential diagram of ' + self.elements[0])

        elif self.dim == 2:
            ax = plt.figure().add_subplot(111)
            num_line = len(self.stable_compounds)
            for i, compound in enumerate(self.stable_compounds):
                vertices \
                    = VerticesList([self._vertices[j]
                                    for j in self._compounds_to_vertex_list[i]])
                vertices_coords = [v.coords for v in vertices]
                x = [v[0] for v in vertices_coords]
                y = [v[1] for v in vertices_coords]
                color = "black"
                plt.plot(x, y, color=color)
                mean = np.mean(vertices_coords, axis=0)
                x_shift = -0.2
                y_shift = -0.2
                if with_label:
                    ax.text(mean[0] + x_shift,
                            mean[1] + y_shift,
                            compound.name,
                            size='smaller',
                            zorder=num_line + i, )
                    if compound.name == remarked_compound:
                        vertices.set_alphabetical_label()
                        for j, v in enumerate(vertices):
                            if v.label:
                                ax.text(v.coords[0] + x_shift,
                                        v.coords[1] + y_shift,
                                        v.label,
                                        size="smaller",
                                        zorder=2*num_line+j,
                                        color="red",
                                        weight="bold")
            if title:
                ax.set_title(title)
            else:
                ax.set_title('Chemical potential diagram of ' + self.elements[0]
                             + " and " + self.elements[1])
            ax.set_xlabel('Chemical potential of ' + self.elements[0])
            ax.set_ylabel('Chemical potential of ' + self.elements[1])
            ax.set_xlim(draw_range, 0)
            ax.set_ylim(draw_range, 0)

        elif self.dim == 3:
            ax = plt.figure().add_subplot(111, projection='3d')
            num_plane = len(self._stable_compounds)
            for i, compound in enumerate(self._stable_compounds):
                vertices_coords = self.get_neighbor_vertices(i)
                sorted_vertices = vertices_coords.sorted_to_loop_in_3d()
                sorted_vertices_coords = [v.coords for v in sorted_vertices]
                face = Poly3DCollection([sorted_vertices_coords])
                color = [(c + 3) / 4 for c in compound.composition]
                face.set_color(color)
                face.set_edgecolor("black")
                ax.add_collection3d(face)
                mean = np.mean(sorted_vertices_coords, axis=0)
                if with_label:
                    ax.text(mean[0], mean[1], mean[2],
                            compound.name,
                            size='smaller',
                            zorder=num_plane+i,
                            ha='center',
                            va='center')
                    if compound.name == remarked_compound:
                        vertices_coords.set_alphabetical_label()
                        for j, v in enumerate(vertices_coords):
                            if v.label:
                                ax.text(v.coords[0],
                                        v.coords[1],
                                        v.coords[2],
                                        v.label,
                                        size="smaller",
                                        zorder=2 * num_plane + j,
                                        weight="bold",
                                        color="red")
            if title:
                ax.set_title(title)
            else:
                ax.set_title('Chemical potential diagram of '
                             + self.elements[0] + ", "
                             + self.elements[1] + ", and "
                             + self.elements[2])
            ax.set_xlabel('Chemical potential of ' + self.elements[0])
            ax.set_ylabel('Chemical potential of ' + self.elements[1])
            ax.set_zlabel('Chemical potential of ' + self.elements[2])
            ax.set_xlim3d(draw_range, 0)
            ax.set_ylim3d(0, draw_range)
            ax.set_zlim3d(draw_range, 0)

        ax.grid(color='b',
                alpha=0.2,
                linestyle='dashed',
                linewidth=0.5)
        plt.tight_layout()  # try
        if save_file_name:
            plt.savefig(save_file_name)
        else:
            plt.show()


def main():
    parser = argparse.ArgumentParser()

    # input
    parser.add_argument("-e", "--energy", dest="energy_file", type=str,
                        default=None,
                        help="Name of text file of energies of compounds")
    parser.add_argument("-v", "--vasp_dirs",
                        dest="vasp_dirs", type=str, nargs='+',
                        default=None,
                        help="Drawing diagram from specified directories"
                             "of vasp calculations")
    parser.add_argument("-p", "--poscar_name",
                        dest="poscar_name", type=str,
                        default="POSCAR",
                        help="Name of POSCAR, like CONTCAR, POSCAR-finish,...")
    parser.add_argument("-o", "--outcar_name",
                        dest="outcar_name", type=str,
                        default="OUTCAR",
                        help="Name of OUTCAR, like OUTCAR-finish")

    # drawing diagram
    parser.add_argument("-w", "--without_label",
                        help="Draw diagram without label.",
                        action="store_true")
    parser.add_argument("-c", "--remarked_compound",
                        dest="remarked_compound", type=str,
                        default=None,
                        help="Name of compound you are remarking."
                             "Outputted equilibrium_points are limited to "
                             "neighboring that compounds, "
                             "and those equilibrium_points are "
                             "labeled in chem_pot_diagram.")
    parser.add_argument("-d", "--draw_range",
                        dest="draw_range", type=float,
                        default=None,
                        help="Drawing range of diagram."
                             "If range is shallower than the deepest vertex,"
                             "ValueError will occur")

    # output
    parser.add_argument("-s", "--save_file",
                        dest="save_file", type=str,
                        default=None,
                        help="File name to save the drawn diagram.")
    parser.add_argument("-y", "--yaml",
                        action="store_const", const=True, default=False,
                        help="Dumps yaml of remarked_compound")

    options = parser.parse_args()
    if options.energy_file and options.vasp_dirs:
        raise ValueError("You can not specify energy_file and vasp_dirs "
                         "simultaneously.")
    if options.energy_file:
        cp = ChemPotDiag.from_file(options.energy_file)
    if options.vasp_dirs:
        poscar_paths = [d + options.poscar_name for d in options.vasp_dirs]
        outcar_paths = [d + options.outcar_name for d in options.vasp_dirs]
        cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths, outcar_paths)
    print("Energies of elements ({0}) : {1}"
          .format(cp.elements, cp.element_energy))
    #  Read options of drawing diagram from parser
    if options.remarked_compound:
        try:
            for vertex in cp.get_neighbor_vertices(options.remarked_compound):
                print(vertex)
        except ValueError:
            print("{0} is unstable. No vertex is labeled."
                  .format(options.remarked_compound))

    kwargs_for_diagram = {}
    if options.remarked_compound:
        kwargs_for_diagram["remarked_compound"] = options.remarked_compound
    if options.save_file:
        kwargs_for_diagram["save_file_name"] = options.save_file
    if options.without_label:
        kwargs_for_diagram["with_label"] = False
    if options.draw_range:
        kwargs_for_diagram["draw_range"] = options.draw_range

    if cp.dim >= 4:
        print("Currently diagram is not available for quaternary or more.")
    else:
        try:
            cp.draw_diagram(**kwargs_for_diagram)
        except ValueError:
            kwargs_for_diagram.pop("remarked_compound")
            cp.draw_diagram(**kwargs_for_diagram)

    if options.yaml:
        if options.remarked_compound is None:
            raise ValueError("remarked_compound is needed to dump yaml")
        cp.dump_yaml(".", options.remarked_compound)


if __name__ == "__main__":
    main()

