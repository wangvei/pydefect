# -*- coding: utf-8 -*-


def get_neighbors_with_distance(structure, coords, neighbor_tolerance=1.1,
                                initial_cutoff=2, is_frac=True):
    """
    Args:
        structure (Structure):
        coords (3x1 array):
        neighbor_tolerance (float):
            get sites whose distance < min_distance * neighbor_tolerance
        initial_cutoff (float): initial cutoff distance when searching neighbors
        is_frac (bool): Is coordinate fractional, not Cartesian
    Returns:
        (list of (Site, distance))
    """

    # conversion to Cartesian
    if is_frac:
        cart_coords = structure.lattice.get_cartesian_coords(coords)
    else:
        cart_coords = coords

    self_threshold = 1e-2

    def exclude_self(distance_list):
        return [site_dist for site_dist in distance_list
                if site_dist[1] > self_threshold]

    # search neighbor
    cutoff = initial_cutoff
    sites_dist = []
    while not exclude_self(sites_dist):
        cutoff = cutoff * 1.2
        sites_dist = structure.get_sites_in_sphere(cart_coords, cutoff)

    min_distance = min(distance for _, distance in exclude_self(sites_dist))
    max_distance = min_distance * neighbor_tolerance
    candidates = structure.get_sites_in_sphere(cart_coords, max_distance)
    sites_dist = [site_dist for site_dist in candidates
                  if site_dist[1] >= min_distance]

    return sites_dist


def get_neighbors(structure, coords, neighbor_tolerance,
                  initial_cutoff=3, is_frac=True, include_self=False):
    """
    Args:
        structure (Structure):
        coords (3x1 array):
        neighbor_tolerance (float):
            get sites whose distance < min_distance * neighbor_tolerance
        initial_cutoff (float):
            initial cutoff distance when searching neighbors
        is_frac (bool):
            Whether the coordinate is in fractional instead of cartesian
        include_self (bool):
            Whether the sites whose distance < 1e-2 are included.
    Returns:
        (list of Site)
    """
    return [site for site, _ in get_neighbors_with_distance(
        structure, coords, neighbor_tolerance,
        initial_cutoff=initial_cutoff, is_frac=is_frac)]