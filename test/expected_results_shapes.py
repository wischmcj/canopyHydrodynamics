from __future__ import annotations

from shapely import Polygon

example_coords = [(1, 2), (3, 5), (5, 7), (7, 9), (9, 11)]

start_poly_coords = [
    (0.75, 0.5),
    (1, 0),
    (0.5, 0.25),
    (0, 0),
    (0.25, 0.5),
    (0, 1),
    (0.5, 0.75),
    (1, 1),
    (0.75, 0.5),
]
star_poly = Polygon(start_poly_coords)

small_tree_wateshed_poly = Polygon(
    [
        (-0.497453, 3.394286),
        (-0.35321, 3.349167),
        (-0.211633, 3.262773),
        (0.027972, 3.048441),
        (0.140587, 2.959396),
        (0.236321, 2.879393),
        (0.263392, 2.859888),
        (0.975505, 2.713707),
        (1.1931, 3.021168),
        (1.220536, 3.144481),
        (1.685335, 3.381566),
        (1.787627, 4.043691),
        (1.926852, 4.01599),
        (2.084833, 3.335167),
        (2.113163, 3.284466),
        (2.065395, 2.640083),
        (2.150888, 2.547923),
        (2.23189, 2.055982),
        (2.326885, 1.930073),
        (2.335058, 1.887999),
        (2.336527, 1.776998),
        (2.247193, 1.623578),
        (2.067047, 1.434109),
        (1.806616, 1.636505),
        (1.455036, 2.239273),
        (1.047176, 2.391643),
        (0.88521, 2.500038),
        (0.825925, 2.52269),
        (0.727159, 2.562484),
        (0.507174, 2.584031),
        (0.456338, 2.624624),
        (0.296308, 2.699969),
        (-0.299115, 2.537844),
        (-0.788138, 2.866206),
        (-0.788044, 3.04045),
        (-0.788937, 3.126915),
        (-0.639431, 3.354953),
        (-0.581124, 3.378361),
        (-0.497453, 3.394286),
    ]
)


funky_squares_overlap_areas = {
    "bottom": {
        "sum_area": 1.0625,
        "effective_area": 1.0,
        "internal_overlap": 0.0625,
        "overlap_with_previous": 1.0,
    },
    "mid": {
        "sum_area": 2.0625,
        "effective_area": 2.0,
        "internal_overlap": 0.0625,
        "overlap_with_previous": 2.0,
    },
    "top": {
        "sum_area": 3.0625,
        "effective_area": 2.75,
        "internal_overlap": 0.3125,
        "overlap_with_previous": 2.75,
    },
}

small_tree_overlap = {
    25: {
        "sum_area": 0.03185684182767693,
        "effective_area": 0.029894983741628172,
        "internal_overlap": 0.0019618580860487587,
        "overlap_with_previous": 0.029894983741628172,
    },
    50: {
        "sum_area": 0.029130503264231784,
        "effective_area": 0.027960517250591167,
        "internal_overlap": 0.0011699860136406177,
        "overlap_with_previous": 0.027960517250591167,
    },
    75: {
        "sum_area": 0.029862028471880096,
        "effective_area": 0.02854408983609465,
        "internal_overlap": 0.0013179386357854463,
        "overlap_with_previous": 0.028544089836094647,
    },
}

squares_projection_overlap = {
    "bottom": {
        "sum_area": 8.0,
        "effective_area": 7.0,
        "internal_overlap": 1.0,
        "overlap_with_previous": 7.0,
    },
    "mid": {
        "sum_area": 8.0,
        "effective_area": 7.0,
        "internal_overlap": 1.0,
        "overlap_with_previous": 7.0,
    },
    "top": {
        "sum_area": 4.0,
        "effective_area": 4.0,
        "internal_overlap": 0.0,
        "overlap_with_previous": 4.0,
    },
}
