
# @jit
# def cyl_vectors(magnitude):

# return aV, bV, a_ortho, b_ortho

# https://stackoverflow.com/questions/39822480/plotting-a-solid-cylinder-centered-on-a-plane-in-matplotlib
def vectorized_get_projection(starts: np.array(), ends: np.array(), radii:np.array()):
    return 'to do'
#     # 50 or so line segments arranged in a polygon
#     # gets us 99% accuracy in approximating the area of a circle
#     # (.5*n*r*r)sin(2*pi/n) = area of n-gon, 2*pi*r*r = area of circle
#     num_approximating_segments = 50
#     magnitude = ends - starts
#     # p0 = np.array([1, 3, 2]) starts
#     # p1 = np.array([8, 5, 9]) ends
#     # R = 5
#     #vector in direction of axis
#     # v = p1 - p0
#     vectors = ends-starts
#     #find magnitude of vector
#     mags = np.apply_along_axis(np.linalg.norm, 1, vectors)
#     #unit vectors for cyl
#     unit_vectors = vectors/mags[:,None]
#     anti_unit__vectors = -unit_vectors
#     #orthogonal vectors
#     rando = np.array([1, 0, 0])
#     if (vectors == rando).all():
#         rando = np.array([0, 1, 0])
#     urando = (rando/np.linalg.norm(rando))
#     first_ortho_vectors = np.cross(unit_vectors, urando)
#     second_ortho_vectors = np.cross(unit_vectors, first_ortho_vectors)

#     num_approx_segments = 50
#     t = np.linspace(0, mags, 50)
#     theta = np.linspace(0, 2 * np.pi, num_approx_segments)

#     rect, cir = np.meshgrid(t, theta)
#     #generate coordinates for surface

#     vx = [unit_vectors]
#     vy =
#     vz =

#     X = [end[0] + start[0] * rect + radii * np.sin(cir) * n1[0] + radii * np.cos(cir) * n2[0] ]


#     X, Y, Z = [ends[i] + starts[i] * rect +
#                 radii[:,None] * np.sin(cir) * first_ortho_vectors[i] +
#                 radii[:,None] * np.cos(cir) * second_ortho_vectors[i] for i in [0, 1, 2]]

#     if (v == not_v).all():
#         not_v = np.array([0, 1, 0])
#     #make vector perpendicular to v
#     n1 = np.cross(v, not_v)
#     #normalize n1
#     n1 /= norm(n1)
#     #make unit vector perpendicular to v and n1
#     n2 = np.cross(v, n1)
#     #surface ranges over t from 0 to length of axis and 0 to 2*pi
#     t = np.linspace(0, mag, 100)
#     theta = np.linspace(0, 2 * np.pi, 100)
#     #use meshgrid to make 2d arrays
#     t, theta = np.meshgrid(t, theta)

#     X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]


#stackoverflow.com/questions/39822480/plotting-a-solid-cylinder-centered-on-a-plane-in-matplotlib
def vectorized_def_cyl(vector, magnitude):
    dim_a = vector[0]
    dim_b = vector[1]
    dim_c = vector[2]
    delt_a = magnitude[0]
    delt_b = magnitude[1]
    delt_c = magnitude[2]

    noCirPoints = 360
    tCir = np.linspace(
        0, 2 * np.pi, noCirPoints
    )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
    a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
    b_ortho = np.sin(tCir)  # y coordinates of the points on a circle
    vNorm = np.sqrt(delt_a**2 + delt_b**2 + delt_c**2)
    # aV =  np.array([0.0,0.0,0.0]) if vNorm == 0 else np.array([delt_a, delt_b, delt_c])/ vNorm
    aV = np.array([delt_a, delt_b, delt_c]) / vNorm
    bV = -aV  # unit vector looking down from top circle (but not translated)
    # unit vector at base of cylinder, pointing up cylinder axis

    # function to find orthgonal vectors
    oVz = lambda v, a, b: ((-v[0] * a - v[1] * b) / v[2])
    # initializing min max arrays+
    pSV = []
    run = math.sqrt(delt_b**2 + delt_a**2)
    rise = delt_c
    if run == 0:
        slope = 1  # straightDown e.g. is in flow
    else:
        slope = rise / run
    ang = np.arctan(slope)
    print(f"projecting vector {aV}")
    # try:

    c1 = np.zeros((360, 2))
    # bBox = stack((
    #                     np.array([0.0,0.0,0.0,0.0]),
    #                     np.array([0.0,0.0,0.0,0.0])
    #                 ),
    #                 True
    #         )
    c2 = np.zeros((360, 2))
    bBox = np.zeros((4, 2))

    # generate coordinates for surface
    # "Tube"
    X, Y, Z = (
        p0[i] + v[i] * t + R * np.sin(theta2) * n1[i] + R * np.cos(theta2) * n2[i]
        for i in [0, 1, 2]
    )
    # "Bottom"
    X2, Y2, Z2 = (
        p0[i] + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i]
        for i in [0, 1, 2]
    )
    # "Top"
    X3, Y3, Z3 = (
        p0[i]
        + v[i] * mag
        + rsample[i] * np.sin(theta) * n1[i]
        + rsample[i] * np.cos(theta) * n2[i]
        for i in [0, 1, 2]
    )

    if not np.isnan(dim_a[0]):
        if np.logical_and(delt_a == 0.0, delt_b == 0.0):
            pX = dim_a[0] + radius * a_ortho
            pY = dim_b[0] + radius * b_ortho
            c2 = stack(np.array((pX, pY)), col=True)
            ang = 1.5708  # 90* straight up and down
        else:
            if aV[2] != 0.0:
                # calculate circular sections (only relevant if cyl is not parallel with the xy plane)
                # calculate set of orthgonal vectors using lambda function
                # That is 360 orthogonal vectors ending at eqidistant points along
                # a circle of radius radius with the starting point of our cylinder
                # at is center
                ZOrtho = oVz(aV[:], a_ortho, b_ortho)
                # unit-ify the orthgonal vectors
                uovd = np.sqrt(a_ortho**2 + b_ortho**2 + ZOrtho**2)
                # Confounded - why does removing the first three [:,None]'s below lead to non-circular projections
                # for XZ?
                uov = (
                    np.hstack(np.array((a_ortho[:, None], b_ortho[:, None], ZOrtho[:, None])))
                    / uovd[:, None]
                )
                # donot re unit-fy, you only want the horizontal component, not the
                # renormalized horizontal component
                # using only the X and Y components, find circle coods in plane of
                # interest
                xaC = dim_a[0] + uov[:, 0] * radius
                yaC = dim_b[0] + uov[:, 1] * radius
                zaC = dim_c[0] + uov[:, 2] * radius
                xbC = dim_a[1] + uov[:, 0] * radius
                ybC = dim_b[1] + uov[:, 1] * radius
                zbC = dim_c[1] + uov[:, 2] * radius

                # c1 =  np.column_stack((
                #            [0 if np.isnan(x) else x for x in xaC],
                #            [0 if np.isnan(y) else y for y in yaC],
                #         ))
                c1 = stack(
                    np.array((
                        np.array([0 if np.isnan(x) else x for x in xaC]),
                        np.array([0 if np.isnan(y) else y for y in yaC]),
                    )),
                    col=True,
                )
                print(c1)
                xC = np.array([0 if np.isnan(x) else x for x in xbC])
                yC = np.array([0 if np.isnan(y) else y for y in ybC])
                # c2 =  np.column_stack((
                #             xC,
                #             yC,
                #         ))
                c2 = stack(np.array((xC, yC)), col=True)
            # calculating the rectangular portion of the projection
            # relevant for all cyls

            # find orthogonal vectors @ endpoints
            # Identifies corners of projected rectangle
            aVp1 = np.array([aV[1], -aV[0]])
            aVp2 = np.array([-aV[1], aV[0]])
            bVp1 = np.array([bV[1], -bV[0]])
            bVp2 = np.array([-bV[1], bV[0]])
            aVp1 = aVp1 / np.linalg.norm(aVp1)
            aVp2 = aVp2 / np.linalg.norm(aVp2)
            bVp1 = bVp1 / np.linalg.norm(bVp1)
            bVp2 = bVp2 / np.linalg.norm(bVp2)
            # from each endpoint, use radius to find vertices of the rectangle
            x1 = dim_a[0] + radius * aVp1[0]
            y1 = dim_b[0] + radius * aVp1[1]
            x2 = dim_a[0] + radius * aVp2[0]
            y2 = dim_b[0] + radius * aVp2[1]
            x3 = dim_a[1] + radius * bVp1[0]
            y3 = dim_b[1] + radius * bVp1[1]
            x4 = dim_a[1] + radius * bVp2[0]
            y4 = dim_b[1] + radius * bVp2[1]

            # x1             aVp1
            # x2 =   radius* aVP2
            # x3             bVp1
            # x4             BvP2

            bBox = stack(np.array((np.array([x1, x2, x3, x4]), np.array([y1, y2, y3, y4]))), True)

        # print(".")
        # print(typeof(bBox))
        # print(".")
        # print(typeof(c1))
        # print(".")
        # print(typeof(c2))
        # print(".")
        # print(typeof(aV))
        return c1, bBox, c2, ang, aV, bV



#Cylinder Collection 


def vectorized_project_cylinders(
    self, plane: str = "XY", force_rerun: bool = False
):
    """Projects cylinders onto the specified plane"""
    if plane not in ("XY", "XZ", "YZ"):
        log.info(f"{plane}: invalid value for plane")
    # elif not force_rerun and self.projections[plane]:
    #     log.info(
    #         "cached projections exist, pass 'force_rerun=True to calculate new projections "
    #     )
    else:
        polys = []
        log.info(f"Projection into {plane} axis begun for file {self.file_name}")
        starts = np.array([cyl.vectors[plane][0] for cyl in self.cylinders])
        ends = np.array([cyl.vectors[plane][1] for cyl in self.cylinders])
        radii = np.array([cyl.radius for cyl in self.cylinders])
        vectorized_get_projection(starts, ends, radii)
        self.projections[plane] = True
        self.pSV = polys