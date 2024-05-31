
def pool_get_projection(cyl, plane):
    """
    Takes in the vector (starting point), magnitude and radius that fully define a cylinder.
    Finds the projection of the cylinder on a plane

    Some linear algebra/diff eq could help us find this for an arbtrary plane.
    """

    noCirPoints = 360
    tCir = np.linspace(
        0, 2 * np.pi, noCirPoints
    )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
    a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
    b_ortho = np.sin(tCir)  # y coordinates of the points on a circle

    radius = cyl.radius

    if plane == "XY":
        magnitude = [cyl.dx, cyl.dy, cyl.dz]
        vector = [np.transpose(cyl.x), np.transpose(cyl.y), np.transpose(cyl.z)]
    elif plane == "XZ":
        magnitude = [cyl.dx, cyl.dz, cyl.dy]
        vector = [np.transpose(cyl.x), np.transpose(cyl.z), np.transpose(cyl.y)]
    else:
        magnitude = [cyl.dy, cyl.dz, cyl.dx]
        vector = [np.transpose(cyl.y), np.transpose(cyl.z), np.transpose(cyl.x)]

    delt_a = magnitude[0]
    delt_b = magnitude[1]
    delt_c = magnitude[2]
    dim_a = vector[0]
    dim_b = vector[1]
    dim_c = vector[2]
    # unit vector at base of cylinder, pointing up cylinder axis
    vNorm = np.sqrt(delt_a**2 + delt_b**2 + delt_c**2)
    aV = np.hstack((delt_a, delt_b, delt_c)) / vNorm
    bV = -aV  # unit vector looking down from top circle (but not translated)
    # function to find the z component of an orthogonal vector in 3D
    # oVz = lambda v, a, b: ((-v[0] * a - v[1] * b) / v[2])

    # initializing min max arrays+
    min_c = np.zeros_like(delt_c)
    max_c = np.zeros_like(delt_c)
    pSV = []
    projection = {
        "cyl_id": cyl.cyl_id,
        "polygon": Polygon(),
        "base_vector": (0, 0, 0),
        "anti_vector": (0, 0, 0),
        "angle": 0,
        "area": 0,
    }
    c1 = Polygon()
    bBox = Polygon()
    c2 = Polygon()
    # coord_list = []
    try:
        # for each cylinder
        if not np.isnan(dim_a[0]):
            if np.logical_and(delt_a == 0, delt_b == 0):
                pX = dim_a[0] + radius * a_ortho
                pY = dim_b[0] + radius * b_ortho
                cPS = Polygon(list(zip(pX, pY)))
                min_c = np.min(dim_c[:])
                max_c = np.max(dim_c[:])
                ang = 0
                projection = {
                    "cyl_id": cyl.cyl_id,
                    "polygon": cPS,
                    "base_vector": aV,
                    "anti_vector": bV,
                    "angle": ang,
                    "area": cPS.area,
                }
                cyl.projected_data[plane] = projection
                
                if plane == "XY":
                    cyl.xy_area = cyl.projected_data["XY"]["area"]
                return projection
            else:
                # find orthogonal vectors @ endpoints
                # Identifies corners of projected rectangle
                aVp1 = np.hstack((aV[1], -aV[0]))
                aVp2 = np.hstack((-aV[1], aV[0]))
                bVp1 = np.hstack((bV[1], -bV[0]))
                bVp2 = np.hstack((-bV[1], bV[0]))
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

                if aV[2] != 0:
                    # calculate set of orthgonal vectors using lambda function
                    # That is 360 orthogonal vectors ending at eqidistant points along
                    # a circle of radius radius with the starting point of our cylinder
                    # at is center
                    ZOrtho = (-aV[0] * a_ortho - aV[1] * b_ortho) / aV[2]
                    # unit-ify the orthgonal vectors
                    uovd = np.sqrt(a_ortho**2 + b_ortho**2 + ZOrtho**2)
                    uov = (
                        np.hstack((a_ortho[:, None], b_ortho[:, None], ZOrtho[:, None]))
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
                    try:
                        c1c = list(
                            zip(
                                [0 if np.isnan(x) else x for x in xaC],
                                [0 if np.isnan(y) else y for y in yaC],
                            )
                        )
                        c2c = list(
                            zip(
                                [0 if np.isnan(x) else x for x in xbC],
                                [0 if np.isnan(y) else y for y in ybC],
                            )
                        )

                        # coord_list.extend(c1c)
                        # coord_list.extend(c2c)
                        c1 = Polygon(c1c)
                        c2 = Polygon(c2c)
                    except Exception as e:
                        log.debug(
                            f"Error creating circular portions of the projections {e}"
                        )

                # assemble total package
                rX = np.vstack((x1, x2, x3, x4))
                rY = np.vstack((y1, y2, y3, y4))
                # test for circle parts in polygon
                try:
                    bBoxc = list(
                        zip(
                            [0 if np.isnan(x) else x for x in rX],
                            [0 if np.isnan(y) else y for y in rY],
                        )
                    )
                    # coord_list.extend(bBoxc)
                    bBox = Polygon(bBoxc)
                    partsPS = [c1, bBox, c2]
                except:
                    log.debug(
                        f"Error creating rectangular portion of the projection: vectors:{vector} magnitudes:{magnitude}"
                    )
                try:
                    if np.max([poly_part.area for poly_part in partsPS]) > 0:
                        to_union = [poly for poly in partsPS if poly.area > 0]
                        cPS = unary_union([part for part in to_union])
                        # cPSc = Polygon(coord_list)
                except:
                    print(np.any(np.isnan(xaC)))
                    log.debug("Error unioning projection polygons ")
                # get angle away from plane projected on to
                run = math.sqrt(delt_b**2 + delt_a**2)
                rise = delt_c
                if run == 0:
                    slope = 1  # straightDown e.g. is in flow
                else:
                    slope = rise / run
                ang = np.arctan(slope)
                projection = {
                    "cyl_id": cyl.cyl_id,
                    "polygon": cPS,
                    "base_vector": aV,
                    "anti_vector": bV,
                    "angle": ang,
                    "area": cPS.area,
                }
        
                cyl.projected_data[plane] = projection

                if plane == "XY":
                    cyl.xy_area = cyl.projected_data["XY"]["area"]
                return projection

        else:
            cyl.projected_data[plane] = projection
            log.debug("dim_a[0] is null, unable to project")
            cyl.xy_area = cyl.projected_data["XY"]["area"]
            return projection
    except UnboundLocalError:
        log.debug(
            f"UnboundLocalError: vector : {vector} magnitude: {magnitude} radius: {radius}"
        ) 
        log.error(
            f"Erroring cyl id {cyl.cyl_id}"
        )
