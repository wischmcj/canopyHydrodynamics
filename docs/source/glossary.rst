Glossary
===========================
.. glossary::
   :sorted:
  Angle
    (angle) Arctan(rise/run). Measured in radians. This is the angle of the branch with the XY plane. This angle is therefore <0 when the branch is tilted away from the tree's trunk and >0 when the branch tilted towards the trunk.
  Branch Order
    (branch_order) The number of junctions away from the trees' trunk. Trunk cyls have order 0, branches sprouting directly from the trunk have order 1, branches forking from order 1 branches have order 2 etc.
  Diameter at Breast Height
    (DBH) The diameter of the tree at a height of 1.35 meters.
  File Name
    (fileName) Can be thought of as the tree name; one files = one full set of cylinders for a given tree.
  Maximum Branch Order
    (maxBo, max_bo) The maximum branch order present in the tree   (in the cyl data).
  Maximum X Coordinate
    (X_max) The maximum X coordinate of the tree.
  Maximum Y Coordinate
    (Y_max) The maximum Y coordinate of the tree.
  Maximum Z Coordinate
    (Z_max) The maximum Z coordinate of the tree.
  Minimum X Coordinate
    (X_min) The minimum X coordinate of the tree.
  Minimum Y Coordinate
    (Y_min) The minimum Y coordinate of the tree.
  Minimum Z Coordinate
    (Z_min) The minimum Z coordinate of the tree.
  Number of Drip Points
    (num_drip_points) The number of points where water drips from the tree.
  Order One Angle Average
    (Order_one_angle_avg) The average angle of the branches with order one.
  Order One Angle Standard Deviation
    (Order_one_angle_std) The standard deviation of the angles of the branches with order one.
  Order Three Angle Average
    (Order_three_angle_avg) The average angle of the branches with order three.
  Order Two Angle Average
    (Order_two_angle_avg) The average angle of the branches with order two.
  Order Two Angle Standard Deviation
    (Order_two_angle_std) The standard deviation of the angles of the branches with order two.
  Order Zero Angle Average
    (Order_zero_angle_avg) The average angle of the branches with order zero.
  Order Zero Angle Standard Deviation
    (Order_zero_angle_std) The standard deviation of the angles of the branches with order zero.
  Projected Surface Area
    (total_psa) Total projected area for all cylinders of the tree.
  Projected Surface Area With Overlap
    (PsaWOverlap, psa_w_overlap) The sum of the projected area of each cylinder calculated independently. This is useful as it is used to exploit the fact that total_psa = psaWOverlap -   (the sum of overlapping areas of the projected cylinders). This allows us to get an estimate of what portion of the canopy is shaded.
  Stem Hull Area
    (stem_hull_area) Area of the alpha shape   (a tight boundary enclosing the stem) calculated for the tree.
  Stem Hull Boundary
    (stem_hull_boundary) The length of the perimeter of the stem alpha shape.
  Stem Projected Surface Area
    (stem_psa) Projected area for all cylinders of the stem.
  Stem Projected Surface Area With Overlap
    (stem_psa_w_overlap) The sum of the projected area of each stem cylinder calculated independently, including overlaps.
  Stem Surface Area
    (stem_surface_area) The sum of the non-projected/3D surface area of stem cylinders minus the areas of the top and bottom of the cylinder.
  Total Hull Area
    (tot_hull_area) Area of the alpha shape   (a tight boundary enclosing the canopy) calculated for the tree.
  Total Hull Boundary
    (tot_hull_boundary) The length of the perimeter of the alpha shape.
  Total Shade
    (TotalShade) PsaWOverlap - total_psa. the sum total of the area where our cylinders overlap.
  Total Surface Area
    (tot_surface_area) The sum of the non-projected/3D surface area of cylinders minus the areas of the top and bottom of the cylinder.
  Top Half Shade
    (topHalfShade, top_half_shade) See Top Quarter Shade.
  Top Half Total Projected Surface Area
    (topHalfTotPsa) The projected area of all cylinders in the top 50% of the canopy.
  Top Quarter Shade
    (topQuarterShade, top_quarter_shade) Total overlap area in the top 25% of the canopy; estimates the total shaded area   (area of diffuse v. direct sunlight) in the highest 25% of the tree.
  Top Quarter Total Projected Surface Area
    (topQuarterTotPsa) The projected area of all cylinders in the top 25% of the canopy.
  Top Three Quarter Shade
    (topThreeQuarterShade, top_three_quarter_shade) See Top Quarter Shade.
  Top Three Quarter Total Projected Surface Area
    (topThreeQuarterTotPsa) The projected area of all cylinders in the top 75% of the canopy.
  Volume
    (volume) The sum of the volumes of each cylinder.