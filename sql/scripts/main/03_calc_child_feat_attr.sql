-- Calculates surface area, tilt, azimuth, and height for each surface polygon.
--
-- Normal orientation:
--   RoofSurface  — flipped so nz > 0 (upward-facing) when the raw cross-product
--                  gives a downward normal (nz < 0).  The horizontal component of
--                  a pitched roof has a large nz contribution, so the nz-sign check
--                  is a reliable proxy for the full normal direction.  Applying the
--                  horizontal dot-product test (used for walls below) is NOT used
--                  here because for low-to-moderate tilt surfaces the horizontal
--                  normal component is small, making the dot product sign sensitive
--                  to where the surface centroid falls relative to the interior point,
--                  which introduces random azimuth errors across ±90°–180°.
--   WallSurface  — all three components flipped so the horizontal component points
--                  OUTWARD, determined by dot-product against the vector from the
--                  building's GroundSurface interior point to the wall-surface centroid.
--                  This corrects for vertex-ordering (CW vs CCW) differences across
--                  datasets; without it the azimuth of a north-facing wall may read
--                  as south (180° error).
--   GroundSurface — orientation unchanged; azimuth is suppressed (see below).
--
-- Reference point: ST_PointOnSurface of the collected GroundSurface polygons,
--   not ST_Centroid.  For L-shaped or U-shaped buildings, the centroid can fall
--   outside the footprint polygon, giving a wrong inward/outward determination.
--   ST_PointOnSurface is guaranteed to lie inside the polygon regardless of shape.
--
-- Azimuth convention:
--   Meaningful for near-vertical surfaces (walls, sloped roofs).
--   Set to -1 for near-horizontal surfaces (|nz| > 0.985, ≈ tilt > 80°) where
--   atan2(ny, nx) becomes ill-defined because both components approach zero.

WITH new_buildings AS (
  SELECT DISTINCT building_feature_id
  FROM {city2tabula_schema}.{lod_schema}_child_feature_geom_dump
  WHERE building_feature_id IN {building_ids}
    AND building_feature_id NOT IN (
      SELECT DISTINCT building_feature_id
      FROM {city2tabula_schema}.{lod_schema}_child_feature_surface
    )
),

-- Interior reference point from GroundSurface polygons.
-- ST_PointOnSurface on the collected geometry guarantees the returned point
-- lies on (and therefore inside) one of the ground polygons, even for
-- non-convex (L-shaped, U-shaped) footprints where ST_Centroid can fall
-- outside the polygon.  ST_Collect (not ST_Union) is used here: Collect
-- simply aggregates geometries with no merge/dissolve work, keeping the
-- cost identical to the original centroid CTE and avoiding lock contention
-- from long-running reads during parallel batch processing.
-- LEFT-joined so buildings without a GroundSurface still proceed (no flip applied).
building_interior_pts AS (
  SELECT
    building_feature_id,
    ST_PointOnSurface(ST_Collect(ST_Force2D(geom))) AS interior_pt
  FROM {city2tabula_schema}.{lod_schema}_child_feature_geom_dump
  WHERE building_feature_id IN (SELECT building_feature_id FROM new_buildings)
    AND classname = 'GroundSurface'
  GROUP BY building_feature_id
),

dumped AS (
  SELECT
    gd.id,
    gd.child_row_id,
    gd.building_feature_id,
    gd.surface_feature_id,
    gd.objectclass_id,
    gd.classname,
    gd.geom AS valid_geom,
    ST_IsPlanar(gd.geom) AS is_planar
  FROM {city2tabula_schema}.{lod_schema}_child_feature_geom_dump gd
  INNER JOIN new_buildings nb ON gd.building_feature_id = nb.building_feature_id
),
points AS (
  SELECT
    *,
    (ST_DumpPoints(valid_geom)).geom AS point_geom,
    (ST_DumpPoints(valid_geom)).path[2] AS pt_idx
  FROM dumped
),
all_points AS (
  SELECT
    id,
    child_row_id,
    building_feature_id,
    surface_feature_id,
    objectclass_id,
    classname,
    valid_geom,
    is_planar,
    ARRAY_AGG(point_geom ORDER BY pt_idx) AS all_pts
  FROM points
  WHERE is_planar
  GROUP BY id, child_row_id, building_feature_id, surface_feature_id, objectclass_id, classname, valid_geom, is_planar
  HAVING COUNT(*) >= 3
),
point_combinations AS (
  SELECT
    id,
    child_row_id,
    building_feature_id,
    surface_feature_id,
    objectclass_id,
    classname,
    valid_geom,
    is_planar,
    all_pts[i] AS p1,
    all_pts[j] AS p2,
    all_pts[k] AS p3,
    i, j, k
  FROM all_points,
    LATERAL generate_series(1, LEAST(array_length(all_pts, 1), 10)) i,
    LATERAL generate_series(i + 1, LEAST(array_length(all_pts, 1), 10)) j,
    LATERAL generate_series(j + 1, LEAST(array_length(all_pts, 1), 10)) k
),
vectors AS (
  SELECT *,
    ST_X(p2) - ST_X(p1) AS a_x,
    ST_Y(p2) - ST_Y(p1) AS a_y,
    ST_Z(p2) - ST_Z(p1) AS a_z,
    ST_X(p3) - ST_X(p1) AS b_x,
    ST_Y(p3) - ST_Y(p1) AS b_y,
    ST_Z(p3) - ST_Z(p1) AS b_z
  FROM point_combinations
),
cross_products AS (
  SELECT *,
    (a_y * b_z - a_z * b_y) AS n_x,
    (a_z * b_x - a_x * b_z) AS n_y,
    (a_x * b_y - a_y * b_x) AS n_z,
    sqrt((a_y * b_z - a_z * b_y)^2 +
         (a_z * b_x - a_x * b_z)^2 +
         (a_x * b_y - a_y * b_x)^2) AS cross_magnitude
  FROM vectors
),
-- Planar surfaces: existing 3-point cross-product approach (fast, exact for flat polygons).
triplet_normals AS (
  SELECT DISTINCT ON (id, child_row_id, building_feature_id, surface_feature_id)
    id,
    child_row_id,
    building_feature_id,
    surface_feature_id,
    objectclass_id,
    classname,
    valid_geom,
    is_planar,
    p1, p2, p3,
    a_x, a_y, a_z,
    b_x, b_y, b_z,
    n_x, n_y, n_z,
    cross_magnitude
  FROM cross_products
  WHERE cross_magnitude > 1e-10
  ORDER BY id, child_row_id, building_feature_id, surface_feature_id, i, j, k
),

-- Non-planar surfaces: Newell's method — accumulates cross-product contributions
-- from every consecutive edge pair around the polygon ring, giving the best-fit
-- normal equivalent to least-squares plane fitting.  Identical to the triplet
-- result for planar polygons; robust for warped (non-planar) CityGML faces where
-- different 3-point triplets yield inconsistent normals.
newell_pts AS (
  SELECT
    id, child_row_id, building_feature_id, surface_feature_id,
    objectclass_id, classname, valid_geom,
    point_geom,
    LEAD(point_geom) OVER (PARTITION BY id ORDER BY pt_idx) AS next_pt
  FROM points
  WHERE NOT is_planar
),
newell_normals AS (
  SELECT
    id, child_row_id, building_feature_id, surface_feature_id,
    objectclass_id, classname, valid_geom,
    FALSE             AS is_planar,
    NULL::geometry    AS p1,     NULL::geometry AS p2,     NULL::geometry AS p3,
    NULL::float8      AS a_x,    NULL::float8   AS a_y,    NULL::float8   AS a_z,
    NULL::float8      AS b_x,    NULL::float8   AS b_y,    NULL::float8   AS b_z,
    n_x, n_y, n_z,
    sqrt(n_x * n_x + n_y * n_y + n_z * n_z) AS cross_magnitude
  FROM (
    SELECT
      id, child_row_id, building_feature_id, surface_feature_id,
      objectclass_id, classname, valid_geom,
      SUM((ST_Y(point_geom) - ST_Y(next_pt)) * (ST_Z(point_geom) + ST_Z(next_pt))) AS n_x,
      SUM((ST_Z(point_geom) - ST_Z(next_pt)) * (ST_X(point_geom) + ST_X(next_pt))) AS n_y,
      SUM((ST_X(point_geom) - ST_X(next_pt)) * (ST_Y(point_geom) + ST_Y(next_pt))) AS n_z
    FROM newell_pts
    WHERE next_pt IS NOT NULL
    GROUP BY id, child_row_id, building_feature_id, surface_feature_id,
             objectclass_id, classname, valid_geom
    HAVING COUNT(*) >= 2
  ) sums
  WHERE sqrt(n_x * n_x + n_y * n_y + n_z * n_z) > 1e-10
),
normals AS (
  SELECT * FROM triplet_normals
  UNION ALL
  SELECT * FROM newell_normals
),

-- Compute surface_flip (+1 or -1) for WallSurface and RoofSurface.
-- dot = (surf_centroid_2d - building_interior_pt) · (nx, ny)
--   >= 0  → horizontal normal already outward  → flip = +1
--    < 0  → horizontal normal points inward    → flip = -1
-- Applied to WallSurface (all three components) and RoofSurface (nx/ny only;
-- nz is handled independently via ABS to guarantee upward-facing).
-- COALESCE to +1 when no GroundSurface interior point exists (no flip, keep as-is).
oriented AS (
  SELECT
    n.*,
    COALESCE(
      CASE
        WHEN (
            (ST_X(ST_Centroid(ST_Force2D(n.valid_geom))) - ST_X(bp.interior_pt))
              * (n.n_x / NULLIF(n.cross_magnitude, 0))
          + (ST_Y(ST_Centroid(ST_Force2D(n.valid_geom))) - ST_Y(bp.interior_pt))
              * (n.n_y / NULLIF(n.cross_magnitude, 0))
        ) >= 0 THEN 1.0
        ELSE -1.0
      END,
      1.0
    ) AS surface_flip
  FROM normals n
  LEFT JOIN building_interior_pts bp ON bp.building_feature_id = n.building_feature_id
),

final AS (
  SELECT
    id,
    child_row_id,
    building_feature_id,
    surface_feature_id,
    objectclass_id,
    classname,
    valid_geom,
    is_planar,
    p1, p2, p3,
    a_x, a_y, a_z,
    b_x, b_y, b_z,
    cross_magnitude AS norm_len,
    -- nx: RoofSurface → flip all three components when nz < 0 (downward normal);
    --     WallSurface → apply outward surface_flip; others → as-is
    CASE
      WHEN classname = 'RoofSurface' AND (n_z / NULLIF(cross_magnitude, 0)) < 0
        THEN -(n_x / NULLIF(cross_magnitude, 0))
      WHEN classname = 'WallSurface'
        THEN surface_flip * (n_x / NULLIF(cross_magnitude, 0))
      ELSE n_x / NULLIF(cross_magnitude, 0)
    END AS nx,
    -- ny
    CASE
      WHEN classname = 'RoofSurface' AND (n_z / NULLIF(cross_magnitude, 0)) < 0
        THEN -(n_y / NULLIF(cross_magnitude, 0))
      WHEN classname = 'WallSurface'
        THEN surface_flip * (n_y / NULLIF(cross_magnitude, 0))
      ELSE n_y / NULLIF(cross_magnitude, 0)
    END AS ny,
    -- nz: RoofSurface → flip when nz < 0 (guarantees upward-facing);
    --     WallSurface → apply surface_flip (nz ≈ 0 for vertical walls,
    --                   negligible effect on tilt);
    --     others → as-is
    CASE
      WHEN classname = 'RoofSurface' AND (n_z / NULLIF(cross_magnitude, 0)) < 0
        THEN -(n_z / NULLIF(cross_magnitude, 0))
      WHEN classname = 'WallSurface'
        THEN surface_flip * (n_z / NULLIF(cross_magnitude, 0))
      ELSE n_z / NULLIF(cross_magnitude, 0)
    END AS nz,
    n_x, n_y, n_z,
    cross_magnitude
  FROM oriented
)
INSERT INTO {city2tabula_schema}.{lod_schema}_child_feature_surface (
    id,
    building_feature_id,
    surface_feature_id,
    objectclass_id,
    classname,
    surface_area,
    surface_area_unit,
    tilt,
    tilt_unit,
    azimuth,
    azimuth_unit,
    is_valid,
    is_planar,
    child_row_id,
    height,
    height_unit,
    geom
)
SELECT
    gen_random_uuid() AS id,
    building_feature_id,
    surface_feature_id,
    objectclass_id,
    classname,
    CASE
      WHEN objectclass_id IN (709, 710, 712)
        THEN {city2tabula_schema}.surface_area_corrected_geom(valid_geom, nx, ny, nz)
      ELSE NULL
    END AS surface_area,
    'sqm' AS surface_area_unit,
    DEGREES(ASIN(ABS(nz))) AS tilt,
    'degrees' AS tilt_unit,
    -- Azimuth: suppress for near-horizontal surfaces (|nz| > 0.985, tilt > ~80°)
    -- where atan2(ny, nx) is ill-defined. Walls (nz ≈ 0) always get a valid azimuth.
    CASE
      WHEN ABS(nz) > 0.985 THEN -1
      ELSE MOD((450.0 - degrees(atan2(ny::numeric, nx::numeric)))::numeric + 360.0, 360.0)
    END AS azimuth,
    'degrees' AS azimuth_unit,
    ST_IsValid(valid_geom) AS is_valid,
    is_planar,
    child_row_id,
    (ST_ZMax(valid_geom) - ST_ZMin(valid_geom)) AS height,
    'm',
    valid_geom AS geom
FROM final;
