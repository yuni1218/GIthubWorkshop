-- Stack object-level science rows with offset rows (no LEFT JOIN)
WITH science AS (
  SELECT
      f.object_id,
      f.ra      AS ra_f,
      f.dec     AS dec_f,
      f.patch,
      f2.g_psfflux_mag,
      f2.r_psfflux_mag,
      f2.i_psfflux_mag,
      f2.z_psfflux_mag,
      f2.y_psfflux_mag,
      f.a_g,
      f.a_r,
      f.a_i,
      f.a_z,
      f.a_y
  FROM pdr3_dud_rev.forced AS f
  JOIN pdr3_dud_rev.forced2 AS f2
    USING (object_id)
  WHERE
      f.isprimary IS TRUE
      AND (
        tractSearch(f.object_id, 10054, 10056)
        OR tractSearch(f.object_id, 9812, 9814)
        OR tractSearch(f.object_id, 9869, 9572)
      )
      AND f2.i_psfflux_mag < 22.5
      AND f.i_extendedness_flag IS FALSE
      AND f.r_extendedness_flag IS FALSE
      AND f.z_extendedness_flag IS FALSE
      AND f.y_extendedness_flag IS FALSE
      AND f.g_extendedness_flag IS FALSE
      AND f.i_pixelflags_saturatedcenter IS FALSE
      AND f.i_pixelflags_interpolatedcenter IS FALSE
      AND f2.i_sdsscentroid_flag IS FALSE
      AND f2.g_psfflux_flag IS FALSE
      AND f2.r_psfflux_flag IS FALSE
      AND f2.i_psfflux_flag IS FALSE
      AND f2.z_psfflux_flag IS FALSE
      AND f2.y_psfflux_flag IS FALSE
)
SELECT
    s.object_id,
    s.ra_f,
    s.dec_f,
    s.patch,
    s.g_psfflux_mag,
    s.r_psfflux_mag,
    s.i_psfflux_mag,
    s.z_psfflux_mag,
    s.y_psfflux_mag,
    s.a_g,
    s.a_r,
    s.a_i,
    s.a_z,
    s.a_y,
    NULL::DOUBLE PRECISION AS g_mag_offset,
    NULL::DOUBLE PRECISION AS r_mag_offset,
    NULL::DOUBLE PRECISION AS i_mag_offset,
    NULL::DOUBLE PRECISION AS z_mag_offset,
    NULL::DOUBLE PRECISION AS y_mag_offset
FROM science AS s

UNION ALL

SELECT
    NULL::BIGINT            AS object_id,
    NULL::DOUBLE PRECISION  AS ra_f,
    NULL::DOUBLE PRECISION  AS dec_f,
    NULL::TEXT              AS patch,
    NULL::DOUBLE PRECISION  AS g_psfflux_mag,
    NULL::DOUBLE PRECISION  AS r_psfflux_mag,
    NULL::DOUBLE PRECISION  AS i_psfflux_mag,
    NULL::DOUBLE PRECISION  AS z_psfflux_mag,
    NULL::DOUBLE PRECISION  AS y_psfflux_mag,
    NULL::DOUBLE PRECISION  AS a_g,
    NULL::DOUBLE PRECISION  AS a_r,
    NULL::DOUBLE PRECISION  AS a_i,
    NULL::DOUBLE PRECISION  AS a_z,
    NULL::DOUBLE PRECISION  AS a_y,
    f4.g_mag_offset,
    f4.r_mag_offset,
    f4.i_mag_offset,
    f4.z_mag_offset,
    f4.y_mag_offset
FROM pdr3_dud_rev.stellar_sequence_offsets AS f4;
