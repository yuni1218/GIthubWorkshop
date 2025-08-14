-- wide
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
    f.a_y,
    (f.skymap_id / 10000)          AS tract,
    ((f.skymap_id % 10000) / 100)  AS patch_x,
    (f.skymap_id % 100)            AS patch_y
  
FROM pdr3_dud_rev.forced AS f
LEFT JOIN pdr3_dud_rev.forced2 AS f2
  USING (object_id)
WHERE
    f.isprimary = TRUE
    AND (
      tractSearch(f.object_id, 10054, 10056)
      OR tractSearch(f.object_id, 9812, 9814)
      OR tractSearch(f.object_id, 9869, 9572)
    )
    AND f2.i_psfflux_mag < 22.5
    AND f.i_extendedness_flag = FALSE
    AND f.r_extendedness_flag = FALSE
    AND f.z_extendedness_flag = FALSE
    AND f.y_extendedness_flag = FALSE
    AND f.g_extendedness_flag = FALSE
    AND f.i_pixelflags_saturatedcenter = FALSE
    AND f.i_pixelflags_interpolatedcenter = FALSE
    AND f2.i_sdsscentroid_flag = FALSE
    AND f2.g_psfflux_flag = FALSE
    AND f2.r_psfflux_flag = FALSE
    AND f2.i_psfflux_flag = FALSE
    AND f2.z_psfflux_flag = FALSE
    AND f2.y_psfflux_flag = FALSE;
