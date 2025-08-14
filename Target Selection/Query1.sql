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
    f3.gr2i2_offset,
    f3.r2i2z_offset,
    f3.i2zy_offset,
    f3.ra     AS ra_patch,
    f3.dec    AS dec_patch,
    f4.g_mag_offset,
    f4.r_mag_offset,
    f4.i_mag_offset,
    f4.z_mag_offset,
    f4.y_mag_offset
FROM pdr3_dud_dud_rev.forced AS f
LEFT JOIN pdr3_dud_rev.forced2 AS f2
  USING (object_id)
LEFT JOIN pdr3_dud_rev.patch_qa AS f3
  ON f.patch = f3.patch
LEFT JOIN pdr3_dud_rev.stellar_sequence_offsets AS f4
  ON f3.patch = f4.patch
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
