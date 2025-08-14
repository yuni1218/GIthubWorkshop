SELECT
    f.g_mag_offset,
    f.r_mag_offset,
    f.i_mag_offset,
    f.z_mag_offset,
    f.y_mag_offset,
    (f.skymap_id / 10000)          AS tract,
    ((f.skymap_id % 10000) / 100)  AS patch_x,
    (f.skymap_id % 100)            AS patch_y
FROM pdr3_dud_rev.stellar_sequence_offsets AS f
WHERE
      (f.skymap_id / 10000) BETWEEN 10054 AND 10056
   OR (f.skymap_id / 10000) BETWEEN 9812  AND 9814
   OR (f.skymap_id / 10000) BETWEEN 9572  AND 9869;
