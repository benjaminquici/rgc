#!/usr/bin/bash

# create_source_new.py \
#     --coords_list positions_new.txt \
#     --mosaic_names G23_VIKING_K-band_regglass_reg G23_MIDAS_215MHz G23_EMU_887MHz_emu_half G23_GLASS_5500MHz_regD \
#     --mosaic_type o r r r \
#     --mean_rms 1 1 0.045 0.025 \
#     --survey_prefix viking midas emu glass \
#     --table_dir /mnt/e/Thesis/Paper_II/Pipeline/Tables/ \
#     --cat_names G23_midas_coords \
#     --cat_prefix midas \
#     --cat_ra_name ra \
#     --cat_dec_name dec \
#     --cat_flux_name int_flux \
#     --cat_err_flux_name err_int_flux \
#     --contour_colour None blue black gray \
#     --xmatch_nvss \
#     --nvss_cat G23_nvss_coords \
#     --nvss_ra RAJ2000 \
#     --nvss_dec DEJ2000 \
#     --nvss_flux int_flux \
#     --nvss_err_flux err_int_flux \
#     --lls emu glass \
#     --core glass emu \
#     --host_cat G23_viking_coords \
#     --host_ra RAcen \
#     --host_dec Deccen \
#     --host_name IAUID \
#     --host_z Z \
#     --output_dir /mnt/e/Thesis/Paper_II/Pipeline/Tables \
#     --out_cat source_cat_new \
#     --out_name Name \
#     --out_ra ra \
#     --out_dec dec \
#     --out_morph Morphology \
#     --out_frclass FR_classification \
#     --out_agnphase AGN_phase\
#     --out_rgphase RG_phase \
#     --out_emu_lls emu_LLS \
#     --out_glass_lls glass_LLS \
#     --out_emu_core_ra ra_core_887 \
#     --out_emu_core_dec dec_core_887 \
#     --out_emu_core_type core_type_887 \
#     --out_glass_core_ra ra_core_5500 \
#     --out_glass_core_dec dec_core_5500 \
#     --out_glass_core_type core_type_5500 \
#     --out_host_ra ra_host \
#     --out_host_dec dec_host \
#     --out_host_name host_iauid \
#     --out_host_z z \
#     --temp_dir /mnt/e/Thesis/Paper_II/Pipeline/tmp/

create_source_new.py \
    --coords_list positions.txt \
    --mosaic_names G23_VIKING_K-band_regglass_reg G23_MIDAS_215MHz G23_EMU_887MHz_emu_half G23_GLASS_5500MHz_regglass_reg \
    --mosaic_type o r r r \
    --mean_rms 1 1 0.045 0.025 \
    --survey_prefix viking midas emu glass \
    --table_dir /mnt/e/Thesis/Paper_II/Pipeline/Tables/ \
    --cat_names G23_midas_coords G23_emu_coords \
    --cat_prefix midas emu \
    --cat_ra_name ra ra \
    --cat_dec_name dec dec \
    --cat_flux_name int_flux int_flux \
    --cat_err_flux_name err_int_flux err_int_flux \
    --contour_colour None blue black gray \
    --xmatch_nvss \
    --nvss_cat G23_nvss_coords \
    --nvss_ra RAJ2000 \
    --nvss_dec DEJ2000 \
    --nvss_flux int_flux \
    --nvss_err_flux err_int_flux \
    --lls emu glass \
    --core glass emu \
    --host_cat G23_viking_coords \
    --host_ra RAcen \
    --host_dec Deccen \
    --host_name IAUID \
    --host_z Z \
    --output_dir /mnt/e/Thesis/Paper_II/Pipeline/Tables \
    --out_cat source_cat_new \
    --out_name Name \
    --out_ra ra \
    --out_dec dec \
    --out_morph Morphology \
    --out_frclass FR_classification \
    --out_agnphase AGN_phase\
    --out_rgphase RG_phase \
    --out_emu_lls emu_LLS \
    --out_glass_lls glass_LLS \
    --out_emu_core_ra ra_core_887 \
    --out_emu_core_dec dec_core_887 \
    --out_emu_core_type core_type_887 \
    --out_glass_core_ra ra_core_5500 \
    --out_glass_core_dec dec_core_5500 \
    --out_glass_core_type core_type_5500 \
    --out_host_ra ra_host \
    --out_host_dec dec_host \
    --out_host_name host_iauid \
    --out_host_z host_z \
    --temp_dir /mnt/e/Thesis/Paper_II/Pipeline/tmp/