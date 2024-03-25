# RSS-SMAP-Salinity-Processing-V6.0
Archive of Processing Code and Support Files for NASA/RSS SMAP Salinity Version 6.0 Validated Release

3/19/2024 - Andrew Manaster and Thomas Meissner - Remote Sensing Systems (RSS)

In this folder, users will find the routines, subroutines, and tables that have been changed in going from the NASA/RSS SMAP Salinity Version 5.0 Validated Release
to the NASA/RSS SMAP Salinity Version 6.0 Validated Release.  If a routine, subroutine, or table is not included in this folder (but is available in this folder's [parent directory](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0)), it has not changed between the V5.0 and V6.0 validated releases.  Please see the [ATBD](https://data.remss.com/smap/SSS/V06.0/documents/Release_V6.0.pdf) for full details on updates and improvements in the NASA/RSS SMAP Salinity Version 6.0 Validated Release.

The folders in this directory contain the following files.  NOTE IN ORDER TO RUN THE NASA/RSS SMAP V6 PROCESSING, THE FILES IN THIS SUBDIRECTORY, MUST REPLACE THE APPROPRIATE FILES (NOTED BELOW) IN THE FOLDERS LOCATED IN THE [PARENT DIRECTORY](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0):


## L1B_ingest
- `MAKE_SMAP_L1B_ingest_V60.F90`: Replaces `MAKE_SMAP_L1B_ingest_V50.F90` in the [V5.0 L1B_ingest.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L1B_ingest.zip) folder.
- `read_smap_l1b.f90`: Shifts the low-rate (LR) collection brightness temperatures to be consistent with the SMAP the high-rate (HR) collection (fixes early mission salinity biases near land; details found in the [SMAP V6 ATBD](https://data.remss.com/smap/SSS/V06.0/documents/Release_V6.0.pdf)).  This file replaces `read_smap_l1b.f90` in the [V5.0 L1B_ingest.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L1B_ingest.zip) folder.
- `smap_l1b_module.f90`: Contains the LR to HR shift values.  This file replaces `smap_l1b_module.f90` in the [V5.0 L1B_ingest.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L1B_ingest.zip) folder.


## L2C
- `MAKE_SMAP_L2C_V60.f90`: Replaces `MAKE_SMAP_L2C_V50.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `fd_delta_temp_refl_V51C.f90`: Corrects for the emissive SMAP reflector.  Correction updated in V6.  This replaces `fd_delta_temp_refl_V3.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `correct_reflector_emissivity_V51C.f90`: This file name was changed from V5.  This replaces `correct_reflector_emissivity_V30.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `find_ta_gal_refl_V51.f90`: File containing the correction for reflected galactic radiation.  Changed in V6 to include an additional correction to remove residual reflected galactic TA biases.  This replaces `find_ta_gal_refl.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `sun_qc_flag_V51.f90`: Code for newly updated, slightly more robust V6 sunglint flag.  This replaces `sun_qc_flag.90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `create_l2_qcflag_V53.f90`: adds the new sunglint correction to the SMAP L2C QC.  This replaces `create_l2_qcflag_V50.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `find_dtb_bias_V53.f90`: Code changed in V6 to both have a new name and use a new table.  This replaces `find_dtb_bias_V50.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `compute_l2c_stats_V53.f90`: This file name was changed from V5.  This replaces `compute_l2c_stats_V50.f90` in the [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.
- `external_files_L2C_module.f90`: Tables and paths updated for V6.  Replaces `external_files_L2C_module.f90` in [V5.0 L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/L2C.zip) folder.


## Ocean Target Calibration
- `PERFORM_OCEAN_TARGET_CALIBRATION.F90`: This file replaces `PERFORM_OCEAN_TARGET_CALIBRATION.F90` in the [V5.0 PERFORM_OCEAN_TARGET_CALIBRATION.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/blob/main/PERFORM_OCEAN_TARGET_CALIBRATION.zip) folder.  This new file uses the updated subroutines listed in the `L2C` section above and the new tables in the `Tables L2C` section below.


## Tables L2C
- `delta_refl_temp_V51C.dat`: Table to be included in `fd_delta_temp_refl_V51C.f90`.  This table replaces `delta_refl_temp_V3A.dat` in the [V5.0 validated release tables_L2C.zip](https://github.com/Remote-Sensing-Systems/RSS-SMAP-Salinity-Processing-V5.0/releases/tag/V5.0-validated-release) folder.
- `delta_TA_gal_refl_res_V51_IQ.dat`: Table used in `find_ta_gal_refl_V51.f90` as an additional correction to the reflected galaxy TA.  This is a new table used in the SMAP V6 processing.

