module  external_files_l2c_module
implicit none   

! APC from Earth via integration over AP
! TM March 2022
character(len=220), parameter   ::   apc_file = '..\tables_L2C\SMAP_AMATRIX_V50_TM_MARCH_2022_EARTH_FOV.txt'

! correction for reflector temperature
character(len=220), parameter   ::   temp_tab_file = '..\tables_L2C\delta_refl_temp_V51C.dat'

! adjustment to galactic map based on Frank's fore-aft analysis
character(len=220), parameter   ::   gal_adj_map_file = '..\tables_L2C\smooth_gal_adjustment_maps.dat'

! further adjustment to the galactic map.  Implemented in V6.0
character(len=220), parameter   ::   gal_adj_map_file_2 = '..\tables_L2C\delta_TA_gal_refl_res_V51_IQ.dat'

! vegetation tau square table
! Aquarius ATBD
character(len=220), parameter   ::   tausq_file='..\tables_L2C\mk_vege_tausq_tab.dat'

! surface roughness emissivity table
character(len=220), parameter   ::   emiss_coeff_harm_file = '..\tables_L2C\dew_phi_VH34_harmonic_tab.dat'
character(len=200), parameter   ::   emiss_coeff_harm_file_AQ = '..\tables_L2C\deW_harm_coeffs_V9A_MI.dat'
character(len=200), parameter   ::   delta_file = '..\tables_L2C\delta_EW_V5_B.dat'		
character(len=220), parameter   ::   demiss_res_file = '..\tables_L2C\dw_res_tab_spline.dat' 

! adjusted land correction
character(len=220), parameter   ::   land_delta_tb_file = '..\tables_L2C\mk_smap_minus_model_tbmap.dat'

! fland gland gradient map
character(len=220), parameter   ::   land_gradient_file = '..\tables_L2C\mk_land_gradient_map.dat'     	

! tb meas - exp table direct access file
character(len=220), parameter   ::   dtb_statfile = '..\tables_L2C\smap_dtb_stats_ocean_V60.dat'

! wind uncertainty estimates
! Carl's file with CCMP - buoy statistics for CCMP RT
character(len=220), parameter   ::   wind_unc_file = '..\tables_L2C\CCMP_RT_minus_Buoy_W_sat_no_sat.txt'

! climatological ice mask
character(len=220), parameter   ::  ICEFILE = '..\tables_L2C\ICE4.DAT'

! sss statistics direct access file
character(len=220), parameter   ::   sss_statfile = '..\tables_L2C\smap_L2C_processing_orbital_stats_V60.dat'

! gland and fland averaged over look direction
character(len=220), parameter   ::   avg_lfrac_file = '..\tables_L2C\avg_fland_gland.dat'


end module external_files_l2c_module  

