      ! updated for V5.3/V6.0
      ! shift LR obs to HR obs

      ! updated for new sim_v2 L1B file by jscott on 2015.02.23
      !
      ! written by jscott on 2015.01.29

      subroutine read_smap_l1b(icase,filename)
      use HDF5
      use smap_l1b_module
      implicit none

      character(*) :: filename
      integer(4) icase                  !icase==0 initialize HDF5 interface
                                        !icase==1 read data file
                                        !icase==2 close HDF5 interface 
	  
      character(len=200)                ::  dataset_name,group_name,attr_name
      integer(4)                        ::  file_id,group_id,dset_id,space_id,type_id,attr_id,error
      
      integer(4)                        ::  iscan, icel, ibuf
      
      integer(HSIZE_T), dimension(2)    ::  dims,maxdims
      integer(HSIZE_T), dimension(2)    ::  data_dims2
      integer(HSIZE_T), dimension(1)    ::  data_dims1
	  
      ! Initialize FORTRAN interface
      if(icase.eq.0) then
      CALL h5open_f(error)
      if(error.ne.0) stop 'error with h5open_f, pgm stopped'
      return
      endif

      !close HDF interface
      if(icase.eq.2) then
      CALL h5close_f(error)
      if(error.ne.0) stop 'error with h5close_f, pgm stopped'
      return
      endif
      
      
      
      !intialize variables to fill value
      antenna_earth_azimuth=-9999.
      antenna_look_angle=-9999.
      antenna_scan_angle=-9999.
      antenna_sidelobe_correction_3=-9999.
      antenna_sidelobe_correction_4=-9999.
      antenna_sidelobe_correction_h=-9999.
      antenna_sidelobe_correction_v=-9999.
      atm_correction_h=-9999.
      atm_correction_v=-9999.
      atm_loss_h=-9999.
      atm_loss_v=-9999.
      earth_boresight_azimuth=-9999.
      earth_boresight_incidence=-9999.
      faraday_rotation_angle=-9999.
      faraday_rotation_correction_h=-9999.
      faraday_rotation_correction_v=-9999.
      galactic_direct_correction_h=-9999.
      galactic_direct_correction_v=-9999.
      galactic_reflected_correction_3=-9999.
      galactic_reflected_correction_4=-9999.
      galactic_reflected_correction_h=-9999.
      galactic_reflected_correction_v=-9999.
      lunar_direct_phi=-9999.
      lunar_direct_theta=-9999.
      lunar_specular_correction_3=-9999.
      lunar_specular_correction_4=-9999.
      lunar_specular_correction_h=-9999.
      lunar_specular_correction_v=-9999.
      lunar_specular_lat=-9999.
      lunar_specular_lon=-9999.
      lunar_specular_phi=-9999.
      lunar_specular_reflection_coefficient_h=-9999.
      lunar_specular_reflection_coefficient_v=-9999.
      lunar_specular_theta=-9999.
      nedt_3=-9999.
      nedt_4=-9999.
      nedt_h=-9999.
      nedt_v=-9999.
      polarization_rotation_angle=-9999.
      solar_direct_correction_h=-9999.
      solar_direct_correction_v=-9999.
      solar_direct_phi=-9999.
      solar_direct_theta=-9999.
      solar_specular_correction_3=-9999.
      solar_specular_correction_4=-9999.
      solar_specular_correction_h=-9999.
      solar_specular_correction_v=-9999.
      solar_specular_lat=-9999.
      solar_specular_lon=-9999.
      solar_specular_phi=-9999.
      solar_specular_reflection_coefficient_h=-9999.
      solar_specular_reflection_coefficient_v=-9999.
      solar_specular_theta=-9999.
      specular_declination=-9999.
      specular_right_ascension=-9999.
      surface_water_fraction_mb=-9999.
      ta_3=-9999.
      ta_4=-9999.
      ta_h=-9999.
      ta_v=-9999.
      ta_filtered_3=-9999.
      ta_filtered_4=-9999.
      ta_filtered_h=-9999.
      ta_filtered_v=-9999.
      tb_3=-9999.
      tb_4=-9999.
      tb_h=-9999.
      tb_v=-9999.
      tb_declination=-9999.
      tb_upwelling=-9999.
      tb_lat=-9999.
      tb_lon=-9999.
      cal_loss12_radome=-9999.
      cal_loss1_reflector=-9999.
      cal_loss2_feed=-9999.
      cal_loss3_omt=-9999.
      cal_loss4_coupler=-9999.
      cal_loss5_diplexer=-9999.
      cal_nd_phase=-9999.
      cal_rx_phase=-9999.
      cal_temp12_radome=-9999.
      cal_temp1_reflector=-9999.
      cal_temp2_feed=-9999.
      cal_temp3_omt=-9999.
      cal_temp4_coupler=-9999.
      cal_temp5_diplexer=-9999.
      cal_temp_nd=-9999.
      cal_temp_ref=-9999.
      cal_temp_xnd=-9999.
      cal_tempref_offset=-9999.
      cal_tnd=-9999.
      cal_tref=-9999.
      cal_txnd=-9999.
      cal_xnd_phase=-9999.
      pitch=-9999.
      roll=-9999.
      sc_alongtrack_velocity=-9999.
      sc_geodetic_alt_ellipsoid=-9999.
      sc_nadir_angle=-9999.
      sc_nadir_lat=-9999.
      sc_nadir_lon=-9999.
      sc_radial_velocity=-9999.
      x_pos=-9999.
      x_vel=-9999.
      y_pos=-9999.
      y_vel=-9999.
      yaw=-9999.
      z_pos=-9999.
      z_vel=-9999.
      
      tb_mode_flag=0
      
      tb_qual_flag_3=0
      tb_qual_flag_4=0
      tb_qual_flag_h=0
      tb_qual_flag_v=0
      
      tb_right_ascension=-9999.
      
      antenna_scan_mode_flag=0
      antenna_scan_qual_flag=0
      footprints_per_scan=0
      tbs_per_scan=0
      
      tb_time_seconds=-9999.
      antenna_scan_time=-9999.
      
      
      
      ! Open an existing file
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
  
      group_name='/Brightness_Temperature'
      CALL h5gopen_f(file_id, group_name, group_id, error)
  
      dataset_name='antenna_scan_angle'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dget_space_f(dset_id, space_id, error)
      CALL h5sget_simple_extent_dims_f(space_id, dims, maxdims, error)
      CALL h5sclose_f(space_id, error)
      CALL h5dclose_f(dset_id, error)
      ncel=dims(1)
      numscans=dims(2)
	  
      dataset_name='antenna_look_angle'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, antenna_look_angle(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='antenna_scan_angle'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, antenna_scan_angle(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='ta_3'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_3(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='ta_4'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_4(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='ta_h'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_h(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='ta_v'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_v(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
      
      dataset_name='ta_filtered_3'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_filtered_3(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='ta_filtered_4'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_filtered_4(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='ta_filtered_h'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_filtered_h(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
      
      dataset_name='ta_filtered_v'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ta_filtered_v(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)

      dataset_name='tb_lat'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, tb_lat(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='tb_lon'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, tb_lon(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
      
      dataset_name='tb_mode_flag'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb_mode_flag(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)

      dataset_name='tb_qual_flag_3'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb_qual_flag_3(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='tb_qual_flag_4'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb_qual_flag_4(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='tb_qual_flag_h'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb_qual_flag_h(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='tb_qual_flag_v'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb_qual_flag_v(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='tb_time_seconds'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tb_time_seconds(1:ncel,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      CALL h5gclose_f(group_id, error)
      
          
      group_name='/Calibration_Data'
      CALL h5gopen_f(file_id, group_name, group_id, error)
      
      dataset_name='cal_loss1_reflector'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, cal_loss1_reflector(1:2,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
      
      dataset_name='cal_temp1_reflector'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, cal_temp1_reflector(1:2,1:numscans), data_dims2, error)
      CALL h5dclose_f(dset_id, error)
      
      CALL h5gclose_f(group_id, error)
	
	
	
      group_name='/Spacecraft_Data'
      CALL h5gopen_f(file_id, group_name, group_id, error)
	  
      dataset_name='antenna_scan_mode_flag'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, antenna_scan_mode_flag(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
      
      
      dataset_name='antenna_scan_time'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, antenna_scan_time(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='footprints_per_scan'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, footprints_per_scan(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='pitch'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, pitch(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='roll'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, roll(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='sc_alongtrack_velocity'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, sc_alongtrack_velocity(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='sc_geodetic_alt_ellipsoid'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, sc_geodetic_alt_ellipsoid(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='sc_nadir_angle'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, sc_nadir_angle(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='sc_nadir_lat'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, sc_nadir_lat(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='sc_nadir_lon'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, sc_nadir_lon(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='sc_radial_velocity'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, sc_radial_velocity(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='x_pos'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, x_pos(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='x_vel'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, x_vel(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='y_pos'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, y_pos(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='y_vel'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, y_vel(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='yaw'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, yaw(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='z_pos'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, z_pos(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      dataset_name='z_vel'
      CALL h5dopen_f(group_id, dataset_name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, z_vel(1:numscans), data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  
      CALL h5gclose_f(group_id, error)
	  
	  
      group_name='/Metadata/OrbitMeasuredLocation'
      CALL h5gopen_f(file_id, group_name, group_id, error)
	  
      attr_name='revNumber'
      CALL h5aopen_f(group_id, attr_name, attr_id, error)
      CALL h5aget_type_f(attr_id, type_id, error)
      CALL h5aread_f(attr_id, type_id, revNumber, data_dims1, error) !int32
      CALL h5aclose_f(attr_id, error)
	  
      CALL h5gclose_f(group_id, error)
      
	  
      !close HDF5 file
      CALL H5Fclose_f(file_id, error)
       
      
      ! added in V5.3/V6.0
      ! shift LR to match HR
      
      do iscan=1,numscans    
      do icel=1,ncel
 
      ibuf = tb_mode_flag(icel,iscan)
    
      if (abs(ta_v(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR   
      ta_v(icel,iscan) = ta_v(icel,iscan) + TA_shift_LR_HR(1)      
      endif
      endif
      
      if (abs(ta_h(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_h(icel,iscan) = ta_h(icel,iscan) + TA_shift_LR_HR(2)
      endif
      endif      
      
      if (abs(ta_3(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_3(icel,iscan) = ta_3(icel,iscan) + TA_shift_LR_HR(3)
      endif
      endif      

      if (abs(ta_4(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_4(icel,iscan) = ta_4(icel,iscan) + TA_shift_LR_HR(4)
      endif
      endif      
      
      if (abs(ta_filtered_v(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_filtered_v(icel,iscan) = ta_filtered_v(icel,iscan) + TA_shift_LR_HR(1)
      endif
      endif
      
      if (abs(ta_filtered_h(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_filtered_h(icel,iscan) = ta_filtered_h(icel,iscan) + TA_shift_LR_HR(2)
      endif
      endif      
      
      if (abs(ta_filtered_3(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_filtered_3(icel,iscan) = ta_filtered_3(icel,iscan) + TA_shift_LR_HR(3)
      endif
      endif      

      if (abs(ta_filtered_4(icel,iscan)-ta_fill)>0.1) then
      if (btest(ibuf,0)) then ! LR
      ta_filtered_4(icel,iscan) = ta_filtered_4(icel,iscan) + TA_shift_LR_HR(4)
      endif
      endif      
       

      enddo ! icel
      enddo ! iscan      
	  
      return
      end subroutine read_smap_l1b
