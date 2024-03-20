! Thomas Meissner
! Remote Sensing Systems
! RSS/NASA SMAP Salinity Processing Code
! Version 6.0 Validated Release
! 03/19/2024
      
! Module 5: Make_SMAP_L2C_V50.F90
! L2C processor  
! 1. Read L2B files.
! 2. DTB correction from ocean target calibration.
! 3. Run RSS SMAP Salinity Algorithm.
! 4. Smooth SSS to 70km.
! 5. Run perturbed retrieval for uncertainty estimates.
! 6. Calculate orbital statistic and update direct access file sss_statfile.
! 7. Create netCDF4 output.
!     

! requires FORTRAN90 dll libraries for netCDF4 and HDF5 and dependencies 
      
! RSS-Compile: Menu -> Tools -> compile64_netcdf

! internal version history
! V6.0 run
! 02 JANUARY 2024
! 1. Ingest V6 L1B
! 2. Early Mission Biases: De-bias high-rate versus low-rate
! 3. Reflected galaxy correction L2C\find_ta_gal_refl_V51.f90
! 4. DELTA_TEMP_REFL revised. use different Trefl for V/H
     !'L2C\fd_delta_temp_refl_V51C.f90'
     !'L2C\correct_reflector_emissivity_V51C.f90'
! 5. sun glint flag
     !'L2C\sun_qc_flag_V51.f90' 
     !'L2C\create_l2_qcflag_V53.f90'

! V5.3 run
! 28 MARCH 2023

! V5.0 run
! 28 JANUARY 2022
! 1. Take out SWH
! 2. Take out fice and gice from RSS AMSR2 iceflag
! 3. Add new AMSR2 discriminant iceflagging and icezones
! 4. Use NCEP 0.25 deg maps for winds, ATM, SM, TMP
! 5. The ATM correction uses NCEP 0.25 deg cloud water instead of IMERG RR.
! 6. No SSS retreival AND no land correction is performed if gland>0.1 or fland>0.1.
!    The TA expected calculation uses fland if no land correction is performed.
!    If a land correction is performed it applies the land correction and sets fland=0.
! 7. No SSS retreival AND no SC correction is performed if icezone=5 or (icezone=6 and iceflag_amsr2(2)=1).
!    The TA expected calculation uses fice if no SIC is performed.
!    If a SIC is performed it applies the SIC and sets fice=0.
! 8. Update Q/C flags fo reflect new sea-ice flagging.
! 9. Update SSS smoothing from 40-km to 70-km:
!    Use new sea-ice flagging.
!    fland t hreshold is now set to 0.005.
!10. Include formal uncertainty estimates.
!11. Use official netCDF FORTRAN library and commands.


! V 4.0 run
! 03 May 2019
! 1. gland and land correction based on 1/8 deg tables (L2A).
! 2. Add DTB in land correction accounting for the difference between model land TB and SMAP land TB.
!    L2C fd_tb_toa_lc.f90
! 3. Algorithm is run on 40km only. The result is called "sss_smap_40km". The sss field is smoothed
!    using next neighbour averaging (9-window). The result is called "sss_smap". All ancillary fields, 
!    in particular gland, gice, RR refer to the 40km cells.
! 4. Use RSS AMSR2 bytemaps (value 252) for computing fice and gice (L2B)
! 5. Use Richard's revised sun glint flag. depends on sun-glint angle and wind speed. (L2C)


! V 3.7 run
! 1. Use V3 GMF
! 2. CCMP wind speed
! 3. IMERG RR used for rain filtering
! 4. TRAN, TBUO TBDW form Liebe O2 and IMERG rain rate
! 5. I use the NOVEMBER 2015 APC as they were WITHOUT the spillover increase that I had done for V2.0
!    That gets the AMAZON right.
! 6. Ocean target calibration

include 'external_files_L2C_module.f90'
include 'l2c_module_smap_V50.f90'
include 'SMAP_ROUGHNESS_GMF_V3B_module.f90'
include 'MATRIX.f90'

include 'sss_module.f90'
include 'uncertainty_module_V50.f90'

include '..\L2B\get_filename_l2b.f90'
include 'get_filename_l2c.f90'
include '..\L2B\check_orbit.f90'

include 'allocate_L2C_arrays_V50.f90'
include 'initialize_APC_V50.f90'
include 'decode_l2b_V50.f90'

include 'fd_delta_temp_refl_V51C.f90'
include 'correct_reflector_emissivity_V51C.f90'

include 'find_dtb_bias_V53.f90'
include 'correct_cal_drift_V30.f90'

include 'find_ta_gal_refl_V51.f90'
include 'adjust_tagal_ref.f90'

include 'fd_ta_earth_V50.f90'
include 'fd_tb_toa_lc_V50.f90'
include 'fd_tb_sur_sic_V50.f90'

include 'fd_ta_expected_V50.f90'

include 'fd_sss_V50.f90'
include 'create_l2_qcflag_V53.f90'
include 'smooth_sss_V50.f90'
include 'perturbed_retrievals_V50.f90'
include 'sun_qc_flag_V51.f90'
include 'create_fill_values_V50.f90'

include 'FDICE4.f90'
include 'climatology_icemask.f90'

include 'compute_l2c_stats_V50.f90'

include 'write_l2c_netCDF_V50.f90'

include 'meissner_wentz_dielectric.f90'
include 'land_corr_step2.f90'
include 'stokes_converters.f90'

include '..\L2B\Fdsec2000.f90'    
include '..\L2A\openbig.f90'
include '..\L2A\find_month_day.f90'      
include '..\L2A\fd_date_2000.f90'	
include '..\L2A\math_routines.f90'



program MAKE_SMAP_L2C_V60 
use l2c_module_smap

implicit none
	
character(len=250)      ::  filename_l2b

integer(4)              ::  iorbit1,iorbit2, ibad
character(len=5)        ::  corbit1, corbit2

integer(4)              ::  iorbit,ires_opt,ires,ierr
	
real(8)                 ::  start_time

real(8)                 ::  secyr,secdy 
integer(4)              ::  lyear,idayjl,imon,idaymo,isecdy
real(4)                 ::  xhour

integer(4)              ::  iflag 
real(8), dimension(4)   ::  dtf_bias
real(8), dimension(2)   ::  tf_ave

logical(4)              ::  lexist


integer(4), parameter   ::  ipublish=1   !=1 for published version with time stamp in filename  

character(len=10),dimension(3)          ::  sbuf
integer(4), dimension(8)                ::  date_time


write(*,*) ' V6.0 SSS L2C processing'

! begin and end orbit are passed through command line
call get_command_argument(1, corbit1)
call get_command_argument(2, corbit2)

read(corbit1,*) iorbit1
read(corbit2,*) iorbit2           


write(*,'(a18,1x,i5.5)') ' begin orbit ',iorbit1
write(*,'(a18,1x,i5.5)') ' end orbit '  ,iorbit2

write(*,*)


! start processing
iopt=(/1,1,1,1,1,  0,0,0,0,0/)
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref

igal_wspd = 2
! =1 use provided wind speed for galaxy roughness
! =2 Frank's new galaxy model

ires_opt=1
ires=40
write(*,*) ' start SMAP FINAL SSS L2C PROCESSING V6.0'
write(*,*)

call allocate_L2C_arrays

call initialize_APC_V5


do iorbit=iorbit1,iorbit2 

    call check_orbit(iorbit, ibad)  
  	if(ibad.eq.1)     then
        write(*,8866) iorbit,' in bad orbit list. skip orbit'
        8866    format(1x,i5.5,a75)
	    cycle      
    endif
 	if(ibad.eq.2)     then
        write(*,8867) iorbit,' in bad ancillary orbit list. skip orbit'
        8867    format(1x,i5.5,a75)
	    cycle      
    endif    

    call get_filename_l2B(iorbit, filename_l2b) 
    inquire(file=filename_l2B,exist=lexist)
    if(.not.(lexist))  then
    write(*,*) filename_l2b
    write(*,*) ' L2B file does not exist. no L2C processing.'
    cycle
    endif

    write(*,*)
    write(*,*) ' processing orbit ',iorbit

    CALL DATE_AND_TIME (sbuf(1),sbuf(2),sbuf(3),date_time)
    write(*,8857) iorbit,' SSS V6.0 L2C processing start     ',&
    date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
    8857  format(1x,i5.5,a35,i4.4,'-',i2.2,'-',i2.2,' T',i2.2,':',i2.2,':',i2.2)    
      
	write(*,*) '     read in L2B file' 
	call decode_l2b(iorbit,filename_l2b, start_time,ierr)
	if (ierr/=0) then
	    write(*,*) iorbit,ierr,' bad orbit'
	    stop
	endif	

	write(*,*) '     calibration and reflector emissivity' 
	call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
	isecdy=nint(secdy)
	xhour = secdy/3600.

    ! correct reflector temperature
    call correct_reflector_emissivity_V51C
    
    ! ocean target calibration and caibration drift correction
    call find_dtb_bias(iorbit,    dtf_bias, tf_ave, iflag)
	if (iflag /=0 ) then
	    write(*,*) iorbit,' has no dtb_bias records. skip'
	    cycle
	endif
	
	dtb_bias_orbit(1:4) = dtf_bias(1:4)  !save	
    tf_ave_orbit(1:2)   = tf_ave(1:2)
  
	call correct_cal_drift_V30  
	
	call climatology_icemask ! remap ice4 climatology on SMAP swath
   
	write(*,*) '     SSS retrieval' 

   	winspd=winspd_anc ! this is the CCMP wind speed (for SMAP SSS V3/V4 files)
  
	call fd_ta_earth
	call fd_tb_toa_lc
	call fd_tb_sur_sic
	
	call fd_ta_expected
	
	call fd_sss
	
	call create_l2_qcflag
	
	call smooth_sss
	
	write(*,*)
	write(*,*) '     retrieval statistics 1 dta/dsss for orbit:' 	
	! save orbital statistics to direct access file
	call compute_l2c_stats(iorbit,start_time)
	write(*,*)	
	
	write(*,*) '     perturbed retrievals and uncertainty estimate' 
	call perturbed_retrievals

	write(*,*)
	write(*,*) '     retrieval statistics 2 dta/dsss for orbit:' 	
	! save orbital statistics to direct access file
	call compute_l2c_stats(iorbit,start_time)
	write(*,*)	
	
	
	602 continue
	write(*,*) '     write out netCDF file' 	
	! fill values
	call create_fill_values
	
	! write out netCDF file
    call write_l2c_netCDF(iorbit,ipublish, start_time)
	
    CALL DATE_AND_TIME (sbuf(1),sbuf(2),sbuf(3),date_time)
    write(*,8859) iorbit,' SSS V6.0 L2C processing complete ',&
    date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
    8859  format(1x,i5.5,a35,i4.4,'-',i2.2,'-',i2.2,' T',i2.2,':',i2.2,':',i2.2) 	
		
enddo ! iorbit	

write(*,*) 
stop ' normal end MAKE_SMAP_L2C_V60 '
end program MAKE_SMAP_L2C_V60 