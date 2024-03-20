! Thomas Meissner
! Remote Sensing Systems
! RSS/NASA SMAP Salinity Processing Code
! Version 6.0 Validated Release
! 03/19/2024
      
! Module 4: PERFORM_OCEAN_TARGET_CALIBRATION
! Perform ocean target calibration for orbit
! This routine should be run BEFORE the L2C processor.  
! It calculates TA measured - expected for open ocean scenes 
! and writes it to the binary direct access file ..\tables_l2c\smap_dtb_stats_ocean_V60.dat.
      
! RSS-Compile: Menu -> Tools -> compile64

! internal version history
! compute TF meas - exp and write to direct access file

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


! V 3.0 run
! 1. I write out Delta_TFV + Delta_TFH + STDDEV DTFV + STDDEV DTVH + TFV_exp + TFH_exp at the TA cal level
! 2. Use V3 GMF
! 3. CCMP wind speed
! 4. IMERG RR used for rain filtering
! 5. TRAN, TBUO TBDW form Liebe O2 and IMERG rain rate
! 6. I use the NOVEMBER 2015 APC as they were WITHOUT the spillover increase that I had doen for V2.0
!    That should get the Amazon about right. The ocean target calibration is done during the L2C processing.


! V2.0 run
! I write out also TA_expected
! 1. use CRID 13030 throughout.
! 2. new JPL thermal model for reflector emissivity 
! 3. emiss refl = 0.01 for V/H 
! 4. add empirical adjustment for reflector temperature
! 5. use November 2015 APC and aincreas A_II 1.108 see C:\SMAP\L2_algo_development\absolute_calibration\IDL\absolute_calibration_v20.pro
! 6. CMC SST
! 7. NCEP 0.25 deg wind speeds
! 8. Frank's new galaxy derived from SMAP for - aft look

! changed 11/17/2015:
! 1. Use WindSat/F17 wind speeds
! 2. Use updated SMAP emissivity model including GMF for S3/S4 
! 3. Include WSPD/SWH and WSPD/SST tables
! 4. Use updated APC matrix from November 2015: adjusted A_Q3 and A_43
! 5. Use running 3-day averages instead of 1-day averages in order to get better coverage 
!
! changed 09/03/2015: include correction for reflector emisivity during L2C calculation
! Frank took out the correction in the L2A processing
! it is now done at this stage.
! I use a reflector emissivity value of 0.011 and I also time interpolate the reflector temperature to get rid of the artificial
! chunkiness in the BETA version.

! updates to reflector emissivity processing: 06/01/2016
! the changes are in subroutine correct_reflector_emissivity_upd1.f90
! 1. starting with Version 3 of the SMAP L1b TA that we get from NSIDC the updated thermal model is used.
! 2. I follow my analysis and use emiss_refl = 0.0054 for both polarizations
! 3. I turn off the eclipse correction.  


include '..\L2C\external_files_L2C_module.f90'
include '..\L2C\l2c_module_smap_V50.f90'
include '..\L2C\SMAP_ROUGHNESS_GMF_V3B_module.f90'
include '..\L2C\MATRIX.f90'

include '..\L2B\get_filename_l2b.f90'
include '..\L2B\check_orbit.f90'

include '..\L2C\allocate_L2C_arrays_V50.f90'
include '..\L2C\initialize_APC_V50.f90'
include '..\L2C\decode_l2b_V50.f90'

include '..\L2C\fd_delta_temp_refl_V51C.f90'
include '..\L2C\correct_reflector_emissivity_V51C.f90'

include '..\L2C\find_ta_gal_refl_V51.f90'
include '..\L2C\adjust_tagal_ref.f90'

include '..\L2C\fd_ta_earth_V50.f90'
include '..\L2C\fd_tb_toa_lc_V50.f90'
include '..\L2C\fd_tb_sur_sic_V50.f90'

include '..\L2C\fd_ta_expected_V50.f90'

! added in V6.0
include '..\L2C\sun_qc_flag_V51.f90'

include '..\L2C\meissner_wentz_dielectric.f90'

include '..\L2C\land_corr_step2.f90'

include '..\L2C\stokes_converters.f90'

include '..\L2B\Fdsec2000.f90'    
include '..\L2A\openbig.f90'
include '..\L2A\find_month_day.f90'      
include '..\L2A\fd_date_2000.f90'	
include '..\L2A\math_routines.f90'


program PERFORM_OCEAN_TARGET_CALIBRATION
use external_files_L2C_module 
use l2c_module_smap

implicit none
	
character(len=250)      ::  filename_l2b

integer(4)              ::  iorbit1,iorbit2
character(len=5)        ::  corbit1, corbit2

integer(4)              ::  iorbit,ires_opt,ires

integer(4)              ::  ilat,ilon,idir,ierr,jerr, iflag_sun  
	
real(8)                 ::  start_time

real(8)                 ::  secyr,secdy 
integer(4)              ::  lyear,idayjl,imon,idaymo,isecdy
real(4)                 ::  xhour, xres


real(8), dimension(2)   ::  dta, dta2, xta
integer(4), dimension(2)::  num

integer(4), parameter   ::  ou=4 

integer(4), parameter   ::  maxtry=100, delay=1000, irecl=76
integer(4)              ::  itry, iostat
 
integer(4)              ::  ibad

logical(4)              ::  lexist

write(*,*) ' V6.0 SSS L2C CAL LOOP'

! begin and end orbit are passed through command line
call get_command_argument(1, corbit1)
call get_command_argument(2, corbit2)

read(corbit1,*) iorbit1
read(corbit2,*) iorbit2           


write(*,'(a18,1x,i5.5)') ' begin orbit ',iorbit1
write(*,'(a18,1x,i5.5)') ' end orbit '  ,iorbit2

write(*,*)

ires_opt=1
ires=40
xres=40.0
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

call allocate_L2C_arrays

call initialize_APC_V5
	
do iorbit=iorbit1,iorbit2 

    ! bad orbit check
    call check_orbit(iorbit, ibad)
    if (ibad ==1) then
    write(*,*) iorbit,ibad,' orbit in bad orbit list. no L2B processing.'
    cycle
    endif 
    if (ibad ==2) then
    write(*,*) iorbit,ibad,' orbit in bad ancillary list. no L2B processing.'
    cycle
    endif 
      
    num=0
    dta=0.d0
    dta2=0.d0
    xta=0.d0
    
    call get_filename_l2B(iorbit, filename_l2b) 
    inquire(file=filename_l2B,exist=lexist)
    if(.not.(lexist))  then
    write(*,*) filename_l2b
    write(*,*) ' L2B file does not exist. no processing.'
    cycle
    endif

    write(*,*)
    write(*,*) ' processing orbit ',iorbit
      
	call decode_l2b(iorbit,filename_l2b, start_time,ierr)
	if (ierr/=0) then
	    write(*,*) iorbit,ierr,' bad orbit'
	    stop
	endif
	

	call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
	isecdy=nint(secdy)
	xhour = secdy/3600.

    ! correct reflector temperature
    call correct_reflector_emissivity_V51 ! new call in V6.0
     
    
   	winspd=winspd_anc ! this is the CCMP wind speed (for SMAP SSS V3/V4 files)
    
  	call fd_ta_earth
	call fd_tb_toa_lc
	
	call fd_tb_sur_sic
	
	call fd_ta_expected
	
    do ilat=1,nlat
	do ilon=1,nlon
	do idir=1,2 
	
	if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
	if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle
	
	! I use a very conservative land/ice mask for the calbration loop
	if (gland(idir,ilon,ilat)>0.0005)  cycle
	if (fland(idir,ilon,ilat)>0.001)   cycle ! changed in V5 from 0.01 to 0.001	
	
	! changed in V5
	if (icezone(ilon,ilat) /= 0) cycle
	
	! cut out sun glint
	! use conservative sun glint from V3.0
	if (sunglt(idir,ilon,ilat)>=0.0 .and. sunglt(idir,ilon,ilat)<50. .and. alpha(idir,ilon,ilat)>30. .and. alpha(idir,ilon,ilat)<150.) cycle
	
	! add new V5.1 sun glint flag
	call compute_sun_qc(sunglt(idir,ilon,ilat), winspd(ilon,ilat), alpha(idir,ilon,ilat),   iflag_sun)
	if (iflag_sun /=0) cycle
	
	! cut out moonglint
	if (abs(monglt(idir,ilon,ilat))<15.) cycle	
	
	! cut out high galaxy
	if (ta_gal_ref(1,idir,ilon,ilat)/2.0 > 1.0) cycle	
	
	! invalid
	if (abs(ta_ant_calibrated(1,idir,ilon,ilat)-missing_val4)<0.1) cycle
	if (abs(ta_ant_calibrated(2,idir,ilon,ilat)-missing_val4)<0.1) cycle
		
	if (abs(ta_ant_exp(1,idir,ilon,ilat)-missing_val4)<0.1) cycle
	if (abs(ta_ant_exp(2,idir,ilon,ilat)-missing_val4)<0.1) cycle
	
	! rain filter
	if (rain(ilon,ilat)>0.1) cycle
 
	num(1:2)  = num(1:2)  + 1
	dta(1:2)  = dta(1:2)  + (ta_ant_calibrated(1:2,idir,ilon,ilat) - ta_ant_exp(1:2,idir,ilon,ilat))
	dta2(1:2) = dta2(1:2) + (ta_ant_calibrated(1:2,idir,ilon,ilat) - ta_ant_exp(1:2,idir,ilon,ilat))**2
	xta(1:2)  = xta(1:2)  +  ta_ant_exp(1:2,idir,ilon,ilat)

	
	enddo
    enddo
    enddo	

	if (num(1) < 2) cycle
	if (num(2) < 2) cycle

    dta=dta/num
    dta2=dta2/num
    dta2 = sqrt(dta2-dta**2)
    xta=xta/num


    write(*,*) ' dta/ta stats:'
    write(*,111) num(1), dta(1), dta2(1), xta(1)
    write(*,111) num(2), dta(2), dta2(2), xta(2)
   	111 format(1x,i12,1x,3(f10.3,1x))
   	
  	    
    ! write to direct access file
    do itry = 1,maxtry
            open(unit=ou,file=dtb_statfile,form='binary',status='old',action='write',access='direct',recl=irecl,IOSTAT=iostat) 
	        write(*,*) iostat
	        if (iostat/=0) then
	            call SLEEPQQ(delay)
	        else ! file is ready to be opened
	            write(unit=ou,rec=iorbit,iostat=jerr) num(1), start_time, lyear, imon, idaymo, xhour, dta(1:2), dta2(1:2), xta(1:2) 
	            if (jerr /=0) then
	                write(*,*) jerr
	                write(*,*) ' error opening dtb_statfile for writing orbital statistics. pgm stop.'
	                stop
	            endif 
                close(ou)
                write(*,*) ' dtb_statfile updated.'
                write(*,*) 
                exit
            endif        
    enddo ! try to open direct access files
	
    109 continue

enddo ! iorbit	

stop ' normal end PERFORM_OCEAN_TARGET_CALIBRATION'
end program PERFORM_OCEAN_TARGET_CALIBRATION