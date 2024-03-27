! Thomas Meissner
! Remote Sensing Systems
! RSS/NASA SMAP Salinity Processing Code
! Version 6.0 Validated Release
! 03/19/2024

! V6.0
! apply shift to LR so that it matches HR
      
! Module 1: Make_SMAP_L1B_ingest_V60.F90
! Prepare ingest      
! This routine extracts the necessary variables and fields from the SMAP L1B HDF5 files (Piepmeier et al.) 
! and writes them out as binary arrays
      
! requires FORTRAN90 dll library for HDF5 and dependencies 
! RSS-Compile: Menu -> Tools -> compile32_hdf5
     
      
      
include 'smap_id_module.f90'
include 'smap_l1b_module.f90'
include 'read_smap_l1b.f90'
include 'get_filename_L1B_ingest.f90'
      
program Make_SMAP_L1B_ingest_V60

use smap_l1b_module
implicit none
      
character(len=200), parameter               ::  pathname='..\sample_data\'
character(len=200), parameter               ::  filename_in_list = 'filename_in_list.txt'
character(len=200)                          ::  ilist
character(len=200)                          ::  filename_in
character(len=200)                          ::  infile
character(len=200)                          ::  filename_out
logical(4)                                  ::  lexist      
integer(4)                                  ::  iasc,vernum,iorbit
character(len=1)                            ::  casc	  
real(4), dimension(3, maxscan)              ::  scpos,scvel,scrpy,sclla
real(4), dimension(maxcel, maxscan, nch)    ::  ta,ta_filtered
integer(4), dimension(maxcel, maxscan, nch) ::  tb_qual_flag

ilist=trim(pathname)//trim(filename_in_list)
open(unit=77,form='formatted',file=ilist,action='read',status='old')
! contains list with L1B files to be ingested

      
! initialize FORTRAN HDF
call read_smap_l1b(0,'')
write(*,*) 'HDF5 interface activated'    
      
! ascending

! initialize
scpos=-9999.
scvel=-9999.
scrpy=-9999.
sclla=-9999.
ta=-9999.
ta_filtered=-9999.

read(77,*) filename_in
read(filename_in(38:42),*) vernum   
read(filename_in(13:17),*) iorbit
read(filename_in(19:19),*) casc
if (casc=='A') then
iasc=1
else if (casc=='D') then
iasc=2
else
write(*,*) casc,' invalid casc'
stop
endif
infile = trim(pathname)//trim(filename_in)
call get_filename_l1B_ingest(iorbit, casc,    filename_out)  
      
inquire(file=infile,exist=lexist)
if (.not.(lexist)) then
write(*,*) infile,' does not exist. pgm stop.'
stop
endif
     
call read_smap_l1b(1,infile)
write(*,*) 'Opened/Read/Closed file: ',trim(infile)

scpos(1,1:numscans)=x_pos(1:numscans)
scpos(2,1:numscans)=y_pos(1:numscans)
scpos(3,1:numscans)=z_pos(1:numscans)
	  
scvel(1,1:numscans)=x_vel(1:numscans)
scvel(2,1:numscans)=y_vel(1:numscans)
scvel(3,1:numscans)=z_vel(1:numscans)
	  
scrpy(1,1:numscans)=roll(1:numscans)
scrpy(2,1:numscans)=pitch(1:numscans)
scrpy(3,1:numscans)=yaw(1:numscans)
	  
sclla(1,1:numscans)=sc_nadir_lat(1:numscans)
sclla(2,1:numscans)=sc_nadir_lon(1:numscans)
sclla(3,1:numscans)=sc_geodetic_alt_ellipsoid(1:numscans)
	  
ta(1:ncel,1:numscans,1)=ta_v(1:ncel,1:numscans)
ta(1:ncel,1:numscans,2)=ta_h(1:ncel,1:numscans)
ta(1:ncel,1:numscans,3)=ta_3(1:ncel,1:numscans)
ta(1:ncel,1:numscans,4)=ta_4(1:ncel,1:numscans)
         
ta_filtered(1:ncel,1:numscans,1)=ta_filtered_v(1:ncel,1:numscans)
ta_filtered(1:ncel,1:numscans,2)=ta_filtered_h(1:ncel,1:numscans)
ta_filtered(1:ncel,1:numscans,3)=ta_filtered_3(1:ncel,1:numscans)
ta_filtered(1:ncel,1:numscans,4)=ta_filtered_4(1:ncel,1:numscans)
         
tb_qual_flag(1:ncel,1:numscans,1)=tb_qual_flag_v(1:ncel,1:numscans)
tb_qual_flag(1:ncel,1:numscans,2)=tb_qual_flag_h(1:ncel,1:numscans)
tb_qual_flag(1:ncel,1:numscans,3)=tb_qual_flag_3(1:ncel,1:numscans)
tb_qual_flag(1:ncel,1:numscans,4)=tb_qual_flag_4(1:ncel,1:numscans)
	  
  
open(26,file=filename_out,form='unformatted',access='stream',action='write')
write(26) revNumber,iasc,ncel,numscans,vernum,                                       &
          antenna_scan_time,footprints_per_scan,scpos,scvel,scrpy,sclla,             &
          tb_time_seconds,antenna_scan_angle,tb_lat,tb_lon,ta,ta_filtered,           &
          cal_loss1_reflector,cal_temp1_reflector,tb_qual_flag
close(26)


! descending

! initialize
scpos=-9999.
scvel=-9999.
scrpy=-9999.
sclla=-9999.
ta=-9999.
ta_filtered=-9999.

read(77,*) filename_in
read(filename_in(38:42),*) vernum   
read(filename_in(13:17),*) iorbit
read(filename_in(19:19),*) casc
if (casc=='A') then
iasc=1
else if (casc=='D') then
iasc=2
else
write(*,*) casc,' invalid casc'
stop
endif
infile = trim(pathname)//trim(filename_in)
call get_filename_l1B_ingest(iorbit, casc,    filename_out)  
      
inquire(file=infile,exist=lexist)
if (.not.(lexist)) then
write(*,*) infile,' does not exist. pgm stop.'
stop
endif
     
call read_smap_l1b(1,infile)
write(*,*) 'Opened/Read/Closed file: ',trim(infile)

scpos(1,1:numscans)=x_pos(1:numscans)
scpos(2,1:numscans)=y_pos(1:numscans)
scpos(3,1:numscans)=z_pos(1:numscans)
	  
scvel(1,1:numscans)=x_vel(1:numscans)
scvel(2,1:numscans)=y_vel(1:numscans)
scvel(3,1:numscans)=z_vel(1:numscans)
	  
scrpy(1,1:numscans)=roll(1:numscans)
scrpy(2,1:numscans)=pitch(1:numscans)
scrpy(3,1:numscans)=yaw(1:numscans)
	  
sclla(1,1:numscans)=sc_nadir_lat(1:numscans)
sclla(2,1:numscans)=sc_nadir_lon(1:numscans)
sclla(3,1:numscans)=sc_geodetic_alt_ellipsoid(1:numscans)
	  
ta(1:ncel,1:numscans,1)=ta_v(1:ncel,1:numscans)
ta(1:ncel,1:numscans,2)=ta_h(1:ncel,1:numscans)
ta(1:ncel,1:numscans,3)=ta_3(1:ncel,1:numscans)
ta(1:ncel,1:numscans,4)=ta_4(1:ncel,1:numscans)
         
ta_filtered(1:ncel,1:numscans,1)=ta_filtered_v(1:ncel,1:numscans)
ta_filtered(1:ncel,1:numscans,2)=ta_filtered_h(1:ncel,1:numscans)
ta_filtered(1:ncel,1:numscans,3)=ta_filtered_3(1:ncel,1:numscans)
ta_filtered(1:ncel,1:numscans,4)=ta_filtered_4(1:ncel,1:numscans)
         
tb_qual_flag(1:ncel,1:numscans,1)=tb_qual_flag_v(1:ncel,1:numscans)
tb_qual_flag(1:ncel,1:numscans,2)=tb_qual_flag_h(1:ncel,1:numscans)
tb_qual_flag(1:ncel,1:numscans,3)=tb_qual_flag_3(1:ncel,1:numscans)
tb_qual_flag(1:ncel,1:numscans,4)=tb_qual_flag_4(1:ncel,1:numscans)
	  
  
open(26,file=filename_out,form='unformatted',access='stream',action='write')
write(26) revNumber,iasc,ncel,numscans,vernum,                                       &
          antenna_scan_time,footprints_per_scan,scpos,scvel,scrpy,sclla,             &
          tb_time_seconds,antenna_scan_angle,tb_lat,tb_lon,ta,ta_filtered,           &
          cal_loss1_reflector,cal_temp1_reflector,tb_qual_flag
close(26)

call read_smap_l1b(2,'')
write(*,*) 'HDF5 interface deactivated'    

close(77)
      
write(*,*) 'normal end to Make_SMAP_L1B_ingest_V60'
stop
end program Make_SMAP_L1B_ingest_V60