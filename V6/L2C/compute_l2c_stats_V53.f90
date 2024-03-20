! V5.3/V6.0
! new stats table
! new sun flag

subroutine compute_l2c_stats(iorbit, start_time)
use external_files_l2c_module
use l2c_module_smap
implicit none

integer(4), intent(in)  ::  iorbit
real(8), intent(in)     ::  start_time

integer(4)              ::  ilat,ilon,idir,jerr  

real(8)                 ::  secyr,secdy 
integer(4)              ::  lyear,idayjl,imon,idaymo,isecdy
real(4)                 ::  xhour

real(8), dimension(3)   ::  delta, delta2
integer(4), dimension(3)::  num
! 1=dta v
! 2=dta h
! 3=sss SMAP - HYCOM

integer(4), parameter   ::  ou=4 

integer(4), parameter   ::  maxtry=30, millisec_delay=60000, irecl=76
integer(4)              ::  itry, iostat

     num=0
    delta=0.d0
    delta2=0.d0
	
	call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
	isecdy=nint(secdy)
	xhour = secdy/3600.	

    do ilat=1,nlat
	do ilon=1,nlon
	do idir=1,2 	
		
	! invalid
	
	if (abs(ta_ant_calibrated(1,idir,ilon,ilat)-missing_val4)<0.1) cycle
	if (abs(ta_ant_calibrated(2,idir,ilon,ilat)-missing_val4)<0.1) cycle
	
	if (abs(ta_ant_exp(1,idir,ilon,ilat)-missing_val4)<0.1) cycle
	if (abs(ta_ant_exp(2,idir,ilon,ilat)-missing_val4)<0.1) cycle

	if (abs(sss_smap(idir,ilon,ilat)-missing_val4)<0.1) cycle
	if (abs(sss_ref (     ilon,ilat)-missing_val4)<0.1) cycle


	if(iqc_flag(idir,ilon,ilat)/=0) cycle
		
	
	num(1:3)    = num(1:3)      + 1
	
	delta(1:2)  = delta(1:2)    + (ta_ant_calibrated(1:2,idir,ilon,ilat) - ta_ant_exp(1:2,idir,ilon,ilat))
	delta(3)    = delta(3)      + (sss_smap(idir,ilon,ilat) - sss_ref(ilon,ilat))

	delta2(1:2)  = delta2(1:2)  + (ta_ant_calibrated(1:2,idir,ilon,ilat) - ta_ant_exp(1:2,idir,ilon,ilat))**2
	delta2(3)    = delta2(3)    + (sss_smap(idir,ilon,ilat) - sss_ref(ilon,ilat))**2
	
	
	enddo
    enddo
    enddo	

	
    where(num>=2) 
    delta=delta/num
    delta2=delta2/num
    delta2 = sqrt(delta2-delta**2)
    endwhere

    write(*,*) '          ',num(1), delta(1), delta2(1)
    write(*,*) '          ',num(2), delta(2), delta2(2)
    write(*,*) '          ',num(3), delta(3), delta2(3)
    write(*,*)

   	    
    ! write to direct access file
    do itry = 1,maxtry
            open(unit=ou,file=sss_statfile,form='binary',status='old',action='write',access='direct',recl=irecl,IOSTAT=iostat) 
	        if (iostat/=0) then
				write(*,*) ' iostat: ', iostat
				if (itry .lt. maxtry) then
					call SLEEPQQ(millisec_delay)
				else
	                write(*,*) ' error opening statfile. pgm stop.'
	                stop
				endif
	        else ! file is ready to be opened
	            write(unit=ou,rec=iorbit,iostat=jerr) num(1), start_time, lyear, imon, idaymo, xhour, delta(1:3), delta2(1:3) 
	            if (jerr /=0) then
	                write(*,*) jerr
	                write(*,*) ' error writing orbital statistics. pgm stop.'
	                stop
	            endif 
                close(ou)
                exit
            endif        
    enddo ! try to open direct access files    


return
end subroutine compute_l2c_stats