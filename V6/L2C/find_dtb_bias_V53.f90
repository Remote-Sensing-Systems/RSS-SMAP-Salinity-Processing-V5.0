! V5.3/V6.0
! new calibration with V5.3/V6.0 file
! the early mission biases in S3 and S4 (oribt 2812 and before) have now been taken out during L1A processing
! I have shifted S3 and S4 for the LR collection to match the HR collection
! In V5.3 I only perform a 1-time overall shift of S3 and S4

! V5.0
! same as V4.0 but redo calibration loop with V5.0 file

! V4.0
! same as V3.0 but redo calibration loop with V4.0 file
! TM 08/12/2019: changed condition for bad orbit for ibad==1 to ibad/=0 to accommodate check for bad_anc orbit

subroutine find_dtb_bias(iorbit,    dtf_bias, tf_ave, iflag)
use external_files_L2C_module 
use l2c_module_smap

implicit none

integer(4), intent(in)                      ::  iorbit
real(8), dimension(4), intent(out)          ::  dtf_bias       !1=V TF meas - exp 2=H TF meas   3=S3 meas - exp  4=S4 meas - exp
real(8), dimension(2), intent(out)          ::  tf_ave         !1=V TF exp  2=H TF exp   
integer(4), intent(out)                     ::  iflag          !0=OK  1= no bias could be computed from bias table for this orbit


integer(4), parameter			            ::  irecl=76

integer(4), parameter                       ::  maxtry=100, delay=1000, iorbit_max=150000
integer(4)                                  ::  itry, iostat, ierr	
	
integer(4), parameter   	                ::  iu=30

real(8), dimension(4)  	                    ::  xsum, ysum
integer(4), dimension(4)                    ::  nsum, msum

integer(4)						            ::  korbit, iorbit_flag

integer(4), parameter                       ::  iwinsize=3*15 ! 3-day running backward window


integer(4)              :: nobs
real(8)                 :: orbit_start_time
!
integer(4)              :: lyear
integer(4)              :: month
integer(4)              :: day
real(4)                 :: hour
real(8), dimension(2)   :: dtf, xtf
real(8), dimension(2)   :: dtf2    
    
    iflag=0 
    
	call check_orbit(iorbit,  iorbit_flag)
    if (iorbit_flag /= 0) then ! bad orbit or bad ancillary
        ! changed fro ibad == 1
        iflag=1
        dtf_bias=0.d0
        tf_ave=0.d0
    return
    endif
    
    
    ! Step 1: Open direct access file with TB meas - exp values for each orbit
     do itry = 1,maxtry
        open(unit=iu,file=dtb_statfile,form='binary',status='old',action='read',access='direct',recl=irecl,IOSTAT=iostat) 
	    if (iostat/=0) then
	        call SLEEPQQ(delay)
	    else ! file is ready for reading 
            exit
        endif        
     enddo ! try to open direct access files        

	
	! Step 2 : Calculate centered running 3-day backward average
	ysum=0.d0
	msum=0
	
    do korbit=iorbit-iwinsize,iorbit
	  
        if (korbit<1)          cycle
	    if (korbit>iorbit_max) cycle	
	    
	    call check_orbit(korbit,  iorbit_flag)
        if (iorbit_flag /= 0) cycle ! bad orbit or bad ancillary. do not include in running average

	    read(unit=iu,rec=korbit,iostat=ierr) nobs, orbit_start_time, lyear, month, day, hour, dtf(1:2), dtf2(1:2), xtf(1:2) 
		if (ierr/=0) stop ' error reading dtb correction table'
	  
	    nsum(1:4)    = nobs
    
	    xsum(1)    = dtf(1) 
	    xsum(2)    = dtf(2) 
	    xsum(3)    = xtf(1) 
	    xsum(4)    = xtf(2) 
	  
        msum = msum + nsum
        ysum = ysum + xsum*nsum

	enddo ! korbit
	
	close(iu) ! close direct access file

    if (msum(1) < 1) then
        iflag=1
        dtf_bias=0.d0
        tf_ave=0.d0
        return
    endif
 
    ysum=ysum/msum
    iflag=0
    
    dtf_bias(1) = ysum(1)
    dtf_bias(2) = ysum(2)
    tf_ave(1)   = ysum(3)
    tf_ave(2)   = ysum(4) 
 
    
    ! V5.3/V6.0
    dtf_bias(3) = offset_S3
    dtf_bias(4) = offset_S4  
   
 
return
end subroutine find_dtb_bias