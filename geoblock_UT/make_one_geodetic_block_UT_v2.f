        program make_one_geodetic_block_UT

c       started with v10f.  small changes to work better for S003,
c       which includes 2 slow moving antennas (CD & WK)

c       V2: added starting azimuths & elevations for each station

c                  old comments below
c       Makes one geodetic bloc in UT
c       v10b allows up to 40 sources per block
c       v10c requires first source to be between azimuth
c            range -50 to +180 at all stations.  
c            This ensures we start in a unique cable warp
c            for AuScope + Ceduna antennas
c       v10d deals with Parkes 30deg elevation limit
c       v10e adds Hobart (26-m; X-Y mount) telescope options:
c            HO for generic 26-m
c            HX for temporary "X-band" only setup
c            writes out a Ceduna (CDDBBC) LST start time 
c       v10f writes out dwell time
c            Increased default "settle-time" to 30 sec to better
c            allow for accelerations/decellarations

c       started with make_geodetic_blocks_lst_10a_grid.f
c       which adapted to AuScope antennas
c            changed PT_LST to Ceduna_LST
c            changed calibrator file to read a hybrid file made by 
c                extract_flux_cat_add_position.f which combines the
c                GSFC's flux.cat with Dave Gordon's IERS positions
c            added AuScope slew and (adjustable) settle times, and
c                if only AuScope antennas, use short (30 sec) scans
c            reject a source if close on the sky to previously
c                selected source, in order to get better azimuth coverage


	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	character*32     geofile, qualfile, srcfile

c       some parameters...
	x_trials = 3000.d0
	timespan =   20.d0                               ! min
        dwell_time =  1.0d0                          ! min
	elev_min =   7.d0                               ! deg

c       Schedule in LST for a central station
c       For AuScope antennas will use Ceduna's east longitude
c        which is the middle of the array
	east_longitude = 8.92d0                          ! hours

	call set_constants 

	write(6,9000)
 9000	format(/'Enter Date and UT start time for geoblock:',
     +         /'yyyy mm dd hh mm')
	read (5,*) iyear, month, iday, ihr, min

	call get_gst (iyear, month, iday, 
     +                date, gmst_0hrUT)

	write (6,9050) date, gmst_0hrUT
 9050	format(/' Calculated JD and GMST @0hr UT:',f15.3,f12.5)

	UT = ihr + min/60.d0                    ! hours
	gmst = gmst_0hrUT + UT*sidereal_rate

	cd_lst_hr = gmst + east_longitude
	if ( cd_lst_hr .gt. 24.d0 ) cd_lst_hr = cd_lst_hr - 24.d0
	if ( cd_lst_hr .gt. 24.d0 ) cd_lst_hr = cd_lst_hr - 24.d0

	write (6,9100) cd_lst_hr
 9100	format(/' Making a geoblock for Ceduna LST =',f7.3,' hours')

	write (geofile,1000) ihr, min
 1000	format('geoblock_',i2.2,i2.2,'.out')

	write (qualfile,2000) ihr, min
 2000	format('qualities_',i2.2,i2.2,'.out')

	write (srcfile,3000) ihr, min
 3000	format('src_catalog_',i2.2,i2.2,'.out')

	write (6,4000) geofile, qualfile, srcfile
 4000	format(/'Starting ', a20,'; ',a20,'; ',a20)

	call make_one_lst_block ( x_trials,timespan,
     +       elev_min, dwell_time,
     +       east_longitude,cd_lst_hr,geofile,qualfile,srcfile)
 

	end

c===========================================================================
	subroutine make_one_lst_block ( x_trials,timespan,
     +             elev_min, dwell_time,
     +             east_longitude,x_start,geofile,qualfile,srcfile)


c       from V8a: uses VLBA skylines.   MJR 2010 March 4

c       This program figures out a good sample of sources to use for geodetic-like
c       observing blocks, used to solve for atmospheric zenith delays.
c       The idea is to observe between 10 to 20 strong compact sources
c       whose positions are known to better than 1 mas.  The sources should
c       be observed with the maximum spread bandwidth and residual multi-band
c       delays are measured.  This program assumes pure tropospheric delays
c       (ie, it uses a tropospheric mapping function).

c       This program generates large numbers of randomly chosen geodetic-like
c       blocks.  It then makes fake data and fits the data for a zenith
c       delay at each station (and a clock parameter).  The block with the
c       best quality is designated the best and information for that block is 
c       output.  Since one obtains a formal uncertainty for each antenna in the 
c       array, I define the quality as the formal uncertainty for the antenna 
c       with the worst (largest) uncertainty.  The quality factor, or formal
c       uncertainty, is in units of cm for an assumed multi-band delay error
c       of 3 cm (0.1 nsec).   

c       The sources should have sub-mas accurate positions and be strong. 
c       Eg, use the ICRF catalog with coordinate accuracies better than 1 mas,
c       X-band structure indexes of 1 or 2, and more than 1000 observations.
c       Alternatively, one may want strong (>0.3 Jy) sources if, eg, AuScope
c       antennas are used.  The program is dimensioned for up to 500 sources.

c       This version works in LST.  There are two input files.  

c         1) "stations_used.inp"    Provides stations used (name/number-code)
c             A large number of VLBI antennas are coded into the program and 
c             are selected via 2-character codes.  See the sample file.  This
c             file should list all the antennas you plan to use.

c         2) "iers_calibrator.inp"  Gives calibrator coordinates in a hybrid format
c             Any number of comment lines, starting with "!", can be included 
c             at the start of this file.

c       The program outputs 3 files:

c          1) qualities.out       which gives the quality factors of the best blocks

c          2) geo_block.out       which gives the sources to observe for the best 
c                                 block.  This is in SCHED input format and can be
c                                 copied directly into the SCHED .keyin file.

c          3) src_catalog.out     which gives the position information for the 
c                                 selected ICRF sources.  This can also be directly
c                                 copied into the source-catalog portion of the 
c                                 SCHED .keyin file.

c       CAVEATS: 

c          The program makes approximate slew time calculations, but SCHED does
c          a more accurate calculation.  Also, this program uses a simpler
c          algorithm than SCHED used to decide if a source "up" or "down" at a 
c          station.  This can lead to different slew times.  Expect differences 
c          in the length of time needed for a block.  Look at the SCHED output
c          carefully.

c       Please report any problems with the program to Mark Reid (reid@cfa.harvard.edu).

c       =============================================================================
C	Version tested with gfortran compiler; Mark J. Reid  July, 2016

	implicit real*8 ( a-h, o-z )

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant

	real*8       elev_saved(40,20), az_saved(40,20), unc_saved(40)
        real*8       sazimuth_initial(40), selev_initial(40) 
        real*8       sazimuth_old(40), selev_old(40) 
	integer      n_cal_saved(40)
	character*8  source_saved(40)
	character*2  stat_codes(40)
	character*1  sign

	character*32 geofile, qualfile, srcfile

        debug    = 0.d0

	n_ref_ant = 2               ! reference antenna #, may want to play with this

	max_num_stns  = 20          ! Dim limit: max number of antennas 
	max_num_scans = 40          ! Dim limit: max number scans  

c	write (6,7770) 
c 7770	format(' Enter control parameters (free format):',/
c     +   ' #_trials timespan CD_LSTstart  min_elev',/
c     +   '            (min)  (dec. hrs)   (deg)    eg,',/
c     +   '   1000.     45.     14.67        8.')
c	read *,x_trials, timespan, x_start, elev_min

	max_trials = int( x_trials + 0.01 )

	if ( timespan .le. 14.99d0 ) then
	   stop ' STOP: Need at least 15 minute time span.'
	endif

	CD_LST_start_hrs = x_start
	if ( CD_LST_start_hrs .lt.  0.d0 .or. 
     +       CD_LST_start_hrs .gt. 24.d0      ) then
	   print *,' CD_LST (hrs) must be between 0 - 24. STOP'
	   stop
	endif

c       Following from UT version...don't need:
c       Now want GST @ 0HR UT (very approximate)
c       call get_gst (iyr, imon, iday,  date_JD, gst0_hr ) 
c 	gast0 = gst0_hr * 15.d0 * deg_rad        ! put in common/geometry/

c        call radians_to_hhmmss ( gast0, gast_hhmmss )
c	print *,' Date (yyyy mm dd): ',iyr,imon,iday
c	print *,' Using GST @0 UT (hhmmss) =',gast_hhmmss

c       ------------------------------------------------------------------------
	seed = 0.76549d0        ! random number generator seed value
	best_quality = 9.d9
	n_trial = 0

        open (unit=10, file=qualfile, status='unknown')
	open (unit=11, file=geofile,  status='unknown')
	open (unit=12, file=srcfile,  status='unknown')

c       Add header lines to output files
	call input_stations_w_starting_azel ( stat_codes,
     +                   sazimuth_initial, selev_initial )

c       Put azimuths between 0 and 360deg
        do n = 1, max_num_stns
           if ( sazimuth_initial(n) .lt. 0.d0 )
     +          sazimuth_initial(n) = sazimuth_initial(n) + 360.d0
        enddo
        write(6,1239) ( stat_codes(n),
     +                  sazimuth_initial(n), selev_initial(n),
     +                  n=1,n_stats )
 1239   format(' Initial azimuths and elevations:',
     +         5(/1x,a2,2f6.0) )

        write (10,9000) ( stat_codes(n), n=1,n_stats )
 9000	format ('  Trial  Quality   ',20(3x,a2))

c       Do many trial geodetic blocks and save the best...
	do n_t = 1, max_trials

	   n_trial = n_trial + 1

c          Transfer initial az/el values to "old" values for each trial
           do i_stn = 1, max_num_stns
              sazimuth_old(i_stn) = sazimuth_initial(i_stn)
              selev_old(i_stn)    = selev_initial(i_stn)
           enddo

	   call fit_geodetic ( debug, seed, 
     +              east_longitude, CD_LST_start_hrs, 
     +              elev_min, dwell_time, 
     +              n_trial, max_num_scans, timespan, 
     +              quality, unc_saved, 
     +              n_scans, n_cal_saved, source_saved, 
     +              sazimuth_old, selev_old,
     +              az_saved, elev_saved )

c          Check if this trial block is better than all previous...
	   if ( quality .lt. best_quality ) then

c             Add line to "qualities.out" file
	      best_quality = quality
	      write (10,9005) n_trial, quality, 
     +                       (unc_saved(i_a),i_a=1,n_stats) 
 9005	      format(i5,f10.3,5x,20f5.2)
c             Print out quality line
	      write (6 ,9005) n_trial, quality, 
     +                       (unc_saved(i_a),i_a=1,n_stats) 
	      
c             Overwrite "geoblock.out" file with new best block
	      rewind (unit=11)

              ihr_CD_LST = x_start                            ! LST hr
              min_CD_LST = (x_start - ihr_CD_LST)*60.d0 + 0.5 ! LST min
              idwell_time_sec = dwell_time*60.d0 + 0.5        ! sec
c              write (11,9003) ihr_CD_LST, min_CD_LST, idwell_time_sec
c 9003         format('! The following gives the CD_LST start time ',
c     +         'for this block:',
c     +        /' LST=CDDBBC  START=',i2.2,':',i2.2,':00  DWELL=',i3)

	      write (11,9001) (stat_codes(n_s), n_s=1,n_stats),
     +                        (stat_codes(n_s), n_s=1,n_stats)
 9001	      format('!',19x,20(4x,a2))
	      do i_s = 1, n_scans
c                add 360deg to negative azimuths before outputting 
		 do i_a = 1, n_stats
		    if ( az_saved(i_s,i_a) .lt. 0.d0 ) then
		       az_saved(i_s,i_a) = az_saved(i_s,i_a) + 360.d0
		    endif
		 enddo
		 write (11,9010) source_saved(i_s),
     +                 (elev_saved(i_s,i_a), i_a=1,n_stats),
     +                 (az_saved(i_s,i_a), i_a=1,n_stats)
 9010		 format(' SOURCE=',1H',a8,1H',' /',2x,40(f5.0,1x))
	      enddo

c             Overwrite "src_catalog.out" file with sources in SCHED input format...
	      rewind (unit=12)
	      do i_s = 1, n_scans

		 n_s = n_cal_saved(i_s)

		 ra_hr = cal_ra(n_s) / deg_rad / 15.d0
		 ihr  = ra_hr
		 imin = (ra_hr - ihr*1.d0) * 60.d0
		 rsec = (ra_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0
		 isec = rsec
		 rsec = rsec - isec*1.d0

		 dec_rad = cal_dec(n_s)
		 if ( dec_rad .ge. 0.d0 ) then
		       sign = '+'
		    else
		       sign = '-'
		 endif
		 dec_deg = abs( dec_rad ) / deg_rad
		 ideg = dec_deg
		 iamin= (dec_deg - ideg*1.d0) * 60.d0
		 dsec = (dec_deg - ideg*1.0d0 - iamin/60.0d0) * 3600.d0
		 iasec= dsec
		 dsec = dsec - iasec*1.d0

		 write (12,9020) cal_src(n_s),ihr,imin,isec,rsec,
     +                           sign,ideg,iamin,iasec,dsec
 9020		 format('  SOURCE=',1H',a8,1H',  
     +                  '  RA=',i2.2,':',i2.2,':',i2.2,f7.6,
     +                  '  DEC=',a1,i2.2,':',i2.2,':',i2.2,f6.5,
     +                  '  EQUINOX=',1H','J2000',1H',' /')

	      enddo

	   endif

	enddo

	close (unit=10)
	close (unit=11)
	close (unit=12)

	return
	end
c==============================================================	
	subroutine fit_geodetic ( debug, seed, 
     +             east_longitude, CD_LST_start_hrs, 
     +             elev_min, dwell_time,
     +             n_trial, max_num_scans, timespan,
     +             quality, unc_saved,
     +             n_scans, n_cal_saved, source_saved,
     +             sazimuth_old, selev_old,
     +             az_saved, elev_saved )


	implicit real*8 ( a-h, o-z )

	real*8		new_params(40),param_sigmas(40),
     +			params(40)
        character*8     parnames(40)
	integer		paramids(40)

        real*8          sazimuth_old(40), selev_old(40)
        real*8          dummy_az(40), dummy_el(40)
        
     	real*8		data(20000), model(20000), resids(20000),
     +			res_err(20000),wtres(20000)

	real*8		partls(20000,40)

	real*8		gst(20000), amplitude(20000)

	integer		n_stat_a(20000), n_stat_b(20000),
     +                  n_source(20000),
     +                  num_per_stat_geo(40)

	logical		print_cor, print_solution, all_in

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant

	real*8       elev_saved(40,20), az_saved(40,20), unc_saved(40)
	integer      n_cal_saved(40)
	character*8  source_saved(40)
	character*2  stat_codes(40)


	luprint = 6

c       Dimension Limits:
	max_num_data	=20000		! number data equations
	num_params	=   40		! number of parameters

	itermax = 1                     ! linear problem
	all_in = .true.

c	------------------------------------------------------------
c       Read "stations_file.inp" to get antenna (X, Y, Z) values
c       which are placed in commons        
	call input_stations_w_starting_azel ( stat_codes,
     +                                dummy_az, dummy_el ) 

	if ( debug .gt. 0.d0 ) then
	   write (luprint,1100) (stat_codes(n),n=1,n_stats)
 1100	   format (//'   Source     ',20(4x,a2))
	endif
	  
c       Read "icrf_calibrator.inp"...
	call get_calibrators ( n_cals )

	call set_params (gain, nth_print, delay_err, 
     +		num_params, params, paramids, parnames ) 

C       --------------------------------------------------------------------
C                    Geodetic-like data input
	num_data = 0
	stm_mid = 0.d0
c       Generate fake delay data for many calibrators
	call generate_geodetic_data ( debug, seed, 
     +                       east_longitude, CD_LST_start_hrs,
     +                       max_num_data, num_data, n_cals, 
     +                       elev_min, dwell_time,
     +                       gst, gst_mid, stat_codes,
     +                       n_stat_a, n_stat_b, n_source,
     +		             num_per_stat_geo, data,
     +                       max_num_scans, n_scans, timespan,
     +                       n_cal_saved, source_saved, 
     +                       sazimuth_old, selev_old,
     +                       az_saved, elev_saved)

	if ( stm_mid .eq. 0.d0 ) stm_mid = gst_mid

C       -------------------------------------------------------------
C       Don't solve for antenna if too little data...                  
	i_0  = 0                          ! start of atmos zenith delays
	i_1  = i_0 + max_num_stns         ! start of mulit-band clock offsets
	do i_ns = 1, max_num_stns
	   if ( num_per_stat_geo(i_ns) .lt. 2 ) then
	      paramids(i_0 + i_ns) = 0                ! atmos zenith delay
	      paramids(i_1 + i_ns) = 0                ! multi-band clock
	      if ( i_ns .le. n_stats) all_in=.false.  ! flag bad
	   endif
	enddo

c       Count up number of solved-for parameters
	num_params_solved_for = 0
	do i_p = 1, num_params
	   if ( paramids(i_p) .gt. 0 ) then
	      num_params_solved_for  = num_params_solved_for + 1
	   endif
	enddo
c	write (6,1200) num_params_solved_for
c 1200	format(/' Solving for',i3,' parameters',/1x)

C       -------------------------------------------------------------
c       Following from UT version...don't need for LST version
c	ut_mid_rads = (stm_mid - gast0)/sidereal_rate   ! UT (radians)
c	call radians_to_hhmmss (ut_mid_rads, ut_mid_hhmmss)
c	write (6,1226) ut_mid_hhmmss
c 1226	format(/' Middle UT time for fits (hhmmss.s):',f12.1)

C	------------------------------------------------------------
C			Iterate Least-Squares Fitting

	idebug    = debug + 0.5
	print_cor = .false.
	iterpass  = 0          

	do while ( iterpass.le.itermax .and. all_in )

	   if ( iterpass .ge. 1 ) then

	      call calc_partials ( debug,
     +	           num_data,
     +             gst, stm_mid, 
     +             n_stat_a, n_stat_b, n_source,
     +		   num_params, params, paramids, 
     +		   partls )

c             Before calling l-s-fit, decide what it will print out...
	      if ( iterpass.eq.itermax ) then
c                   Don't print out last solution and correlations
	            print_cor = .false.
		    print_solution = .false.
		 else
		    if ( mod(iterpass-1,nth_print).eq.0 ) then
		         print_solution = .true.
		      else
			 print_solution = .false.
		    endif
	      endif

	      call least_squares_fit( idebug,print_solution,
     +                       print_cor,max_num_data,
     +		             num_params,num_data,
     +			     parnames,params,paramids,resids,
     +			     partls,res_err,new_params,param_sigmas )

	      call update_params ( debug, num_params,
     +		                   gain, params, new_params )

	   endif

	   call calc_residuals ( debug, params, 
     +		              num_data, 
     +		              gst, stm_mid, 
     +                        n_stat_a, n_stat_b, n_source,  
     +			      data, model, resids )

	   call weight_data ( num_data, delay_err,
     +                        res_err )

	   call calc_sigsq ( num_data, num_params_solved_for,  
     +                  std_err, resids, res_err, 
     +			sigsq, sigma_pdf )

	   iterpass = iterpass + 1

	enddo

c	==============================================================
c       Determine quality factor and pass station atm uncertainties back

	if ( all_in ) then
	      quality = 0.d0
	   else
	      quality = 999.d0          ! so "solution" won't be used
	endif

	do i_s = 1, n_stats

c          Quality = Max uncertainty
	   if (quality.lt.param_sigmas(i_s)) quality=param_sigmas(i_s)
	   unc_saved(i_s) = param_sigmas(i_s)

	enddo

	return
	end

c==============================================================
	subroutine set_constants 

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	c 		= 2.9979248d10			! speed light (cm/sec)
	pi 		= 4.d0*atan(1.d0)		
	twopi		= 2.d0*pi 
	secday		= 86400.d0			! seconds per solar day
	sidereal_rate	= 1.002737923d0			! solar to sidereal conv

	deg_rad  	= pi/180.d0			! convert deg's to rad's
	asec_rad	= deg_rad/3600.d0		! arcsec to radians

	return
	end

c==============================================================
        subroutine input_stations_w_starting_azel ( stat_codes,
     +                              sazimuth_old, selev_old ) 

C       Reads "station_file.inp" for station codes; the puts station coords
C            in common/geometry/

	implicit real*8 (a-h,o-z)

	character*32    station_file

	character*2     stat2, stat_codes(40)

        real*8  sazimuth_old(40), selev_old(40)
        
	logical         used

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant


c       ----------------------------------------------------------------------------
c                         Get station location info organized
	lu_stat = 9
	station_file = 'station_file.inp'              ! Hardwired station input file
	open (lu_stat, file=station_file, status='old')

c       Initialize station coordinates to zero, allowing test if station not used
	do i_st = 1, max_num_stns
	   x(i_st) = 0.d0
	   y(i_st) = 0.d0
	   z(i_st) = 0.d0
	enddo

	ieof    = 0
	n_stats = 0
	do while ( ieof .ge. 0 )

	   read (lu_stat, *, iostat=ieof)  stat2, i_stat_num,
     +           starting_az, starting_el
	   
	   if ( ieof .ge. 0 ) then

	      used = .false.

c             Check if station number is between 1 and max_num_stns...
	      if ( i_stat_num .lt. 1   .or.
     +             i_stat_num .gt. max_num_stns ) then
		 write (6,2000) stat2, i_stat_num
 2000		 format(/' Bad station # (1-20 only): ',a2,i3)
		 stop
	      endif

c             Store starting az/el values for antenna number "i_stat_num"
              sazimuth_old(i_stat_num) = starting_az
              selev_old(i_stat_num)    = starting_el
              
C	      Set station coordinates in standard *left-handed* geodetic coordinate system:
C             X toward Greenwich longitude; Y toward 90 West longitude; Z toward N pole
	      if ( stat2 .eq. 'FD' ) then
		 x(i_stat_num) = -1324009.0026d0	! FD
		 y(i_stat_num) =  5332182.0834d0	! (meters)
		 z(i_stat_num) =  3231962.4355d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'KP' ) then
		 x(i_stat_num) = -1995678.4969d0	! KP
		 y(i_stat_num) =  5037317.8209d0
		 z(i_stat_num) =  3357328.0825d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'LA' ) then
		 x(i_stat_num) = -1449752.2303d0	! LA
		 y(i_stat_num) =  4975298.7034d0
		 z(i_stat_num) =  3709123.8860d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'OV' ) then
		 x(i_stat_num) = -2409149.9782d0        ! OV
		 y(i_stat_num) =  4478573.3221d0
		 z(i_stat_num) =  3838617.3390d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'PT' ) then
		 x(i_stat_num) = -1640953.5776d0	! PT
		 y(i_stat_num) =  5014816.1165d0
		 z(i_stat_num) =  3575411.8292d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'BR' ) then
		 x(i_stat_num) = -2112064.8515d0	! BR
		 y(i_stat_num) =  3705356.6129d0
		 z(i_stat_num) =  4726813.7587d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'NL' ) then
		 x(i_stat_num) =  -130872.1216d0	! NL
		 y(i_stat_num) =  4762317.2264d0
		 z(i_stat_num) =  4226850.9983d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'HN' ) then
		 x(i_stat_num) = +1446375.2506d0	! HN
		 y(i_stat_num) =  4447939.7520d0
		 z(i_stat_num) =  4322306.0648d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'SC' ) then
		 x(i_stat_num) = +2607848.6047d0	! SC 
		 y(i_stat_num) =  5488069.8283d0
		 z(i_stat_num) =  1932739.4616d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'MK' ) then
		 x(i_stat_num) = -5464074.8410d0	! MK 
		 y(i_stat_num) =  2495249.3104d0
		 z(i_stat_num) =  2148296.7183d0
		 used          = .true.
	      endif

              if ( stat2 .eq. 'VL' ) then
                 x(i_stat_num) = -1601185.3039d0        ! VLA
                 y(i_stat_num) =  5041977.1810d0
                 z(i_stat_num) =  3554875.6407d0
                 used          = .true.
              endif

	      if ( stat2 .eq. 'EB' ) then
		 x(i_stat_num) =  4033947.46164d0	! Effelsberg 
		 y(i_stat_num) =  -486990.51504d0
		 z(i_stat_num) =  4900430.80011d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'GB' ) then
		 x(i_stat_num) =   882589.64360d0       ! Green Bank 
		 y(i_stat_num) =  4924872.32087d0
		 z(i_stat_num) =  3943729.36253d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'JV' ) then
		 x(i_stat_num) =  3822626.4970d0	! Jodrell Bank 
		 y(i_stat_num) =   154105.5889d0
		 z(i_stat_num) =  5086486.2618d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'MC' ) then
		 x(i_stat_num) =  4461369.9883d0	! Medicina 
		 y(i_stat_num) =  -919596.8324d0
		 z(i_stat_num) =  4449559.1894d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'NT' ) then
		 x(i_stat_num) =  4934563.1230d0	! Noto
		 y(i_stat_num) = -1321201.2698d0
		 z(i_stat_num) =  3806484.4778d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'ON' ) then
		 x(i_stat_num) =  3370968.1810d0	! Onsala 85
		 y(i_stat_num) =  -711464.9170d0
		 z(i_stat_num) =  5349664.1130d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'TR' ) then
		 x(i_stat_num) =  3638558.0000d0	! Torun
		 y(i_stat_num) = -1221967.0000d0
		 z(i_stat_num) =  5077041.0000d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'WB' ) then
		 x(i_stat_num) =  3828440.6400d0	! Westerbork
		 y(i_stat_num) =  -445226.0300d0
		 z(i_stat_num) =  5064923.0800d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'HH' ) then
		 x(i_stat_num) =  5085442.7805d0	! Hartebeestok
		 y(i_stat_num) = -2668263.4908d0
		 z(i_stat_num) = -2768697.0345d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'UR' ) then
		 x(i_stat_num) =   228310.726d0	        ! Urumuqui
		 y(i_stat_num) = -4631922.805d0
		 z(i_stat_num) =  4367063.964d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'CM' ) then
		 x(i_stat_num) =  3920354.8000d0	! Cambridge 32m
		 y(i_stat_num) =    -2545.7000d0
		 z(i_stat_num) =  5014285.0000d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'MH' ) then
		 x(i_stat_num) =  2892579.9681d0	! Metsahovi
		 y(i_stat_num) = -1311719.0699d0
		 z(i_stat_num) =  5512640.6897d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'HO' .or.
     +             stat2 .eq. 'HX'     ) then
		 x(i_stat_num) = -3950237.404d0  	! Hobart 26-m
		 y(i_stat_num) = -2522347.682d0
		 z(i_stat_num) = -4311561.833d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'HB' ) then
		 x(i_stat_num) = -3949990.687d0  	! Hobart 12-m
		 y(i_stat_num) = -2522421.191d0
		 z(i_stat_num) = -4311708.164d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'KE' ) then
		 x(i_stat_num) = -4147354.453d0 	! Katherine 12-m
		 y(i_stat_num) = -4581542.384d0
		 z(i_stat_num) = -1573303.310d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'YG' ) then
		 x(i_stat_num) = -2388896.040d0 	! Yarragadee 12-m
		 y(i_stat_num) = -5043349.973d0
		 z(i_stat_num) = -3078590.989d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'WA' ) then
		 x(i_stat_num) = -5115425.60d0  	! Warkworth 30-m
		 y(i_stat_num) =  -477880.31d0
		 z(i_stat_num) = -3767042.81d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'WW' ) then
		 x(i_stat_num) = -5115425.60d0  	! Temporarily using Warkworth 30-m
		 y(i_stat_num) =  -477880.31d0          ! coordinates for their 12-m
		 z(i_stat_num) = -3767042.81d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'CD' ) then
		 x(i_stat_num) = -3753443.168d0 	! Ceduna 30-m
		 y(i_stat_num) = -3912709.794d0
		 z(i_stat_num) = -3348067.060d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'MP' ) then
		 x(i_stat_num) = -4682768.1d0           ! Mopra 
		 y(i_stat_num) = -2802618.8d0
		 z(i_stat_num) = -3291759.2d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'AT' ) then
		 x(i_stat_num) = -4752447.5000d0        ! ATCA 01
		 y(i_stat_num) = -2790326.6000d0
		 z(i_stat_num) = -3200491.2900d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'PA' ) then
		 x(i_stat_num) = -4554232.0016d0        ! Parkes
		 y(i_stat_num) = -2816758.9592d0
		 z(i_stat_num) = -3454035.8457d0
		 used          = .true.
	      endif

c             Check if the station code matched one of the known station names...
	      if ( used ) then
		  n_stats = n_stats + 1
		  stat_codes(n_stats) = stat2
c		  write (6,2008) stat2, i_stat_num,
c     +                x(i_stat_num), y(i_stat_num), z(i_stat_num)
c 2008		  format(1x,a2,i5,3f15.4)
		else
		  write (6,2010) stat2, i_stat_num
 2010		  format(/' Sub input_stations: Station not a',
     +                   ' recognized 2-char code: ',a2,i3)
	      endif

	      if ( n_stats .gt. max_num_stns ) then
		 STOP ' STOP: Exceeded dimension for number stations.'
	      endif

	   endif

	enddo

	close (unit=lu_stat)

	return
	end

c==============================================================
	subroutine get_calibrators ( n_cals )

c       Reads "iers_calibrator.inp" for calibrator positions.  
c       File is in "IERS" format (Dave Gordon's internal format)
c       Info is put in common/geometry/
c       Fixed "-00" Declination problem (2012/08/20)

	implicit real*8 (a-h,o-z)

	character*32    calibrator_file

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant

	character*8  cal_name
	character*1  c1, blank1, asign

c       -----------------------------------------------------------------------
c       Input calibrator file name and positions for geodetic like data

	lu_stat         = 9
	calibrator_file = 'extracted_calibrators.inp' ! Hardwired source input file
	open (lu_stat, file=calibrator_file, status='old')
c	write (6, 3999) calibrator_file
c 3999	format(/' Opening "',a32,'" info...')

c       ------------------------------------------------------------------
c       Read in cal-source coordinates...

c       First, skip comment cards (starting with non-blank character)
	blank1 = ' '
	c1     = '!'
	do while ( c1 .eq. '!' )
	   read(lu_stat,3711) c1
 3711	   format(a1)
	enddo
	backspace (unit=lu_stat)

	n_cals = 0
	ieof   = 0
	do while ( ieof .ge. 0 .and. n_cals .lt. 500 )

           read (lu_stat,3712,iostat=ieof) cal_name, 
     +                       ira_hh,ira_mm,xra_ss,
     +                       asign,idec_dd,idec_mm,xdec_ss,
     +                       x_err, y_err, n_sessions, n_delays,
     +                       flux
 3712	   format(1x,a8,1x,2i3.2,f12.8,1x,a1,i2.2,i3.2,f11.7,
     +            1x,f7.3,1x,f7.3,2i8,f8.2)

	   if ( ieof .ge. 0 ) then

	      ra_hr  = ira_hh*1.d0 + ira_mm/60.d0 + xra_ss/3600.d0
	      ra_rad = ra_hr*pi/12.d0

	      i_sign  = 1
	      if ( asign .eq. '-' ) i_sign = -1

	      dec_deg = abs(idec_dd) + idec_mm/60.d0 + xdec_ss/3600.d0
	      dec_rad = i_sign*dec_deg*pi/180.d0
 		 
	      n_cals   = n_cals + 1
	      cal_src(n_cals) = cal_name
	      cal_ra(n_cals)  = ra_rad ! radians
	      cal_dec(n_cals) = dec_rad	! radians

c		 write (6,3720) n_cals, cal_name, ira_hh,ira_mm,xra_ss,
c     +                          asign,idec_dd,idec_mm,xdec_ss
c 3720		 format(' IERS Calibrator:',i3,1x,a8,
c     +                  5x,i2.2,i3.2,f10.6,3x,a1,i2.2,i3.2,f9.5 )
	   
	   endif

	enddo

c	print*,'DEBUG: ',n_cals
	
	close (unit=lu_stat)

	return
	end

c==============================================================
        subroutine set_params (gain, nth_print, delay_err, 
     +		num_params, params, paramids, parnames )

        implicit real*8 (a-h,o-z)

	real*8		params(40)
	character*8     parnames(40)
	integer		paramids(40)

	real*8          del_tau(40), tau_id(40), 
     +                  delay_offset(40), delay_id(40)

        character*32    control_file

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant


	nth_print= 10               
	num_iter = 1
	gain     = 1.d0

	delay_err = 3.d0             ! cm

	do i=1, max_num_stns	                 ! Station zenith delay offsets
	    del_tau(i) = 0.d0                    ! cm (any value OK to start)
	    tau_id(i)  = 1.d0                    ! solve for code set true
	enddo

	do i=1, max_num_stns                     ! Station m-b delay offsets
	    delay_offset(i) = 0.d0
	    delay_id(i)     = 1.d0               ! solve for code set false
	enddo

c       Can't solve for all antenna clocks simultaneously...
	delay_id(n_ref_ant) = 0.d0

c       ----------------------------------------------------------------
c                    Transfer values to params array...

C	Set zenith path length parameters 
	i_0  = 0
	do i = 1, max_num_stns
		params(i_0 + i) = del_tau(i)	       ! cm
	enddo

	i_1  = i_0 + max_num_stns
C	Set multi-band delay offsets parameters
	do i = 1, max_num_stns
		params(i_1 + i) = delay_offset(i)	! cm
	enddo

c       --------------------------------------------------------------
c                  Set parameter "IDs"...
        do i = 1, num_params
                paramids(i) = 1                      ! first set all to "solve for"
        enddo

	i_0  = 0
	i_1  = i_0 + max_num_stns

	do i = 1, max_num_stns	                     ! station parameters
	   if ( tau_id(i)       .eq. 0.0 ) paramids(i_0+i) = 0 ! don't solve for
	   if ( delay_id(i)     .eq. 0.0 ) paramids(i_1+i) = 0

	   if ( i .gt. n_stats ) then
	      paramids(i_0+i) = 0                    ! don't solve for
	      paramids(i_1+i) = 0                    ! don't solve for
	   endif
	enddo

	parnames(1) = 'A1 T(cm)'
	parnames(2) = 'A2 T(cm)'
	parnames(3) = 'A3 T(cm)'
	parnames(4) = 'A4 T(cm)'
	parnames(5) = 'A5 T(cm)'
	parnames(6) = 'A6 T(cm)'
	parnames(7) = 'A7 T(cm)'
	parnames(8) = 'A8 T(cm)'
	parnames(9) = 'A9 T(cm)'
	parnames(10)= 'A10T(cm)'
	parnames(11)= 'A11T(cm)'
	parnames(12)= 'A12T(cm)'
	parnames(13)= 'A13T(cm)'
	parnames(14)= 'A14T(cm)'
	parnames(15)= 'A15T(cm)'
	parnames(16)= 'A16T(cm)'
	parnames(17)= 'A17T(cm)'
	parnames(18)= 'A18T(cm)'
	parnames(19)= 'A19T(cm)'
	parnames(20)= 'A20T(cm)'

	parnames(21)= 'A1 C(cm)'
	parnames(22)= 'A2 C(cm)'
	parnames(23)= 'A3 C(cm)'
	parnames(24)= 'A4 C(cm)'
	parnames(25)= 'A5 C(cm)'
	parnames(26)= 'A6 C(cm)'
	parnames(27)= 'A7 C(cm)'
	parnames(28)= 'A8 C(cm)'
	parnames(29)= 'A9 C(cm)'
	parnames(30)= 'A10C(cm)'
	parnames(31)= 'A11C(cm)'
	parnames(32)= 'A12C(cm)'
	parnames(33)= 'A13C(cm)'
	parnames(34)= 'A14C(cm)'
	parnames(35)= 'A15C(cm)'
	parnames(36)= 'A16C(cm)'
	parnames(37)= 'A17C(cm)'
	parnames(38)= 'A18C(cm)'
	parnames(39)= 'A19C(cm)'
	parnames(40)= 'A20C(cm)'

        return
        end

c==============================================================
	subroutine generate_geodetic_data ( debug, seed, 
     +                  east_longitude, CD_LST_start_hrs,
     +                  max_num_data, num_data, n_cals, 
     +                  elev_min, dwell_time,
     +	                gst, gst_mid, stat_codes, 
     +                  n_stat_a, n_stat_b, n_source,
     +			num_per_stat_geo, data,
     +                  max_num_scans, n_scans, timespan,
     +                  n_cal_saved, source_saved, 
     +                  sazimuth_old, selev_old,
     +                  az_saved, elev_saved )

C	Generate a trial geodetic-like data set...

	implicit real*8 (a-h,o-z)

	real*8		gst(20000), data(20000)
	integer		n_stat_a(20000), n_stat_b(20000), 
     +                  n_source(20000)

	character*8     sname

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	integer         num_per_stat_geo(40)
	real*8          gst_geo(1000),delay_mb(1000),rate(1000)
	integer         n_stat_geo(1000), n_src_geo(1000)
	data            max_n_geo/1000/         ! Subroutine internal dim limit

	logical         used_stat

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant

	real*8       elev_saved(40,20), az_saved(40,20)
	character*8  source_saved(40)
	character*2  stat_codes(40)
	real*8       selev(40), sazimuth(40),
     +               selev_old(40), sazimuth_old(40)
	integer      num_dropped(40), n_cal_saved(40)
	logical      max_dropped_failed, good_angle_failed

	real*8       ra_saved(40), dec_saved(40)

	min_stats_per_scan = 0.75d0 * n_stats
	if ( min_stats_per_scan .lt. 3 ) min_stats_per_scan = 3

	max_dropped = 1       ! maximum number sources not observed at a station
	angle_min   = 10.d0   ! [deg] minimum angular separation among sources
c                             ! (in order to try and get a big azimuth range)   
c                             ! If too large, however, can put program into infinite loop!

c       -----------------------------------------------------------------
c	Create geodetic-like data...

	n_geo = 0 
c	UT_hrs = UT_start_hrs              ! starting UT [hrs]
c	ut_days = UT_hrs / 24.d0           ! [days]
	cd_lst_hrs = cd_lst_start_hrs      ! starting LST time at PT [hrs]

c       Zero counters for dropped scans for each station
	do i_ant = 1, n_stats
	   num_dropped(i_ant) = 0
	enddo

c       Randomize where starting in list of calibrators
	call gauss_random ( seed, u1, u2, g1, g2 )
	n_src_picked = u2 * (n_cals - 1) + 1

c       Generate number of geo scans...
	n_scans  = 0
	timeused = 0.d0                      ! min
	do while ( n_scans  .lt. max_num_scans .and.
     +             timeused .lt. (timespan-1.5d0)   )

 	   max_dropped_failed = .true.
	   good_angle_failed  = .true.
	   do while ( max_dropped_failed .or. good_angle_failed )

	      n_scans = n_scans + 1
	   
	      call pick_a_source (seed, cd_lst_hrs, n_scans, timespan, 
     +             n_src_picked, n_cals, stat_codes,
     +             min_stats_per_scan, 
     +             elev_min, dwell_time,
     +             sname, sra, sdec, 
     +             selev, sazimuth, selev_old, sazimuth_old,
     +             slew_dwell_time )

c             Check if this source is close on the sky to any other
c             source previously selected.   If so, skip this source.
	      call angle_on_sky_check ( sname, sra, sdec, angle_min,
     +                          n_scans, ra_saved, dec_saved,
     +                          good_angle_failed )

c             Make sure no antenna has exceeded "max_dropped"
	      max_dropped_failed = .false.
	      do i_ant = 1, n_stats
		 if ( num_dropped(i_ant) + 1 .gt. max_dropped .and.
     +                selev(i_ant) .lt. 0.d0 ) then
		    max_dropped_failed = .true.
		 endif
	      enddo

c             If exceeded "max_dropped" or not a "good_angle" 
c             (too close to a previous source), reset counter and try again
	      if ( max_dropped_failed .or. good_angle_failed) 
     +            n_scans = n_scans - 1

	   enddo

c          DEBUG printout...
c           if ( n_scans .eq. 1 ) then
c              write (6,1238) sname, slew_dwell_time
c 1238         format('DEBUG: src, slew+dwell time(min): ',a8,f8.1)
c           endif

	   timeused = timeused + slew_dwell_time ! min of UT
	   timeused_sidereal = timeused*sidereal_rate   ! min of sidereal time
	   slew_dwell_sidereal = slew_dwell_time*sidereal_rate

c          Update start time of next scan...
	   cd_lst_hrs = cd_lst_hrs + slew_dwell_sidereal/60.d0   ! hrs
	   gst_hrs    = cd_lst_hrs - 8.92d0             ! Greenwich Sidereal time (hrs)

c	   UT_hrs     = UT_hrs + slew_dwell_time/60.d0  ! hrs 
c	   ut_days    = UT_hrs / 24.d0                  ! days
	      
c          Update count of number of dropped scans for each station
	   do i_ant = 1, n_stats
	      if ( selev(i_ant) .lt. 0.d0 ) then
		 num_dropped(i_ant) = num_dropped(i_ant) + 1
	      endif
	   enddo

	   if ( debug .gt. 0.d0 ) then
	      write (6,1011) sname, (selev(n),n=1,n_stats)
 1011	      format(2x,a8,5x,20f6.1)
	   endif

c          Save source name and number to ultimately make SCHED file
	   source_saved(n_scans) = sname
	   n_cal_saved(n_scans)  = n_src_picked
	   do i_ant = 1, n_stats
	      az_saved(n_scans,i_ant)   = sazimuth(i_ant)
	      elev_saved(n_scans,i_ant) = selev(i_ant)
	   enddo

	   do i_ant = 1, n_stats

	      if ( selev( i_ant ) .gt. 0.d0 ) then

c                Source elevation not flagged...proceed to use data

c                Generate fake data.  Since we only care about the errors and
c                not the actual antenna zenith delays, we can make the data with zero
c                mean and appropriate rms to simulate real data (as if the correlator
c                model were perfect).  Very simple!
	         if ( i_ant .ne. n_ref_ant ) then
		       call gauss_random ( seed, u1,u2, g1,g2 )
		       del_mb_sec = g1 * 1.d-10          ! rms of 10^-10 sec=3 cm 

		    else
c                      Have reference antanna...
		       del_mb_sec = 0.d0
	         endif

	         n_geo = n_geo + 1

	         n_stat_geo(n_geo) = i_ant                             ! Station number
	         n_src_geo(n_geo)  = n_src_picked                      ! Cal source number

c	         ut_rads = ut_days * twopi
c	         gst_geo(n_geo)  = gast0 + ut_rads*sidereal_rate       ! GST (rads)
		 gst_geo(n_geo)  = gst_hrs * 15.d0 * deg_rad           ! GST (rads)

	         delay_mb(n_geo) = del_mb_sec * c                      ! MBdelay (cm)

	      endif     ! elevation good or flagged

	   enddo    !  stations
		 
	enddo    ! geo scans

c	print*,'DEBUG: n_scans, timeused=',n_scans,timeused

	n_geo_pts = n_geo

c       ----------------------------------------------------------------
c                  Make baseline data from station oriented data

	do n = 1, n_stats
	   num_per_stat_geo(n) = 0
	enddo
	gst_mid = 0.d0

	t_diff_max = (10.d0/86400.d0)*twopi              ! 10 sec in radians

	n_1 = n_geo_pts - 1
	do n_geo = 1, n_1

	   n_src = n_src_geo(n_geo)
	   gst_n = gst_geo(n_geo)

	   m_0 = n_geo + 1
	   do m_geo = m_0, n_geo_pts

	      m_src = n_src_geo(m_geo)
	      gst_m = gst_geo(m_geo)

c             Must match source and time to make baseline data 
	      if ( n_src .eq. m_src   .and. 
     +             abs(gst_n - gst_m) .lt. t_diff_max ) then

c                Put baseline "n_geo,m_geo" geodetic like data in large "data" array
c                multi-band delay...

		 num_data = num_data + 1

		 n_source(num_data) = n_src

		 i_a =  n_stat_geo(n_geo)
		 i_b =  n_stat_geo(m_geo)
		 n_stat_a(num_data) = i_a
		 n_stat_b(num_data) = i_b
		 num_per_stat_geo(i_a) = num_per_stat_geo(i_a) + 1
		 num_per_stat_geo(i_b) = num_per_stat_geo(i_b) + 1

		 gst(num_data)      = gst_n
		 gst_mid            = gst_mid + gst_n

c                Form delay data from station tables (sense "B minus A")
		 data(num_data)     = delay_mb(m_geo) - delay_mb(n_geo)   ! cm

	      endif

	   enddo        ! inner loop

	enddo       ! outer loop

	num_geos = num_data 
	if ( num_geos .ge. 1 ) gst_mid = gst_mid/num_geos

c	write (6,2010) num_geos
c 2010	format(' Generated',i5,' baseline geodetic-like delays')

	return
	end
c==============================================================
	subroutine pick_a_source (seed, cd_lst_hrs, n_scan, timespan,
     +             n_src_picked, n_cals, stat_codes,
     +             min_stats_per_scan,  
     +             elev_min, dwell_time,
     +             sname, sra, sdec, 
     +             selev, sazimuth, selev_old, sazimuth_old,
     +             slew_dwell_time )

c       This routine picks a source from the array "cal_src" in
c       a quasi-random manner.  It starts where it left off on the
c       previous call and skips either randomly down the list to
c       pick the next one.  

c       A source is required to be up at 75% or more stations, where
c       "up" means above 8 degrees elevation.  But a minimum of 3
c       stations is required (eg, for AuScope)
 	
	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant

	character*8  sname

	real*8       selev(40), sazimuth(40),
     +               selev_old(40), sazimuth_old(40)
	character*2  stat_codes(40), stat_code

	logical      found_bad_slew, bad_starting_az, AuScope

c       Set maximim allowed slew time; if exceeds this for any station
c       skip the source.  NB: do not set too low if you have a slow
c       moving antenna such as EB !
c	slew_time_max_allowed = 5.5d0                 ! minutes
c	slew_time_max_allowed = timespan / 8.0d0      ! minutes
	slew_time_max_allowed = timespan /10.0d0      ! eg, about 3 minutes max

c       Set up parameter for maximum number sources to skip randomly...
c	skip = 200.d0 * sqrt(1.d0/timespan)
c	skip =  97.d0 * sqrt(1.d0/timespan)
	skip =  50.d0 * sqrt(1.d0/timespan)
	
c       If first scan and initial ("old") azimuth/elevation values not set,
c       use 180/45 degrees...
	if ( n_scan .eq. 1 ) then
	   do i_a = 1, n_stats
              if ( sazimuth_old(i_a).eq.0.d0 .and.
     +             selev_old(i_a)   .eq.0.d0 ) then 
                 sazimuth_old(i_a) =  180.d0 ! deg
                 selev_old(i_a)    =   45.d0 ! deg
                 if ( stat_codes(i_a) .eq. 'EB' )
     +               sazimuth_old(i_a) =  257.d0  ! center of wrap
              endif
	   enddo
	endif

c	ut_rads = ut_days * twopi
c	stm     = gast0 + ut_rads * sidereal_rate      ! GST (rads)

c       Hardwired for Ceduna LST ***************************************************
        stm = ( cd_lst_hrs - 8.92d0 ) * 15.d0 * deg_rad ! GST (rads)

c       Pick next calibrator from list to test elevations and slews.
c       Will require source available at minimum number of stations.  
	n_good_stns    = 0
	found_bad_slew = .false.
	do while  (n_good_stns.lt.min_stats_per_scan .or. 
     +             found_bad_slew .or. bad_starting_az )

c          Pick a trial calibrator
	   call gauss_random ( seed, u1, u2, g1, g2 )
c	   if ( u1 .ge. 0.d0     .and. u1 .lt. 0.3333d0 ) i_add = 1
c	   if ( u1 .ge. 0.3333d0 .and. u1 .lt. 0.6667d0 ) i_add = 2
c	   if ( u1 .ge. 0.6667d0 .and. u1 .le. 1.d0     ) i_add = 3

	   i_add = int( u1 * skip ) + 1
           n_src_picked = n_src_picked + i_add
	   if ( n_src_picked .gt. n_cals ) n_src_picked = i_add

	   sname= cal_src( n_src_picked )
	   sra  = cal_ra ( n_src_picked )                  ! rad
	   sdec = cal_dec( n_src_picked )                  ! rad

c          Cycle through all stations, checking slews and elevations...
	   biggest_slew_time = 0.d0
	   n_good_stns       = 0
	   found_bad_slew    = .false.
	   bad_starting_az   = .false.

	   do i_a = 1, n_stats

	      xm = x(i_a)	! station coords [meters]
	      ym = y(i_a)
	      zm = z(i_a)

	      call calc_az_el ( xm, ym, zm, sra, sdec, stm, 
     +                          az, el )

c             Put az in proper range for oddball "EB" station
	      stat_code = stat_codes(i_a)
	      if ( stat_code.eq.'EB' .and. az.lt.34.d0 ) 
     +               az = az + 360.d0

c             If first scan, flag if azimuth in ambiguous range 
c             for AuScope+Ceduna...
c 	      if ( n_scan.eq.1 .and. 
c      +            (stat_code.eq.'HB' .or. stat_code.eq.'KE' .or.
c      +             stat_code.eq.'YG' .or. stat_code.eq.'CD') ) then
c 		 az_check = az
c c                force between -180 and +180 deg
c 		 if ( az.lt.-180.d0 ) az_check = az+360.d0
c 		 if ( az.gt.+180.d0 ) az_check = az-360.d0
c 		 if ( az_check.lt.-49.d0 .or. 
c      +                az_check.gt.179.d0      ) then
c 		    bad_starting_az = .true.
c 		 endif
c 	      endif

c             Check if slew time too long...
	      az_old = sazimuth_old(i_a)
	      az_new = az
	      el_old = selev_old(i_a)
	      el_new = el

c             Calculate slew time (may change "az_new" if wrapped)
	      call calc_slew_time ( az_old, az_new, el_old, el_new,
     +                              stat_code,
     +                              slew_time )

	      if ( slew_time .le. slew_time_max_allowed ) then

c                   Store az,el in output arrays 
c                   NB: az_new possibly updated in "calc_slew_time"
		    sazimuth( i_a ) = az_new ! degrees
		    selev( i_a )    = el_new

c                   Next, check if elevation allowed by "skylines" and
c                   input elev_min parameter...

c                   V8a Uses "elev_min" parameter, unless have station
c                   with "skylines" calculations.  Flag source as down
c                   if below skyline limits:  
c                       iflag=1 source up; iflag=0 source down 
		    iflag = 1

c                   Start checking with VERA antennas...
		    if ( stat_code .eq. 'VM' .or. 
     +                   stat_code .eq. 'VR' .or.  
     +                   stat_code .eq. 'VO' .or. 
     +                   stat_code .eq. 'VS'      ) then

		          call VERA_skylines( stat_code, az_new, el_new, 
     +                                        iflag)
		       else
c                         For other stations use SCHED's numbers
			  if ( el_new .gt. elev_min ) then
                             call SCHED_skylines( stat_code, 
     +                                            az_new, el_new,
     +                                            iflag )
			  endif
		    endif 

		    if ( iflag.eq.1         .and.
     +                   el_new.gt.elev_min .and. 
     +                   el_new.lt.85.d0           ) then

		          n_good_stns  = n_good_stns + 1
			  if ( slew_time .gt. biggest_slew_time )
     +                          biggest_slew_time = slew_time

		       else

		          selev( i_a )  = -9999. ! flag bad elev

		    endif	! elevation check

c                   Special case for Parkes 30deg elation limit...
		    if ( stat_code.eq.'PA' .and. 
     +                   el_new.lt.30.6d0 )  selev(i_a) = -9999.       

		 else

c                   If any station slew time bad, flag this trial source...
		    found_bad_slew = .true.
c		    print *,' PICKONE: bad slew at',stat_code,' ',
c     +                  az_old, az_new, el_old, el_new, slew_time

	      endif

	   enddo       ! run through stations for one trial source

	enddo       ! got enough stations with good elevations & slews

c       Transfer az,el for sources that were up to "old" storage array...
	do i_a = 1, n_stats
c          Check that source is not flagged (el<0) this station...
	   if ( selev(i_a) .gt. 0.d0 ) then
	      sazimuth_old(i_a) = sazimuth(i_a)
	      selev_old(i_a)    = selev(i_a)
	   endif
c	   print *,' PICKONE: Az,El: ',stat_codes(i_a),' ',
c     +                     sname,sazimuth(i_a),selev(i_a)
	enddo

c       Previously I set dwell_time to default to 1 min, 
c       but for strong sources (Scorr>0.3 Jy and desired delay 
c       precision of 0.1 nsec), reset it to 30 seconds...
c	dwell_time      = 0.5d0                                ! min

c       If any antenna is not an AuScope one, then flag is set false.
c	AuScope = .true.
c	do i_a = 1, n_stats
c	   if (stat_codes(i_a).ne.'KE'.and.stat_codes(i_a).ne.'YG'.and.
c     +         stat_codes(i_a).ne.'HB'.and.stat_codes(i_a).ne.'WW'    ) 
c     +         AuScope = .false.
c	enddo
c	if ( AuScope ) dwell_time = 0.5d0              ! min for AuScope

c       Pass out slew+dwell time
	slew_dwell_time = biggest_slew_time + dwell_time       ! min

	return
	end

c==============================================================
	subroutine angle_on_sky_check ( sname, sra, sdec, angle_min,
     +                          n_scans, ra_saved, dec_saved,
     +                          good_angle_failed )

	implicit real*8 (a-h,o-z)

	real*8       ra_saved(40), dec_saved(40)
	character*8  sname
	logical      good_angle_failed

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	good_angle_failed = .false.

	if ( n_scans .ge. 2 ) then

c          Check that all previous sources are not close in angle
	   n1 = n_scans - 1
	   do n_s = 1, n1

	      cos_dec = cos(sdec)

	      del_x = (sra - ra_saved(n_s))*cos_dec      ! rad
	      del_y = sdec - dec_saved(n_s)              ! rad

	      del_angle = sqrt( del_x**2 + del_y**2 )/deg_rad    ! deg

	      if ( del_angle .lt. angle_min ) then

c		 print*,n_s,n_scans,sname
c		 print*,sra,ra_saved(n_s),sdec,dec_saved(n_s)

		 good_angle_failed = .true.

	      endif

	   enddo

	endif

	if ( .not.good_angle_failed ) then
c          Transfer coordinates to saved arrays...
	   ra_saved(n_scans)  = sra
	   dec_saved(n_scans) = sdec
	endif

	return
	end

c==============================================================
	subroutine calc_az_el ( xm,ym,zm, ra,dec, stm, 
     +                          azimuth, elevation )

	implicit real*8 ( a-h,o-z )

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

c       A-station elevation
	stat_eq = sqrt( xm**2 + ym**2 )
	stat_lat = atan2( zm, stat_eq )       ! station latitude (rad)
	stat_w_long = atan2( ym, xm )          ! station west longitude (rad)

	sinlat  = sin(stat_lat)
	coslat  = cos(stat_lat)

c	local sidereal time of observation...
	sha  = stm - ra - stat_w_long        ! source hour angle (rad)
	cossha  = cos(sha) 
	sinsha  = sin(sha)

c	calculate source zenith angle...
	cosdec  = cos(dec)
	sindec  = sin(dec)
	cosza = sinlat*sindec + coslat*cosdec*cossha
	sinza = sin( acos(cosza) )

	zenith_angle = acos( cosza )/deg_rad    ! degrees
	elevation    = 90.d0 - zenith_angle     !    "

c       Calculate source azimuth angle...
	cosaz = ( sindec - cosza*sinlat ) / ( sinza*coslat )
	sinaz = -sinsha*cosdec / sinza
	az = atan2( sinaz, cosaz )/deg_rad      ! degrees

c       Put between 0 and 360 degrees
	if ( az .lt. 0.d0 ) az = az + 360.d0
	azimuth = az

	return
	end
c==============================================================
	subroutine calc_slew_time ( az_old, az_new, el_old, el_new,
     +                              stat_code,
     +                              slew_time )

	implicit real*8 ( a-h,o-z )

	character*2   stat_code

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad


c       Calculate nearest new azimuth to slew to
c       Az wrap numbers for VLBA and GBT...
	az_min = -90.d0
	az_max = 450.d0
	if ( stat_code .eq. 'EB' ) then
	   az_min =  34.d0
	   az_max = 480.d0
	endif
	if ( stat_code .eq. 'AT' ) then
	   az_min =-150.d0
	   az_max = 330.d0
	endif
	if ( stat_code .eq. 'CD' ) then
	   az_min =-143.d0
	   az_max = 299.d0
	endif
	if ( stat_code .eq. 'MP' ) then
	   az_min =-150.d0
	   az_max = 330.d0
	endif
	if ( stat_code .eq. 'PA' ) then
	   az_min =-155.d0
	   az_max = 295.d0
	endif

	if ( stat_code .eq. 'HB' .or.
     +       stat_code .eq. 'KE' .or.
     +       stat_code .eq. 'YG' .or.
     +       stat_code .eq. 'WA' .or.
     +       stat_code .eq. 'WW'      ) then
	   az_min =-270.d0
	   az_max = 270.d0
	endif

C       The following are guesses only...
	if ( stat_code .eq. 'HO' .or.
     +       stat_code .eq. 'HX'     ) then
	   az_min =-150.d0
	   az_max = 330.d0
	endif
	if ( stat_code .eq. 'HH' ) then
	   az_min =-150.d0
	   az_max = 330.d0
	endif

c       Put azimuth within limits...
        call check_az_limits ( az_min, az_max, 
     +                         az_new )
c       Save new azimuth value...
	az_new_out = az_new

c       Calculate 3 azimuth differences...
c       First, no-wrap difference
	az_dif = abs( az_new - az_old )

C       Next try lower azimuth wrap...
	az_new_low = az_new - 360.d0
        call check_az_limits ( az_min, az_max, 
     +                         az_new_low )
	if ( az_new_low .gt. az_min ) then

	   az_dif_low = abs( az_new_low - az_old )

	   if ( az_dif_low .lt. az_dif ) then
	      az_dif     = az_dif_low
	      az_new_out = az_new_low
	   endif

	endif

C       Finally try higher azimuth wrap...
	az_new_high = az_new + 360.d0
        call check_az_limits ( az_min, az_max, 
     +                         az_new_high )
	if ( az_new_high .lt. az_max ) then

	   az_dif_high = abs( az_new_high - az_old )

	   if ( az_dif_high .lt. az_dif ) then
	      az_dif     = az_dif_high
	      az_new_out = az_new_high
	   endif

	endif

c       Settle time (also allowing for accel/decel) 
	settle = 0.5d0                           ! minutes (default)
	if ( stat_code .eq. 'HB' .or.
     +       stat_code .eq. 'KE' .or.
     +       stat_code .eq. 'WA' .or.
     +       stat_code .eq. 'YG' .or.
     +       stat_code .eq. 'WW' ) settle = 5.d0/60.d0     ! min
c
	if ( stat_code .eq. 'CD' ) settle = 15.d0/60.d0     ! min

c       Replace az_new with possible new wrap value...
	az_new = az_new_out

	az_rate = 90.d0	            ! deg/min  VLBA default value 
	if ( stat_code .eq. 'EB' ) az_rate = 20.d0
	if ( stat_code .eq. 'GB' ) az_rate = 30.d0
	if ( stat_code .eq. 'VL' ) az_rate = 40.d0
	if ( stat_code .eq. 'HH' ) az_rate = 30.d0
	if ( stat_code .eq. 'AT' ) az_rate = 40.d0
	if ( stat_code .eq. 'WA' ) az_rate = 20.d0
	if ( stat_code .eq. 'CD' ) az_rate = 40.d0
	if ( stat_code .eq. 'HO' .or. 
     +       stat_code .eq. 'HX' ) az_rate = 40.d0
	if ( stat_code .eq. 'MP' ) az_rate = 38.d0
	if ( stat_code .eq. 'PA' ) az_rate = 15.d0
	if ( stat_code .eq. 'VM' .or.
     +       stat_code .eq. 'VR' .or.
     +       stat_code .eq. 'VO' .or.
     +       stat_code .eq. 'VS' ) az_rate = 120.d0
	if ( stat_code .eq. 'HB' .or.
     +       stat_code .eq. 'KE' .or.
     +       stat_code .eq. 'YG' .or.
     +       stat_code .eq. 'WW' ) az_rate = 300.d0

c       Calc slew time and add "start,acceleration,settle" time
	az_slew_time = az_dif / az_rate + settle     ! minutes
 
	el_dif = abs( el_new - el_old )

	el_rate = 40.d0		! deg/min  VLBA default value 
	if ( stat_code .eq. 'EB' ) el_rate = 15.d0
	if ( stat_code .eq. 'GB' ) el_rate = 15.d0
	if ( stat_code .eq. 'VL' ) el_rate = 20.d0
	if ( stat_code .eq. 'HH' ) el_rate = 30.d0
	if ( stat_code .eq. 'AT' ) el_rate = 20.d0
	if ( stat_code .eq. 'WA' ) el_rate = 20.d0
	if ( stat_code .eq. 'CD' ) el_rate = 30.d0
	if ( stat_code .eq. 'HO' .or.
     +       stat_code .eq. 'HX' ) el_rate = 40.d0
	if ( stat_code .eq. 'MP' ) el_rate = 19.d0
	if ( stat_code .eq. 'PA' ) el_rate = 15.d0
	if ( stat_code .eq. 'VM' .or.
     +       stat_code .eq. 'VR' .or.
     +       stat_code .eq. 'VO' .or.
     +       stat_code .eq. 'VS' ) el_rate = 120.d0
	if ( stat_code .eq. 'HB' .or.
     +       stat_code .eq. 'KE' .or.
     +       stat_code .eq. 'YG' .or.
     +       stat_code .eq. 'WW' ) el_rate = 75.d0

	el_slew_time = el_dif / el_rate + settle     ! minutes

	slew_time = max( az_slew_time, el_slew_time )   ! min

c       Recalculate for Hobart 26-m X,Y mount...
        if ( stat_code .eq. 'HO' .or.
     +       stat_code .eq. 'HX'     ) then

           call azel_to_xy ( az_old, el_old,
     +                        x_old,  y_old )
           call azel_to_xy ( az_new, el_new,
     +                        x_new,  y_new )
c          assigning x,y-rates as az,el-rates
           x_slew_time = abs(x_new-x_old)/az_rate + settle
           y_slew_time = abs(y_new-y_old)/el_rate + settle
           slew_time = max( x_slew_time, y_slew_time ) ! min

        endif

	return
	end

      subroutine check_az_limits ( az_min, az_max, 
     +                             az )

      implicit real*8 (a-h,o-z)

      if ( az .lt. az_min ) az = az + 360.d0
      if ( az .gt. az_max ) az = az - 360.d0

      return
      end

      subroutine azel_to_xy (az, el,
     +                        x, y  )

c     converts (az,el) in degrees to
c              (x,y) in degrees

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan(1.d0)
      deg_rad = pi/180.d0

      cosaz = cos(az*deg_rad)
      sinaz = sin(az*deg_rad)
      cosel = cos(el*deg_rad)
      sinel = cos(el*deg_rad)

      r0 = cosaz*cosel
      r1 = sinaz*cosel
      r2 = sinel

      s0 = r2
      s1 = r0
      s2 = r1

      x_rad = atan2(s1,s0)
c      if ( x_rad.lt.0.0 ) x_rad = x_rad + 2.d0*pi

      d = sqrt( s0**2 + s1**2 + s2**2 )
      y_rad = asin(s2/d)

      x = x_rad / deg_rad                    ! deg
      y = y_rad / deg_rad                    ! deg

      return
      end

c==============================================================cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     skyline checker for VERA
c
c     written by M Honma on 2007/4/27
c       adapted for make_geodetic...f by M. Reid
c
c     input : STAT2 (2 chars): 
c                'VM'=MIZ, 'VR'=IRK, 'VO'=OGA or 'VS'=ISG
c             AZ (in degree)
c             EL (in degree)
c     output: flag (integer) =1 observable, =0 unobservable
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine VERA_skylines (STAT2, AZ, EL, flag)
      Implicit none

c inputs and outputs
      character*2 STAT2
      real*8 AZ, EL
      integer flag
c for interpolation
      integer i_sk1, i_sk2
      real*8 AZ_sk1, AZ_sk2, EL_sk1, EL_sk2
      real*8 EL_sk
      real*8 sk_curr(72)
c general variables
      integer i,j

c skiline data
      real*8 skmz(72),skir(72),skis(72),skog(72)
c
      data (skmz(i),i=1,72)/
     * 1.3d0,1.3d0,1.3d0,1.3d0,1.4d0,1.5d0,1.6d0,1.6d0,1.7d0,1.7d0,
     * 1.8d0,2.0d0,2.3d0,2.2d0,2.0d0,2.0d0,2.2d0,2.6d0,2.5d0,2.3d0,
     * 2.3d0,2.3d0,2.4d0,2.4d0,2.3d0,1.9d0,1.7d0,1.5d0,1.2d0,1.2d0,
     * 1.5d0,2.0d0,2.5d0,2.4d0,2.0d0,1.3d0,1.5d0,1.0d0,2.0d0,1.7d0,
     * 1.5d0,1.4d0,1.5d0,1.7d0,2.0d0,2.3d0,2.3d0,2.2d0,2.0d0,2.0d0,
     * 1.9d0,2.0d0,2.3d0,2.6d0,2.8d0,3.0d0,2.9d0,2.7d0,2.6d0,2.4d0,
     * 2.3d0,1.9d0,1.7d0,1.7d0,1.7d0,1.5d0,1.2d0,1.1d0,1.1d0,1.0d0,
     * 1.0d0,1.0d0/
c
      data (skir(i),i=1,72)/
     * 4.5d0,5.5d0,6.5d0,7.8d0,8.5d0,9.0d0,10.4d0,11.6d0,12.0d0,11.7d0,
     * 11.5d0,11.5d0,10.0d0,6.0d0,2.0d0,1.3d0,2.7d0,2.7d0,2.7d0,1.5d0,
     * 1.0d0,0.5d0,0.0d0,0.0d0,0.7d0,2.1d0,3.0d0,3.7d0,4.0d0,4.8d0,
     * 6.2d0,7.3d0,7.2d0,6.7d0,6.0d0,6.2d0,6.2d0,6.1d0,4.8d0,2.9d0,
     * 2.1d0,1.3d0,1.2d0,1.4d0,1.6d0,1.6d0,1.5d0,1.4d0,1.3d0,1.3d0,
     * 1.8d0,1.8d0,1.8d0,1.7d0,1.5d0,2.2d0,2.2d0,1.5d0,0.3d0,0.0d0,
     * 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.3d0,1.0d0,2.4d0,
     * 3.0d0,3.6d0/
c
      data (skog(i),i=1,72)/
     * 9.5d0,11.d0,12.5d0,14.d0,15.d0,16.d0,17.d0,18.d0,18.1d0,18.2d0,
     * 17.5d0,17.5d0,17.d0,16.d0,15.d0,14.d0,13.2d0,13.4d0,14.1d0,16.d0,
     * 17.d0,18.d0,18.3d0,18.2d0,17.6d0,16.d0,13.d0,11.d0,10.2d0,11.2d0,
     * 13.d0,14.1d0,15.5d0,16.6d0,17.d0,17.d0,17.d0,16.8d0,16.8d0,
     * 16.2d0,
     * 15.7d0,15.d0,13.7d0,13.d0,12.d0,
     * 11.5d0,11.d0,11.d0,11.d0,11.2d0,
     * 11.5d0,11.8d0,12.d0,12.1d0,12.2d0,12.d0,11.8d0,11.7d0,11.5d0,
     * 11.d0,
     * 10.5d0,10.d0,9.3d0,8.7d0,8.3d0,8.d0,7.5d0,7.4d0,7.2d0,7.2d0,
     * 7.5d0,8.1d0/
c
      data (skis(i),i=1,72)/
     * 13.d0,14.3d0,15.5d0,16.7d0,17.1d0,
     * 17.7d0,19.d0,19.d0,19.d0,18.4d0,
     * 17.8d0,16.2d0,14.7d0,13.d0,10.4d0,9.2d0,9.d0,10.6d0,9.9d0,8.5d0,
     * 7.0d0,8.0d0,8.8d0,9.0d0,9.0d0,8.3d0,7.8d0,8.1d0,8.0d0,7.5d0,
     * 8.0d0,9.0d0,9.0d0,8.5d0,7.4d0,6.0d0,4.9d0,4.5d0,4.0d0,3.3d0,
     * 4.0d0,5.0d0,5.1d0,5.0d0,4.4d0,3.2d0,2.0d0,2.0d0,1.4d0,1.0d0,
     * 1.6d0,2.1d0,2.0d0,2.1d0,5.0d0,7.2d0,7.8d0,8.0d0,8.0d0,7.6d0,
     * 8.0d0,7.6d0,7.0d0,6.6d0,7.2d0,8.0d0,8.5d0,9.0d0,8.5d0,8.9d0,
     * 9.0d0,11.1d0/

c station name identification
      if (STAT2 .eq. 'VM') then
         do i=1, 72
            sk_curr(i) = skmz(i)
         end do
      end if

      if (STAT2 .eq. 'VR') then
         do i=1, 72
            sk_curr(i) = skir(i)
         end do
      end if

      if (STAT2 .eq. 'VO') then
         do i=1, 72
            sk_curr(i) = skog(i)
         end do
      end if

      if (STAT2 .eq. 'VS') then
         do i=1, 72
            sk_curr(i) = skis(i)
         end do
      end if

c     correction for 360-deg modulo in AZ
      if (AZ.lt.0.0d0) AZ = AZ + 360.0d0
      AZ = mod(AZ,360.0d0)

c     lookup for nearest skyline data point
      i_sk1 = int(AZ/5.0d0) + 1
      i_sk2 = i_sk1 + 1

c     correction for 355 < AZ < 360
      if (i_sk2.gt.72) then
         i_sk1 = i_sk1 - 1
         i_sk2 = i_sk2 - 1
      end if

c     getting skyline AZ, EL data
      AZ_sk1 = 5.0d0*(i_sk1-1)
      EL_sk1 = sk_curr(i_sk1)
      AZ_sk2 = 5.0d0*(i_sk2-1)
      EL_sk2 = sk_curr(i_sk2)
      
c     interpolation
      EL_sk = EL_sk1 + (EL_sk2-EL_sk1)/(AZ_sk2-AZ_sk1)*(AZ-AZ_sk1)

c     skyline check
      if (EL.gt.EL_sk) then
         flag = 1              ! source up
      else
         flag = 0              ! source down
      end if
c
      return
      end
c==============================================================
      Subroutine SCHED_skylines (STAT2, AZ, EL, flag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     skyline checker for VLBA (and later other) antennas
c
c     input : STAT2 (2 chars)
c             AZ (in degree)
c             EL (in degree)
c     output: flag (integer) =1 observable, =0 unobservable
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8(a-h,o-z)

c inputs and outputs
      character*2 STAT2
      real*8 AZ, EL
      integer flag

c for interpolation
      integer i_sk1, i_sk2
      real*8 AZ_sk1, AZ_sk2, EL_sk1, EL_sk2
      real*8 EL_sk

c general variables
      integer i,j
      logical found

c skyline data
      real*8 sky_az(90), sky_el(90)

c     --------------------------------------------------
c                       Station information
c     VLBA_PT antenna
      real*8 PT_az(25), PT_el(25)
      data n_PT /25/
      data PT_az/   0,  5, 60, 65, 70, 75, 80, 85,165,170,180,185,190,
     +             195,200,240,245,250,255,265,270,275,280,285,360 /
      data PT_el/   2,  2,  2,  3,  3,  2,  3,  2,  2,  3,  3,  4,  4,  
     +              3,  4,  4,  3,  4,  3,  3,  4,  3,  3,  2,  2  /

c     VLBA_PT antenna
      real*8 KP_az(34), KP_el(34)
      data n_KP /34/
      data KP_az/   0,  5, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 
     +            105,110,115,120,150,155,165,170,180,185,215,220,225,
     +            230,235,240,255,260,265,270,360 /
      data KP_el/   2,  2,  2,  5,  6,  7,  7,  8,  8,  9,  6,  3,  2,  
     +              2,  3,  3,  2,  2,  3,  3,  2,  2,  3,  3,  4,  4,
     +              5,  5,  4,  4,  3,  3,  2,  2 /

c     VLBA_LA antenna
      real*8 LA_az(18), LA_el(18)
      data n_LA /18/
      data LA_az/   0,  5, 75, 80, 85,130,135,145,150,250,255,300,305,
     +            315,320,340,345,360  /
      data LA_el/   2,  2,  2,  3,  2,  2,  3,  3,  2,  2,  3,  3,  4,
     +              4,  3,  3,  2,  2 /

c     VLBA_BR antenna
      real*8 BR_az(39), BR_el(39)
      data n_BR /39/
      data BR_az/   0,  5, 10, 15, 25, 30, 40, 45, 70, 75,120,125,130,
     +            135,155,160,185,190,195,220,225,235,240,245,250,255,
     +            265,270,275,300,305,310,315,330,335,340,345,350,360 /
      data BR_el/   2,  2,  3,  2,  2,  3,  3,  4,  4,  5,  5,  4,  4,
     +              3,  3,  2,  2,  3,  4,  4,  3,  3,  4,  4,  5,  6,
     +              6,  5,  6,  6,  5,  6,  5,  5,  4,  4,  3,  2,  2 /

c     VLBA_FD antenna
      real*8 FD_az(51), FD_el(51)
      data n_FD /51/
      data FD_az/   0,  5, 10, 15, 20, 30, 35, 40, 45, 50, 55, 60, 65,
     +             70, 75, 80, 85, 90, 95,100,105,110,115,150,155,160,
     +            220,225,230,240,245,250,255,260,265,270,275,280,285,
     +            290,295,300,305,310,315,325,330,335,340,345,360  /
      data FD_el/   5,  4,  5,  5,  3,  3,  2,  3,  2,  2,  3,  4,  7,
     +              5,  4,  4,  5,  6,  6,  5,  4,  3,  2,  2,  3,  2,
     +              2,  4,  2,  2,  3,  3,  4,  5,  5,  4,  4,  3,  3,
     +              2,  2,  3,  4,  5,  4,  4,  5,  6,  6,  5,  5 /

c     VLBA_SC antenna
      real*8 SC_az(50), SC_el(50)
      data n_SC /50/
      data SC_az/   0,  5, 10, 20, 25, 40, 45, 50, 55, 60, 65, 70, 75,
     +             80, 85, 95,100,105,110,115,120,125,130,135,140,145,
     +            150,155,160,165,175,180,185,190,200,205,210,215,220,
     +            230,235,240,245,250,260,265,270,275,280,360 /
      data SC_el/   2,  2,  3,  3,  2,  2,  3,  3,  4,  6,  6,  8,  9,
     +              9,  8,  8,  9, 10, 12, 14, 16, 16, 15, 13, 13, 12,
     +             11, 11, 10,  9,  9, 11, 13, 14, 14, 15, 13, 12, 10,
     +             10,  9,  8,  8,  7,  7,  6,  4,  3,  2,  2 /

c     VLBA_NL antenna
      real*8 NL_az(26), NL_el(26)
      data n_NL /26/
      data NL_az/   0,  5, 75, 80, 85,100,105,110,115,120,125,130,135,
     +            140,145,150,155,160,165,170,190,195,200,220,225,360 /
      data NL_el/   2,  2,  2,  3,  6,  6,  8,  7,  7,  6,  7,  7,  6,
     +              6,  7,  7,  6,  5,  4,  3,  3,  2,  3,  3,  2,  2 /

c     VLBA_OV antenna
      real*8 OV_az(55), OV_el(55)
      data n_OV /55/
      data OV_az/   0,  5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 65,
     +             70, 75, 80, 85, 90, 95,100,105,110,115,120,125,130,
     +            145,150,155,175,180,185,190,195,200,205,210,230,235,
     +            240,245,250,260,265,270,280,285,290,295,300,305,310,
     +            350,355,360 /
      data OV_el/   3,  5,  5,  7,  9, 10, 12, 13, 13, 14, 15, 15, 13,
     +             12, 11, 10,  9,  8,  7,  6,  6,  5,  4,  3,  3,  4,
     +              4,  3,  2,  2,  3,  3,  4,  5,  5,  6,  7,  7,  6,
     +              7,  8,  9,  9,  8,  7,  7,  5,  4,  4,  3,  3,  2,
     +              2,  3,  3 /

c     VLBA_MK antenna
      real*8 MK_az(46), MK_el(46)
      data n_MK /46/
      data MK_az/   0,  5, 10, 15, 20,120,125,130,135,140,145,150,155,
     +            160,165,170,175,185,190,195,200,205,210,215,220,255,
     +            260,270,275,280,285,290,295,300,305,310,315,320,325,
     +            330,335,340,345,350,355,360 /
      data MK_el/   5,  4,  3,  3,  2,  2,  4,  5,  5,  4,  4,  6,  8,
     +              8, 11, 12, 13, 13, 11, 11,  9,  7,  5,  3,  2,  2,
     +              3,  3,  5,  6,  8, 10, 12, 14, 12, 11,  9, 10, 11,
     +             10, 12, 14, 12,  9,  7,  5 /

c     VLBA_HN antenna
      real*8 HN_az(53), HN_el(53)
      data n_HN /53/
      data HN_az/   0,  5, 30, 35, 40, 45, 65, 70, 80, 85, 90, 95,100,
     +            105,110,115,120,125,130,135,140,145,150,155,160,165,
     +            170,190,195,200,205,210,220,225,230,235,240,245,250,
     +            255,270,275,290,295,315,320,325,330,335,345,350,355,
     +            360 /
      data HN_el/   6,  6,  6,  4,  5,  4,  4,  5,  5,  4,  5,  4,  4,
     +              5,  3,  4,  4,  5,  4,  6,  5,  7,  7,  5,  3,  5,
     +              4,  4,  2,  5,  5,  6,  6,  5,  6,  4,  5,  5,  4,
     +              5,  5,  4,  4,  5,  5,  6,  5,  5,  6,  6,  5,  5,
     +              6 /

c     EB_VLBA
      real*8 EB_az(36), EB_el(36)
      data n_EB /36/
      data EB_az/   0, 10, 20, 30, 40, 50, 60, 70, 80, 90,100,110,120,
     +            130,140,145,150,155,190,195,200,210,220,245,255,260,
     +            270,280,290,300,310,320,330,340,350,360 /
      data EB_el/   10, 11, 13, 15, 17, 17, 15, 13, 12, 11, 11, 11, 11,
     +              11, 10,  9,  8,  7,  7,  8,  9, 10, 10, 10, 11, 12,
     +              12, 13, 13, 12, 11, 11, 10,  9, 10, 10 /


c Get station info into working arrays
      npts = 0
      if (STAT2 .eq. 'PT') then
	 npts = n_PT
         do n=1, npts
            sky_az(n) = PT_az(n)
	    sky_el(n) = PT_el(n)
         end do
      end if
      if (STAT2 .eq. 'KP') then
	 npts = n_KP
         do n=1, npts
            sky_az(n) = KP_az(n)
	    sky_el(n) = KP_el(n)
         end do
      end if
      if (STAT2 .eq. 'LA') then
	 npts = n_LA
         do n=1, npts
            sky_az(n) = LA_az(n)
	    sky_el(n) = LA_el(n)
         end do
      end if
      if (STAT2 .eq. 'BR') then
	 npts = n_BR
         do n=1, npts
            sky_az(n) = BR_az(n)
	    sky_el(n) = BR_el(n)
         end do
      end if
      if (STAT2 .eq. 'FD') then
	 npts = n_FD
         do n=1, npts
            sky_az(n) = FD_az(n)
	    sky_el(n) = FD_el(n)
         end do
      end if
      if (STAT2 .eq. 'SC') then
	 npts = n_SC
         do n=1, npts
            sky_az(n) = SC_az(n)
	    sky_el(n) = SC_el(n)
         end do
      end if
      if (STAT2 .eq. 'NL') then
	 npts = n_NL
         do n=1, npts
            sky_az(n) = NL_az(n)
	    sky_el(n) = NL_el(n)
         end do
      end if
      if (STAT2 .eq. 'OV') then
	 npts = n_OV
         do n=1, npts
            sky_az(n) = OV_az(n)
	    sky_el(n) = OV_el(n)
         end do
      end if
      if (STAT2 .eq. 'MK') then
	 npts = n_MK
         do n=1, npts
            sky_az(n) = MK_az(n)
	    sky_el(n) = MK_el(n)
         end do
      end if
      if (STAT2 .eq. 'HN') then
	 npts = n_HN
         do n=1, npts
            sky_az(n) = HN_az(n)
	    sky_el(n) = HN_el(n)
         end do
      end if
      if (STAT2 .eq. 'EB') then
	 npts = n_EB
         do n=1, npts
            sky_az(n) = EB_az(n)
	    sky_el(n) = EB_el(n)
         end do
      end if

c     Make sure that station has skyline data...
      if ( npts .gt. 0 ) then

c        correction for 360-deg modulo in AZ
         if ( AZ .lt. 0.0d0 ) AZ = AZ + 360.0d0
         AZ = mod(AZ,360.0d0)

c        find az entry just below observed AZ
         found = .false.
         n = 0
         do while ( .not.found .and. n.lt.npts )
	    n = n + 1
	    if ( sky_az(n) .ge. AZ ) found = .true.
         enddo

         if ( .not.found ) print*,'DEBUG: no skyline for ',STAT2

	 if ( n .ge. 1 ) then
	       n1 = n - 1
	       n2 = n
	    else
	       n1 = 1
	       n2 = npts
	 endif

c        getting skyline data that bracket source AZ,EL
         AZ_sk1 = sky_az(n1)
         EL_sk1 = sky_el(n1)
         AZ_sk2 = sky_az(n2)
         EL_sk2 = sky_el(n2)
      
c        interpolation
	 EL_sk = EL_sk1+(EL_sk2-EL_sk1)*(AZ-AZ_sk1)/(AZ_sk2-AZ_sk1)

c        skyline check
	 if ( EL .gt. EL_sk ) then
	       flag = 1		! source up
	    else
	       flag = 0		! source down
c	       print *,'DEBUG: skyline flagged: ',stat2,az,el
         end if

      endif    ! checking that station has skyline data

      return
      end
c==============================================================
	subroutine calc_residuals ( debug, params, 
     +		              num_data, 
     +		              gst, stm_mid, 
     +                        n_stat_a, n_stat_b, n_source,  
     +			      data, model, resids )

 
C	Calculates models and resduals for entire data set.
C       This now includes phase data and geodetic-like data.

	implicit real*8 ( a-h,o-z )

	real*8	params(40)
	real*8	gst(20000),data(20000),model(20000),resids(20000)

	integer	n_stat_a(20000), n_stat_b(20000), n_source(20000)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant


	n_d = 0
	n_geos = num_data
	if ( n_geos .gt. 0 ) then

c          Go through geodetic data 
	   do n = 1, n_geos

	      n_d = n_d + 1
	      m_data = n_d

	      n_geo_src = n_source(n_d)

	      call set_geo_parameters (  debug, params, 
     +                  n_d, n_geo_src,
     +			n_stat_a, n_stat_b, gst,
     +			stm, ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b,
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b)


	      call calc_geo_model ( debug, stm, stm_mid, 
     +                  ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b, 
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +			tau_tot_cm )

c             Calculate delay model and residual.
	      model(n_d)  = tau_tot_cm                    ! cm
	      resids(n_d) = data(n_d) - model(n_d)

	   enddo

	endif

	return
	end
c==============================================================
	subroutine set_geo_parameters ( debug, params, 
     +                  n_d, n_geo_src,
     +			n_stat_a, n_stat_b, gst,
     +			stm, ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b,
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b )

c	Given data number (n_d), cal source number (n_geo_src), and "params" array, 
c       setup parameters needed for the geodetic model calculations
	
	implicit real*8 ( a-h,o-z )

	real*8		params(40)

	real*8		gst(20000)
	integer		n_stat_a(20000), n_stat_b(20000)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant


	stm	 = gst(n_d)

C       Transfer source positions to pass back to calling routine 
	ra  = cal_ra(n_geo_src)		   ! calibrator position (rad)
	dec = cal_dec(n_geo_src)

	i_stat_a = n_stat_a(n_d)
	i_stat_b = n_stat_b(n_d)

	x_a = x(i_stat_a)		   ! meters
	y_a = y(i_stat_a)
	z_a = z(i_stat_a)

	x_b = x(i_stat_b)
	y_b = y(i_stat_b)
	z_b = z(i_stat_b)

	n0  = 0                            ! start of zenith atmos delays
	n1  = n0 + max_num_stns            ! start of multi-band clock offsets

	del_tau_z_a	= params(n0 + i_stat_a)	            ! cm of delay
	del_tau_z_b	= params(n0 + i_stat_b)

	clock_a         = params(n1 + i_stat_a)             ! cm (delay/clight)
	clock_b         = params(n1 + i_stat_b)     

	return
	end
c==============================================================
	subroutine calc_geo_model ( debug, stm, stm_mid, 
     +                          ra, dec, 
     +				x_a,y_a,z_a, x_b,y_b,z_b, 
     +				del_tau_z_a, del_tau_z_b,
     +                          clock_a, clock_b,
     +				tau_tot_cm )

	implicit real*8 ( a-h, o-z )

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad


c       -----------------------------------------------------------
c	Calculate atmospheric delay error for geodetic cal data:

c	For station A...	
	del_tau_dot_a = 0.d0
	call atmosphere_shift (debug, stm, stm_mid,
     +                  ra,dec, x_a,y_a,z_a,		        
     +		        del_tau_z_a,del_tau_dot_a,
     +			del_tau_a_cm)                                   
	
c	For station B...	
	del_tau_dot_b = 0.d0
	call atmosphere_shift (debug, stm, stm_mid, 
     +                  ra,dec, x_b,y_b,z_b, 
     +			del_tau_z_b,del_tau_dot_b,
     +			del_tau_b_cm)                                  

c	Difference for baseline delay...	
	tau_atm_cm = del_tau_b_cm - del_tau_a_cm          ! cm

c       Calculate clock offset...
	tau_clock_cm= clock_b - clock_a                   ! cm

c       Now add the atmospheric and clock delays...
	tau_tot_cm = tau_atm_cm + tau_clock_cm            ! cm

	return 
	end
c==============================================================
	subroutine atmosphere_shift(debug,stm,stm_mid,ra,dec,	
     +				xm,ym,zm,del_tau_z,del_tau_dot,
     +				tau_cm)

c	Returns excess atmospheric delay in cm
c       for given station owing to a vertical path error in cm.  

c       This version assumes the excesss path is from a dry
c	atmospheric term.  (Commented out equations shows what
c       form it would be for a wet component instead.)

c       See Thompson, Moran, Swenson Eq 13.41 

	implicit real*8 (a-h, o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	real*8	x(20), y(20), z(20), 
     +          cal_ra(500), cal_dec(500)
	character*8  cal_src(500)
	common /geometry/ x,y,z, gast0, 
     +     cal_src, cal_ra, cal_dec,
     +     max_num_stns, n_stats, n_ref_ant


	rad_hr	= 12.d0/pi 

	call calc_az_el ( xm, ym, zm, ra, dec, stm, 
     +                          az, el )

	za = 90.d0 - el                                ! degrees
	cos_za = cos( za * deg_rad )
	sec_za = 1.d0 / cos_za
	tan_za = tan( za * deg_rad )

	del_tau_z_total = del_tau_z + del_tau_dot*(stm-stm_mid)*rad_hr    ! cm

c       For consistency with TMS equations, convert to mBars (and back to cm)...
	Po      = del_tau_z_total  / 0.228d0		                  ! mBars
	tau_dry = 0.228d0*Po*sec_za * (1.d0-0.0013d0*tan_za*tan_za)	  ! cm

c       Were we to think the excess delay was from a wet atmosphere, then...
c	pVo     = ??			! partial pressure H2O at surface in mBar
c					! 5 mBar ~ 0.72 cm prec H20 => approx 5 cm path delay
c	T       = 280.d0		! Surface temperature (K)         
c	tau_wet = 7.5d4*pVo/(T*T)       ! wet delay at Zenith          
c       tau_wet = 7.5d4*pVo*sec_za*(1.d0-0.0003d0*tan_za*tan_za)/(T*T)	  ! cm
c       NB: the only difference is the smaller coefficient of tan^2(za)

	tau_cm	  = tau_dry	                        ! cm

	return
	end
c==============================================================
	subroutine calc_partials ( debug,
     +	           num_data,
     +             gst, stm_mid, 
     +             n_stat_a, n_stat_b, n_source,
     +		   num_params, params, paramids, 
     +		   partls )

	implicit real*8 ( a-h,o-z )

	real*8	params(40), params_wiggled(40)
	real*8	partls(20000,40), gst(20000)

	integer	n_stat_a(20000), n_stat_b(20000), n_source(20000)
	integer	paramids(40)

c	Set secondary parameter array and zero partials array...
	do i_p = 1, num_params
		params_wiggled(i_p) = params(i_p)
		do i_d = 1, num_data
			partls(i_d,i_p) = 0.d0
		enddo
	enddo

c       ----------------------------------------------------------------
c       Now run through geodetic-like data...
	if ( num_data .gt. 0 ) then

	   do n_d = 1, num_data

	      n_geo_src = n_source(n_d)

	      call set_geo_parameters (  debug, params, 
     +                  n_d, n_geo_src,
     +			n_stat_a, n_stat_b, gst,
     +			stm, ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b,
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b )

	      call calc_geo_model ( debug, stm, stm_mid, 
     +                  ra, dec, 
     +		        x_a,y_a,z_a, x_b,y_b,z_b, 
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +			tau )

	      do i_p = 1, num_params

		 if ( paramids(i_p) .eq. 1 ) then

C                   Change parameters by 1 part in 10^4; but also never
C			by less than a value of 10^-4 (arcsec or cm)...
		    del_param = abs( params(i_p) ) * 1.d-04
		    if ( del_param .lt. 1.d-04 ) del_param = 1.d-04

		       params_wiggled(i_p) = params(i_p) + del_param
			
		       call set_geo_parameters (  debug, 
     +                          params_wiggled, 
     +                          n_d, n_geo_src,
     +                          n_stat_a, n_stat_b, gst,
     +			        stm, ra, dec, 
     +			        x_a,y_a,z_a, x_b,y_b,z_b,
     +			        del_tau_z_a, del_tau_z_b,
     +                          clock_a, clock_b )

		       call calc_geo_model ( debug, stm, stm_mid, 
     +                          ra, dec, 
     +				x_a,y_a,z_a, x_b,y_b,z_b, 
     +				del_tau_z_a, del_tau_z_b,
     +                          clock_a, clock_b,
     +				tau_wiggled )

c                       First, calculate multi-band delay partial
			partls(n_d,i_p) = (tau_wiggled-tau)/
     +                                     del_param

c                       Go back to original parameter value...
			params_wiggled(i_p)= params(i_p) 

		   endif

		enddo

	     enddo

	endif

	return
	end
c==============================================================
	subroutine update_params (debug, num_params, 
     +                            gain, params, new_params )

	implicit real*8 (a-h, o-z)

	real*8		params(40), new_params(40)
	real*8		param_old(40)

	do j = 1, num_params

		param_old(j) = params(j)

		delta_param = new_params(j) - params(j) 

		sign = 1.d0
		if ( delta_param .lt. 0.0 ) sign = -1.d0

		trial_param = params(j) + gain * delta_param

		params(j) = trial_param

	enddo

	return
	end
c==============================================================
	subroutine calc_sigsq (num_data, num_params_solved_for,
     +                  std_err, resids, res_err, 
     +			sigsq, sigma_pdf )

C	Calculate new "sigsq"...

	implicit real*8	(a-h,o-z)

	real*8	 resids(20000), res_err(20000)


	sigsq = 0.d0
	do i_d = 1, num_data
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data - num_params_solved_for

	if ( num_deg_freedom .gt. 1 )
     +			sigsq_pdf = sigsq / float(num_deg_freedom) 
	sigma_pdf = dsqrt( sigsq_pdf )

c	write (6,1300) sigma_pdf, num_deg_freedom
c1300	format(/' Sigma_pdf =',f10.3,' for all data (',i6,
c     +	        ' degrees of freedom)',/1x)

	return
	end
c==============================================================
	subroutine weight_data ( num_data, delay_err, 
     +                           res_err)

C       Assigns data errors 

	implicit real*8 (a-h,o-z)

	real*8   res_err(1)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad


	n_geos = num_data 
	if ( n_geos .gt. 0 ) then

	   i_d = 0
	   do n = 1, n_geos

		 i_d = i_d + 1
		 res_err(i_d) = delay_err                        ! cm

	   enddo

	endif

	return
	end
c==============================================================
      subroutine calc_gmst0 ( day, gmst0 )

c     Calculates approximate Greenwich mean sidereal time at 0 UT
c     Input is Julian day number (ends in 0.5!)
c     Output is GMST0 in decimal hours

      implicit real*8 (a-h,o-z)

      day_number = day - 2451545.0d0
      JD = int( day_number )

      d = float(JD) + 0.5d0             ! enforce 0 UT

      t = d / 36525.d0

      gmst0 = 6.697374558d0 + 0.06570982441908d0 * d +
     +                        0.000026d0 * t

      gmst0 = mod ( gmst0, 24.d0 )                     ! hrs

      if ( gmst0 .lt. 0.d0 ) gmst0 = gmst0 + 24.d0

      return
      end
c==============================================================
	subroutine hmsrad ( hms, hms_rad )

C	Converts hhmmss.sss format to radians

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	xhms = abs(hms)

	ih = xhms/10000.d0 + 1.d-10
	im = dmod(xhms,10000.d0)/100.d0 + 1.d-10 
	s  = dmod(xhms,100.d0)

	hms_rad = dfloat(ih) + dfloat(im)/60.d0 + s/3600.d0
	hms_rad = hms_rad * pi / 12.d0
	if ( hms .lt. 0.d0 )  hms_rad = -hms_rad

	return
	end
c       ========================================================
	subroutine dmsrad ( dms, dms_rad )

c	Converts hhmmss.sss format to radians

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

	xdms = abs(dms)

	id = xdms/10000.d0
	im = dmod(xdms,10000.d0)/100.d0
	s  = dmod(xdms,100.d0)

	dms_rad = dfloat(id) + dfloat(im)/60.d0 + s/3600.d0
	dms_rad = dms_rad*deg_rad
	if ( dms .lt. 0.d0 )  dms_rad = -dms_rad

	return
	end
c       =========================================================
        subroutine radians_to_hhmmss ( ra_rad, ra_hhmmss)

c       Input :  ra_rad       in radians
c       Output:  ra_hhmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

        rad_to_hr = 12.d0/ pi

	ra_hr = ra_rad * rad_to_hr

        ihr  = ra_hr
        imin = (ra_hr - ihr*1.d0) * 60.d0
        sec  = (ra_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0

        ra_hhmmss  = ihr*10000.d0 + imin*100.d0 + sec

        return
        end
c       =========================================================
        subroutine radians_to_ddmmss ( dec_rad, dec_ddmmss)

c       Input :  dec_rad       in radians
c       Output:  dec_ddmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad

        rad_to_deg = 180.d0/ pi

	if ( dec_rad .ge. 0.d0 ) then
	     dec_sign =  1.d0
	  else
	     dec_sign = -1.d0
	endif

	dec_deg = abs( dec_rad ) / deg_rad

        ideg = dec_deg
        imin = (dec_deg - ideg*1.d0) * 60.d0
        sec  = (dec_deg - ideg*1.0d0 - imin/60.0d0) * 3600.d0

        dec_ddmmss  = ideg*10000.d0 + imin*100.d0 + sec
	dec_ddmmss  = dec_ddmmss * dec_sign

        return
        end
c==============================================================
        subroutine gauss_random ( x, u1,u2, g1,g2 )
 
C       Routine to return two uniformly distributed random numbers (0-1)
C       called "u1" and "u2", and two Gaussian distributed random numbers
C       called "g1" and "g2" with zero mean and unity sigma (I think).
 
C       Input a seed value of "x"...
C       a random number between 0 and 1 and ending
C       in an odd digit.  EG, seed =  0.76549 and not  0.76548
 
C       Every subsequent call should use the updated value of "x"
 
C       Uses the Box-Muller transformation

	implicit real*8 (a-h,o-z)
 
        twopi = 6.283185308d0
 
        x = 987.d0*x
        x = x - int(x)
        u1 = x
 
        x = 987.d0*x
        x = x - int(x)
        u2 = x
 
        amp = sqrt( -2.d0*log(u1) )
        arg = twopi*u2
        g1 = amp*cos(arg)
        g2 = amp*sin(arg)
         
        return
        end

c==============================================================
      SUBROUTINE get_gst (YEAR, MONTH, DAY, DATE, gmst_hr)

C Returns DATE = Julian Date and GMST = Greenwich Mean sidereal time
C at 0 hrs U.T. on day DAY/MONTH/YEAR. Accuracy not tested; GMST is
C correct at least to nearest second
C
C History:
C  1977 Aug  5 - TJP.
C  1991 May 18 - TJP.
c  2006 Sep 11 - MJR.
C-----------------------------------------------------------------------

      INTEGER YEAR
      INTEGER MONTH
      INTEGER DAY

      DOUBLE PRECISION DATE
      INTEGER MOFF(12), IC, NYRM1, JD, NDAYS
      DOUBLE PRECISION DU,SMD,T, gmst_hr
      DATA MOFF/0,31,59,90,120,151,181,212,243,273,304,334/
C
C JD number at 12 hrs UT on Jan 0 of year YEAR (Gregorian Calendar).
C
      NYRM1 = YEAR-1
      IC = NYRM1/100
      JD = 1721425 + 365*NYRM1 + NYRM1/4 - IC + IC/4
C
C Number of days from Standard Epoch 1900 Jan 0.5 (JD 2415020.0) to
C Jan 0.0 of YEAR.
C
      DU = DBLE(JD-2415020) - 0.5D0
C
C Day number; is it a leap year?
C
      NDAYS = MOFF(MONTH) + DAY
      IF (MONTH.GT.2 .AND. (
     1     ( MOD(YEAR,4).EQ.0 .AND. MOD(YEAR,100).NE.0 ) .OR.
     2     MOD(YEAR,400).EQ.0 ) ) NDAYS = NDAYS+1
C
C UT from Epoch to today (days), (centuries).
C
      SMD = DU+DBLE(NDAYS)
      T = SMD/36525.D0
C
C Greenwich Mean Sidereal Time.
C
      GMST = DBLE(6*3600 + 38*60) +45.836D0
     1     + 8 640 184.542D0*T  +  0.0929D0*T**2
      GMST = MOD(GMST,86400D0)

      gmst_hr = GMST / 3600.d0
C
C Julian Date.
C
      DATE = 2415020.D0+SMD
C
      return
      end
c==============================================================
      SUBROUTINE LEAST_SQUARES_FIT ( IDEBUG,print_solution,
     +                          PRINT_COR,MAXDAT,NUMXES,
     +				N,VNAME,X,ID,S,C,E,XHAT2,ERROR )

C     SUBROUTINE FOR LINEAR LEAST SQUARE FITTING...
C     IT REQUIRES INPUT OF THE FOLLOWING:
C        IDEBUG    = DEBUGGING PARAMETER
C                    0 = SUPRESS PRINT OUT BEYOND SOLUTION AND CORRELATI
C                    3 = PRINT OUT ALL DATA AND MATRICES
c        print_solution = logical flag: print out least-squares solution
C	 PRINT_COR = LOGICAL FLAG: PRINT CORRELATION MATRIX
C	 MAXDAT    = MAX NUM DATA POINTS.  DIMENSION LIMITS IN MAIN PROG.
C        NUMXES    = NUMBER OF PARAMETERS IN MODEL (WHETHER SOLVED FOR O
C        N         = NUMBER OF EQUATIONS
C        VNAME     = ALPHANUMERIC 'NAME' FOR THE PARAMETERS
C        X         = 1-D ARRAY OF INITIAL PARAMETER VALUES
C        ID        = 1-D ARRAY OF 'SOLVE FOR' CODES
C                  0 = HOLD PARAMETER CONSTANT
C                  1 = SOLVE FOR PARAMETER
C        S         = 1-D ARRAY OF DATA (E.G. VLB DIFFERENCED PHASES)
C        C         = 2-D ARRAY OF PARTIAL DERIVATIVES OF MODEL
C        E         = 1-D ARRAY OF ERRORS ON DATA POINTS
C
C     THE OUTPUT INCLUDES:
C        XHAT2     = 1-D ARRAY OF FINAL LEAST SQUARE VALUES FOR ALL PARA
C        THE SOLUTION AND CORRELATIONS ARE AUTOMATICALLY PRINTED OUT

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(MAXDAT,1)
      DIMENSION S(1),E(1),X(1),ID(1),XHAT2(1)

	character*8  VNAME(1),QUEST(2)
C
        REAL*8          B(40,40),BSTORE(40,40),XIDENT(40,40),
     +                  COR(40,40),ERROR(40),STDEV(40),PR(40),
     +                  XHAT(40),SNEW(20000),PRODCT(20000)

        DIMENSION LL(40),MM(40)
        DATA      MSIZE /40/

C
      LOGICAL PRINT_COR, print_solution

      DATA QUEST/'YES     ','NO      '/

 9999 FORMAT (1X,10(1PD13.5))
 9998 FORMAT (1X,20F6.2)
 9417 FORMAT (/1X)

      IF (IDEBUG.ge.3) then

	 WRITE (6,5522)
 5522	 FORMAT (////50X,' C MATRIX')
	 do I=1,N
	    WRITE (6,9999) (C(I,J),J=1,NUMXES)
	 enddo
	 WRITE (6,5523)
 5523	 FORMAT (////40X,'        S(I)           E(I)')
	 do I=1,N
	    WRITE (6,9999) S(I),E(I)
	 enddo

      endif

C     Count number of solve-for parameters (M)...
      M=0
      DO I=1,NUMXES
	 M=M+ID(I)
      enddo

c     Check number of data .ge. number of parameters...
      IF ( N.lt.M )  then
	 WRITE (6,2094)  N,M
 2094	 FORMAT (////' LEAST_SQUARES_FIT: # DATA POINTS (',I4,
     +               ') < # PARAMETERS (',I2,').  STOP' )
	 STOP
      endif

c     Compress partial derivative matrix to retain only those
c     associated with solve-for parameters...
      JNEW=0
      do J=1,NUMXES
	  if ( ID(J) .ne. 0 ) then
	     JNEW=JNEW+1
	     do I=1,N
   		 C(I,JNEW)=C(I,J)
	     enddo
	  endif
      enddo

C     Weight equations by dividing each by their error: E(I)
      do I=1,N
	  SNEW(I)=S(I)/E(I)
	  do J=1,M
	     C(I,J)=C(I,J)/E(I)
	  enddo
      enddo

C     Debug printout...
      IF (IDEBUG.ge.3) then
	 WRITE (6,2006)
 2006	 FORMAT ('1'/////50X,' CSTAR MATRIX')
	 do I=1,N
 	    WRITE (6,9999) (C(I,J),J=1,M)
	 enddo
	 WRITE (6,2009)
 2009	 FORMAT (/////50X,'     SSTAR')
	 do I=1,N
	    WRITE (6,9999) SNEW(I)
	 enddo
      endif

C     Square up partial deravitive matrix...
      ITEST2=0
      BMIN10=0.D0
      BMAX10=0.D0
      DO I=1,M
	DO J=1,M
	   B(I,J)=0.D0
	enddo
      enddo

      do I=1,M
	 do J=1,M
	      do L=1,N
		 B(I,J)=B(I,J) + C(L,I)*C(L,J)
              enddo
c             Test for big range of exponents...
	      if (B(I,J).ne.0.D0) then
		 B1=DLOG10( DABS( B(I,J) ) )
		 IF (BMIN10.GT.B1) BMIN10=B1
		 IF (BMAX10.LT.B1) BMAX10=B1
              endif
 	      BSTORE(I,J)=B(I,J)
	  enddo
      enddo

c     More debug printout...
      if (IDEBUG.ge.3) then
	 WRITE (6,2010)
 2010	 FORMAT ('1'/////50X,' C*TRANSPOSE C')
	 do I=1,M
	    WRITE (6,9999) (B(I,J),J=1,M)
	 enddo
      endif


c     Warn user if the matrix has a large range of exponents...
      IF (DABS(BMAX10-BMIN10).ge.16.D0) then
	 WRITE (6,9622) BMIN10,BMAX10
 9622	 FORMAT (///1X,'********   BMIN10 = ',F6.1,
     +          'BMAX10 = ',F6.1,'   *************',///1X)
      endif

c     "Center" exponents about zero...
      BADJST=10.D0**( (BMAX10+BMIN10)/2.D0 )
      do I=1,M
	 do J=1,M
 	    B(I,J)=B(I,J)/BADJST
	 enddo
      enddo

c     About to invert the matrix...
C        THE SUBROUTINE 'MINV' IS A DOUBLE PRECISION MATRIX INVERSION
C        IT INVERTS THE MATRIX 'BB' AND STORES IT AS 'BB' (I.E. THE ORIGIN
C        MATIRX IS DESTROYED

      CALL MINV(B,M,DETERM,LL,MM,MSIZE)

c     Re-scale inverted matrix for "Centering"
      do I=1,M
	 do J=1,M
	    B(I,J)=B(I,J)/BADJST
	 enddo
      enddo

c     More debug printout...
      IF (IDEBUG.ge.3) then
	 WRITE (6,2011) DETERM
 2011	 FORMAT ('1'////'  THE DETERMINANT IS',1PD13.5)
	 WRITE (6,2022)
 2022	 FORMAT (////45X,' (C*TRANSPOSE C*) INVERSE')
	 do I=1,M
 	    WRITE (6,9999) (B(I,J),J=1,M)
	 enddo
      endif

c     Prepare to check matrix inversion by multiplying inverse matrix
c     by the original martix to obtain (hopefully) the indentity matrix...
      do I=1,M
	 do J=1,M
	    XIDENT(I,J)=0.D0
	 enddo
      enddo

      do I=1,M
	 do J=1,M
	    do L=1,M
	       XIDENT(I,J)=XIDENT(I,J) + B(I,L)*BSTORE(L,J)
	    enddo
	    if (I.ne.J) then
c               Off-diagonal elements should be near zero...
		IF (DABS(XIDENT(I,J)).GT.1.D-06) ITEST2=1
	      else
c               Diagonal elements should be near unity...
		IF (DABS(XIDENT(I,I)-1.D0).GT.1.D-06) ITEST2=-1
	    endif
	 enddo
      enddo

      if (ITEST2.eq.0) then

c        Matrix seemed to invert.
c        Calculate corrections (XHAT's) to parameters...
	 do I=1,M
	    XHAT(I)=0.D0
            do J=1,N
	       PRODCT(J)=0.D0
	       do L=1,M
   		  PRODCT(J)=PRODCT(J) + B(I,L)*C(J,L)
	       enddo
	       XHAT(I)=XHAT(I)+PRODCT(J)*SNEW(J)
	    enddo
         enddo

C        XHAT'S ARE (IF ADDED TO THE X'S) THE UNCONSTRAINED LEAST SQUARE
C        SOLUTION OF THE 'SOLVE FOR' PARAMETERS ONLY.
C        XHAT2(J) = THE LEAST SQUARES SOLTIONS (INCLUDING NON SOLVE FOR
C        PARAMETERS)
	 IN=0
	 do J=1,NUMXES
	    if (ID(J).ne.0) then
	          IN=IN+1
	          XHAT2(J)=XHAT(IN) + X(J)
	       else
    		  XHAT2(J)=X(J)
	    endif
	 enddo

c        Calculate correlation coefficients...
	 do I=1,M
	    do J=1,M
	       COR(I,J)=-BSTORE(I,J)/DSQRT(BSTORE(I,I)*BSTORE(J,J))
	    enddo
	 enddo

c       Check for unphysical negative variances...
	INOTE=0
	do I=1,M
	      IF (B(I,I).le.0.D0) then
		 B(I,I)=-B(I,I)
		 INOTE = INOTE + 1
	      endif
	enddo
	if (INOTE.gt.0) then
	   WRITE (6,2071) INOTE
 2071	   FORMAT (///' ***** THERE ARE ',I2,' NEGATIVE VARIANCES; ',
     +            'SOLUTION IS UNPHYSICAL')
	endif

c       Convert variances to standard deviations...
	do J=1,M
	   STDEV(J)=DSQRT(B(J,J))
	enddo


C       REARRANGE CALCULATED 1 SIGMA ERRORS TO APPLY TO THE SOLVED FOR
C       PARAMETERS
	IN=0
	do J=1,NUMXES
	   if (ID(J).ne.0) then
	         IN=IN+1
		 ERROR(J)=STDEV(IN)
	      else
		 ERROR(J)=0.D0
	   endif
	enddo

        if ( print_solution ) then
C          OUTPUT FOR THE UNCONSTRAINED LEAST SQUARES SOLUTION
           WRITE (6,2040)
 2040	   FORMAT (/1X,'PARAMETER    ORIGINAL VALUE   LEAST SQUARE',
     +               ' VALUE  1 SIGMA ERRORS   SOLVED FOR?')

	   do J=1,NUMXES
	      L=1
	      IF (ID(J).EQ.0) L=2

	      if ( L.eq.1 ) then
 	         WRITE (6,2041) VNAME(J),X(J),XHAT2(J),ERROR(J),QUEST(L)
 2041	         FORMAT (2X,A8,5X,3(F13.6,5X),4X,A8)
	      endif

	   enddo
        endif

c       Print correlations?
	if ( PRINT_COR )  then
C	      CALCULATE THE CORRELATION COEFFICIENTS (COR)... 
	      do I=1,M
	         do J=1,M
		    COR(I,J)=B(I,J)/DSQRT(B(I,I)*B(J,J))
		 enddo
	      enddo

c	      WRITE (6,2056)
c 2056	      FORMAT(/10X,' THE CORRELATION COEFFICIENTS')
c	      do I=1,M
c   		  WRITE (6,9998) (COR(I,J),J=1,M)
c	      enddo

C  	      THE MULTIPLE CORRELATION COEFFICIENTS (PR) ARE CALCULATED
	      do I=1,M
		 PR(I)= 1.D0 - (1.d0/(BSTORE(I,I)*B(I,I)))
	      enddo
c	      WRITE (6,2060)
c 2060         FORMAT (///10X,'THE MULTIPLE CORRELATION COEFFICIENTS'/)
              I = 0
              do J=1,NUMXES
              	if ( ID(J).ne.0 ) then
              		I = I + 1
C                       Not sure if the true PR is the sqrt(PR)??
c                        WRITE (6,2061) VNAME(J),PR(I)
c 2061                   FORMAT (10X,A8,2X,F10.5)
		endif
	      enddo

         endif

       else

c         Print out warning...
	  WRITE (6,2095)  ITEST2
 2095	  FORMAT (////'  ***** ITEST2 =',I2,' ***** '
     +               /'       MATRIX INVERSION FAILED.')
	  STOP

      endif

      RETURN
      END
c==============================================================
      SUBROUTINE MINV(A,N,D,L,M,MSIZE)

      IMPLICIT REAL*8 (A-H,O-Z)
C     IBM SCIENTIFIC SUBROUTINE PACAKAGE PAGE 118
C     ******************************************************************
C
C     PURPOSE  INVERT A MATRIX
C
C     DESCRIPTION OF PARAMETERS
C        A  INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY INVER
C        N  ORDER OF MATRIX A
C        D  RESULTANT DETERMINANT
C        L  WORK VECTOR OF LENGTH N
C        M  WORK VECTOR OF LENGTH N
C     MSIZE ORDER OF TWO DIMENSIONAL MATRIX A IN MAIN PROGRAM
C
C     METHOD
C        THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT IS A
C        CALCULATED. A DETERMINANT OF ZERO INDICATES THAT THE MATRIX IS
C        SINGULAR.
C
C     ******************************************************************
C
      DIMENSION A(1),L(1),M(1)
C
C     STORE MATRIX A IN VECTOR FORM (COLUMN BY COLUMN)
C     SEE SSP PROGRAM ARRAY  P.98
C
      IF(MSIZE.EQ.N) GO TO 2
      NI=MSIZE-N
      IJ=0
      NM=0
      DO 230 K=1,N
      DO 225 LL=1,N
      IJ=IJ+1
      NM=NM+1
  225 A(IJ)=A(NM)
  230 NM=NM+NI
    2 CONTINUE
C
C     SEARCH FOR LARGEST ELEMENT
C
      D=1.0D0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
   10 IF(ABS(BIGA)-ABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C     INTERCHANGE ROWS
C
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C
C     INTERCHANGE COLUMNS
C
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
C
C     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS CONTAINED
C     BIGA)
C
   45 IF(BIGA) 48,46,48
   46 D=0.0D0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
C     REDUCE MATRIX
C
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
C
C     DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C     PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C     REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.0D0/BIGA
   80 CONTINUE
C
C     FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  100 K=K-1
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 CONTINUE
C
C     PUT MATRIX BACK INTO SQUARE FORM
C
      IF(MSIZE.EQ.N) GO TO 4
      IJ=N*N+1
      NM=N*MSIZE+1
      DO 210 K=1,N
      NM=NM-NI
      DO 210 LL=1,N
      IJ=IJ-1
      NM=NM-1
  210 A(NM)=A(IJ)
    4 CONTINUE
      RETURN
      END
