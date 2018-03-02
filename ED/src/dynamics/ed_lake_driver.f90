!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the fluxes between water and the air.  ED solves only   !
! the land regions, so this subroutine is a simplified version of LEAF-3 without consider- !
! ing land, only water.  It's called lake model but it is used for oceans and rivers too   !
! (using prescribed surface temperature.)                                                  !
!------------------------------------------------------------------------------------------!
subroutine simple_lake_model()
   use consts_coms            , only : stefan            & ! intent(in)
                                     , vonk              & ! intent(in)
                                     , grav              & ! intent(in)
                                     , rdry              & ! intent(in)
                                     , t00               & ! intent(in)
                                     , p00               & ! intent(in)
                                     , p00i              & ! intent(in)
                                     , rocp              & ! intent(in)
                                     , cpor              & ! intent(in)
                                     , epim1             & ! intent(in)
                                     , mmdryi            & ! intent(in)
                                     , cliq8             & ! intent(in)
                                     , tsupercool_liq8   & ! intent(in)
                                     , day_sec           & ! intent(in)
                                     , mmdry             ! ! intent(in)
   use met_driver_coms        , only : met_frq           & ! intent(in)
                                     , nformats          & ! intent(in)
                                     , met_nv            & ! intent(in)
                                     , met_interp        & ! intent(in)
                                     , met_vars          & ! intent(in)
                                     , met_driv_state    ! ! intent(in)
   use ed_lake_stepper        , only : integrate_lake    ! ! subroutine
   use ed_misc_coms           , only : current_time      & ! intent(in)
                                     , dtlsm             ! ! intent(in)
   use ed_state_vars          , only : edgrid_g          & ! intent(in)
                                     , edtype            & ! intent(in)
                                     , polygontype       & ! intent(in)
                                     , sitetype          & ! structure
                                     , patchtype         ! ! intent(in)
   use grid_coms              , only : ngrids            & ! intent(in)
                                     , nzg               & ! intent(in)
                                     , nzs               ! ! intent(in)
   use rk4_coms               , only : rk4min_can_shv    & ! intent(in)
                                     , rk4max_can_shv    & ! intent(in)
                                     , rk4min_can_temp   & ! intent(in)
                                     , rk4max_can_temp   & ! intent(in)
                                     , rk4min_can_co2    & ! intent(in)
                                     , rk4max_can_co2    & ! intent(in)
                                     , checkbudget       ! ! intent(in)
   use therm_lib              , only : rhovsil           & ! function
                                     , reducedpress      & ! function
                                     , thetaeiv          & ! function
                                     , idealdenssh       ! ! function
   use canopy_struct_dynamics , only : ed_stars          & ! subroutine
                                     , vertical_vel_flux ! ! function
   use lake_coms              , only : lake_buff         & ! intent(out)
                                     , dtlakei           & ! intent(in)
                                     , dtlake            & ! intent(in)
                                     , initial_lake_buff & ! subroutine
                                     , zero_lakesite     ! ! subroutine
   use budget_utils           , only : update_budget     & ! function
                                     , compute_budget    ! ! function
   implicit none
   !----- External functions --------------------------------------------------------------!
   logical                   , external    :: isleap

   !----- Local variables -----------------------------------------------------------------!
   integer                                 :: ifm, ipy, isi
   integer                                 :: i, iformat, iv
   integer                                 :: j
   integer                                 :: points_per_day, nday, np
   integer                                 :: nseconds_elapsed, ndays_elapsed
   integer                                 :: mnext
   integer                                 :: mprev
   logical                                 :: success
   real                                    :: timefac_sst
   real(kind=8)                            :: dsst_dt ! Derivative of lake surface temperature
   type(edtype)              , pointer     :: cgrid
   type(polygontype)         , pointer     :: cpoly
   type(sitetype)            , pointer     :: csite
   type(patchtype)           , pointer     :: cpatch
   type(met_driv_state)      , pointer     :: cmet
   real                                    :: wcurr_loss2atm
   real                                    :: ecurr_loss2atm
   real                                    :: co2curr_loss2atm
   real                                    :: wcurr_loss2drainage
   real                                    :: ecurr_loss2drainage
   real                                    :: wcurr_loss2runoff
   real                                    :: ecurr_loss2runoff
   real                                    :: ecurr_netrad
   real                                    :: old_can_theiv
   real                                    :: old_can_shv
   real                                    :: old_can_co2
   real                                    :: old_can_rhos
   real                                    :: old_can_temp
   real                                    :: old_can_prss
   integer                                 :: nsteps
   
   !----- Locally saved variables. --------------------------------------------------------!
   logical                   , save        :: first_time = .true. 
   real(kind=8)                            :: wfreeb              ! Free water
   real(kind=8)                            :: qwfree            ! Free water internal energy
   !---------------------------------------------------------------------------------------!


   !----- If this is the first time, nullify the structure. -------------------------------!
   if (first_time) then
      call initial_lake_buff()
      first_time = .false.
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the index of tmp and write down it's time step                               !
   !---------------------------------------------------------------------------------------!
   success=.false.
   do i=1, nformats
      do j=1, met_nv(i)
         if (trim(met_vars(i,j)) == 'tmp') then
            iformat=i
            iv=j
            success=.true.
            exit
         end if
      end do
      if (success) exit
   end do
   
   if (.NOT. success) then
      call fatal_error('Cannot find tmp from the meteorological input!','simple_lake_model'&
                                        ,'ed_lake_driver.f90')
   end if
   !------------------------------------------------------------------------------------------!
   
   !------------------------------------------------------------------------------------------!
   !    Big domain loops start here.                                                          !
   !------------------------------------------------------------------------------------------!
   gridloop: do ifm=1,ngrids
      cgrid => edgrid_g(ifm)
      
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         select case(met_interp(iformat,iv))
         case (2,4)
            !---------------------------------------------------------------------------------!
            !     Time independent variables, set timefac_sst to 0.                           !
            !---------------------------------------------------------------------------------!
            timefac_sst = 0.
            dsst_dt     = 0.
         case default
            !----- First, get the number of points per day and per month. --------------------!
            points_per_day = nint(day_sec/met_frq(iformat,iv))
            select case (current_time%month)
            case (1,3,5,7,8,10,12)
               nday = 31
            case (4,6,9,11)
               nday = 30
            case (2)
               if (isleap(current_time%year)) then
                  nday = 29
               else
                  nday = 28
               end if
            end select
            np = nday * points_per_day
            !---------------------------------------------------------------------------------!
            
            !----- Find the number of seconds since the beginning of this month. -------------!
            ndays_elapsed    = current_time%date - 1
            nseconds_elapsed = nint(current_time%time) + ndays_elapsed * nint(day_sec)
            !---------------------------------------------------------------------------------!

            !---------------------------------------------------------------------------------!
            ! mprev  - Index of the element of the previous metinput data.                    !
            ! mnext  - Index of the element of the next metinput data.                        !
            !---------------------------------------------------------------------------------!
            mprev  = int(real(nseconds_elapsed)/met_frq(iformat,iv)) + 1
            mnext  = mprev + 1
            !---------------------------------------------------------------------------------!
            
            !---------------------------------------------------------------------------------!
            !     For the particular case where we at the last day of the month, after the    !
            ! last met driver data of that month, the next time is the first day of the       !
            ! next month.  This will always be in the next position after day 31, even        !
            ! when the current month is shorter.                                              !
            !---------------------------------------------------------------------------------!
            if(mnext > np)then
               mnext = 1 + nint(day_sec / met_frq(iformat,iv)) * 31
            endif            
            !---------------------------------------------------------------------------------!
            
            timefac_sst = dtlsm / met_frq(iformat,iv)
            dsst_dt     = ( cgrid%metinput(ipy)%tmp(mnext) - cgrid%metinput(ipy)%tmp(mprev)&
                          ) / met_frq(iformat,iv) / 5.0
         end select


         siteloop: do isi = 1,cpoly%nsites
            if (cpoly%islakesite(isi) /= 1) cycle siteloop
            
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)
            
            !----- Save the previous thermodynamic state. ------------------------------------!
            old_can_theiv    = csite%can_theiv(1)
            old_can_shv      = csite%can_shv(1)
            old_can_co2      = csite%can_co2(1)
            old_can_rhos     = csite%can_rhos(1)
            old_can_temp     = csite%can_temp(1)
            old_can_prss     = csite%can_prss(1)
            !---------------------------------------------------------------------------------!



            !---------------------------------------------------------------------------------!
            !     Flush all buffers to zero, for a clean start.                               !
            !---------------------------------------------------------------------------------!
            call zero_lakesite(lake_buff%initp)
            call zero_lakesite(lake_buff%yscal)
            call zero_lakesite(lake_buff%y    )
            call zero_lakesite(lake_buff%dydx )
            call zero_lakesite(lake_buff%yerr )
            call zero_lakesite(lake_buff%ytemp)
            call zero_lakesite(lake_buff%ak2  )
            call zero_lakesite(lake_buff%ak3  )
            call zero_lakesite(lake_buff%ak4  )
            call zero_lakesite(lake_buff%ak5  )
            call zero_lakesite(lake_buff%ak6  )
            call zero_lakesite(lake_buff%ak7  )
            !---------------------------------------------------------------------------------!



            !----- Set up the meteorological forcing structure. ------------------------------!
            call copy_met_2_lakesite(cmet,dsst_dt,cgrid%cosz(ipy),cgrid%lat(ipy),             &
                 cgrid%lon(ipy),ifm,ipy)
            !---------------------------------------------------------------------------------!

            !----- Compute current storage terms. --------------------------------------------!
            call update_budget(csite,cpoly%lsl(isi),1,1)

            !---------------------------------------------------------------------------------!
            !     Update lake/ground temperature.                                             !
            !---------------------------------------------------------------------------------!
            csite%ground_temp(1) = csite%ground_temp(1) + dsst_dt*dtlsm
            !---------------------------------------------------------------------------------!



            !---------------------------------------------------------------------------------!
            !     Assign the initial values for "canopy" properties.                          !
            !---------------------------------------------------------------------------------!
            call copy_lake_init(csite,lake_buff%initp)
            !---------------------------------------------------------------------------------!



            !---------------------------------------------------------------------------------!
            !     Integrate the state variables and fluxes.                                   !
            !---------------------------------------------------------------------------------!
            call integrate_lake(dtlsm,csite%htry(1),nsteps,csite)
            !---------------------------------------------------------------------------------!
            
            
            
            !---------------------------------------------------------------------------------!
            !        Normalize canopy-atmosphere flux values                                  !
            !---------------------------------------------------------------------------------!
            lake_buff%initp%upwp = lake_buff%initp%can_rhos * lake_buff%initp%upwp * dtlakei
            lake_buff%initp%tpwp = lake_buff%initp%can_rhos * lake_buff%initp%tpwp * dtlakei
            lake_buff%initp%qpwp = lake_buff%initp%can_rhos * lake_buff%initp%qpwp * dtlakei
            lake_buff%initp%cpwp = lake_buff%initp%can_rhos * lake_buff%initp%cpwp * dtlakei
            lake_buff%initp%wpwp = lake_buff%initp%can_rhos * lake_buff%initp%wpwp * dtlakei

            !----- Compute runoff for output -------------------------------------------------!
            !      Because it is lake, all rainfall becomes runoff.                           !
            wfreeb                  = cmet%pcpg  * dtlake
            qwfree                  = cmet%qpcpg * dtlake !wfreeb * cliq8 * (csite%ground_temp(1)-tsupercool_liq8)
            csite%runoff(1)         = wfreeb * dtlakei
            csite%avg_runoff(1)     = csite%avg_runoff(1) + wfreeb
            csite%avg_runoff_heat(1) = csite%avg_runoff_heat(1) + qwfree
            
            lake_buff%initp%wbudget_loss2runoff = wfreeb
            lake_buff%initp%ebudget_loss2runoff = qwfree
            lake_buff%initp%wbudget_storage     = lake_buff%initp%wbudget_storage - wfreeb
            lake_buff%initp%ebudget_storage     = lake_buff%initp%ebudget_storage - qwfree
            lake_buff%initp%wbudget_loss2drainage = 0.
            lake_buff%initp%ebudget_loss2drainage = 0.
            
            !---------------------------------------------------------------------------------!
            !     Copy the state variables and the fluxes to the site output arrays.          !
            !---------------------------------------------------------------------------------!
            call initp2lakesite(csite,nzg,nzs,lake_buff%initp,wcurr_loss2atm,ecurr_loss2atm   &
                       ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage              &
                       ,wcurr_loss2runoff,ecurr_loss2runoff)
            !---------------------------------------------------------------------------------!
            
            !---------------------------------------------------------------------------------!
            !     Make final adjustment to lake temperature and energy.                       !
            !---------------------------------------------------------------------------------!
            call final_adjust_ebudget(csite,cpoly%lsl(isi),cmet%qpcpg,1,ecurr_netrad          &
                                     ,ecurr_loss2atm,ecurr_loss2drainage ,ecurr_loss2runoff   &
                                     ,old_can_theiv,old_can_shv,old_can_co2,old_can_rhos      &
                                     ,old_can_temp)
            !---------------------------------------------------------------------------------!
            
            
            !----- Add the number of steps into the step counter. ----------------------------!
            cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)
            
            
            !---------------------------------------------------------------------------------!
            !    Update the minimum monthly temperature, based on canopy temperature.         !
            !---------------------------------------------------------------------------------!
            if (cpoly%site(isi)%can_temp(1) < cpoly%min_monthly_temp(isi)) then
               cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(1)
            end if



            !---------------------------------------------------------------------------!
            !     Compute the residuals.                                                !
            !---------------------------------------------------------------------------!
            if (.not. checkbudget) then
               co2curr_loss2atm             = 0.
               ecurr_loss2atm               = 0.
               ecurr_loss2drainage          = 0.
               ecurr_loss2runoff            = 0.
               wcurr_loss2atm               = 0.
               wcurr_loss2drainage          = 0.
               wcurr_loss2runoff            = 0.
            end if
            
            call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,1             &
                               ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm              &
                               ,co2curr_loss2atm,wcurr_loss2drainage                    &
                               ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff &
                               ,cpoly%area(isi),cgrid%cbudget_nep(ipy),old_can_theiv    &
                               ,old_can_shv,old_can_co2,old_can_rhos,old_can_temp       &
                               ,old_can_prss)
         end do siteloop
      end do polyloop
   end do gridloop

   return
end subroutine simple_lake_model
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This subroutine will initialise the surface fields for the simple lake model.        !
!------------------------------------------------------------------------------------------!
subroutine copy_met_2_lakesite(cmet,dsst_dt,cosz,lat,lon,ifm,ipy)
   use therm_lib8             , only : thetaeiv8         & ! function
                                     , idealdenssh8      & ! function
                                     , rehuil8           ! ! function
   use lake_coms              , only : lakemet           ! ! intent(out)
   use consts_coms            , only : cpdryi8           & ! intent(in)
                                     , p00i8             & ! intent(in)
                                     , p008              & ! intent(in)
                                     , cpor8             ! ! intent(in)
   use canopy_air_coms        , only : ubmin8            ! ! intent(in)
   use met_driver_coms        , only : met_driv_state    ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(met_driv_state)      , target          :: cmet
   real(kind=8)              , intent(in)      :: dsst_dt
   real                      , intent(inout)   :: cosz
   real                      , intent(in)      :: lat
   real                      , intent(in)      :: lon
   integer                   , intent(in)      :: ifm
   integer                   , intent(in)      :: ipy

   !----- External functions. -------------------------------------------------------------!
   logical                   , external :: is_finite
   logical                   , external :: is_finite8
   !---------------------------------------------------------------------------------------!

   !----- Local variables -----------------------------------------------------------------!
   logical                              :: ok_flpoint



   !---------------------------------------------------------------------------------------!
   !     Copy the values to the meteorological buffer.  Start with those that require only !
   ! a simple copy with the conversion from single to double precision.                    !
   !---------------------------------------------------------------------------------------!
   cosz              = min(cosz, 0.001)       
   lakemet%rshort    = dble(cmet%rshort)      
   lakemet%rlong     = dble(cmet%rlong)       
   lakemet%tanz      = tan(acos(dble(cosz)))  
   lakemet%lon       = dble(lon)              
   lakemet%lat       = dble(lat)              
   
   lakemet%geoht     = dble(cmet%geoht    )
   lakemet%atm_theta = dble(cmet%atm_theta)
   lakemet%atm_co2   = dble(cmet%atm_co2  )
   lakemet%atm_exner = dble(cmet%exner    )
   lakemet%atm_shv   = dble(cmet%atm_shv  )
   lakemet%atm_vels  = dble(cmet%vels )
   lakemet%atm_rvap  = dble(cmet%atm_shv/(1.0-cmet%atm_shv))
   lakemet%pcpg      = dble(cmet%pcpg     )
   lakemet%qpcpg     = dble(cmet%qpcpg    )
   lakemet%dpcpg     = dble(cmet%dpcpg    )
   !----- SST derivative is already in double precision, just copy it. --------------------!
   lakemet%dsst_dt   = dsst_dt
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Sanity check.                                                                     !
   !---------------------------------------------------------------------------------------!
   ok_flpoint = is_finite8(lakemet%rshort)    .and. is_finite8(lakemet%rlong)    .and.     &
                is_finite8(lakemet%tanz)      .and. is_finite8(lakemet%lon)      .and.     &
                is_finite8(lakemet%lat)       .and. is_finite8(lakemet%geoht)    .and.     &
                is_finite8(lakemet%atm_theta) .and. is_finite8(lakemet%atm_co2)  .and.     &
                is_finite8(lakemet%atm_exner) .and. is_finite8(lakemet%atm_rvap) .and.     &
                is_finite8(lakemet%dsst_dt)   .and. is_finite8(lakemet%atm_shv)  .and.     &
                is_finite8(lakemet%atm_vels)
                
   if (.not. ok_flpoint) then
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) '  Something went wrong...                        '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' Meteorological forcing.'
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rshort         :',lakemet%rshort
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rlong          :',lakemet%rlong
      write(unit=*,fmt='(a,1x,es12.5)') ' - Tanz           :',lakemet%tanz
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lon            :',lakemet%lon
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lat            :',lakemet%lat
      write(unit=*,fmt='(a,1x,es12.5)') ' - Geoht          :',lakemet%geoht
      write(unit=*,fmt='(a,1x,es12.5)') ' - Theta          :',lakemet%atm_theta
      write(unit=*,fmt='(a,1x,es12.5)') ' - CO2            :',lakemet%atm_co2
      write(unit=*,fmt='(a,1x,es12.5)') ' - Exner          :',lakemet%atm_exner
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rvap           :',lakemet%atm_rvap
      write(unit=*,fmt='(a,1x,es12.5)') ' - Shv            :',lakemet%atm_shv
      write(unit=*,fmt='(a,1x,es12.5)') ' - Vels           :',lakemet%atm_vels
      write(unit=*,fmt='(a,1x,es12.5)') ' - d(SST)/dt      :',lakemet%dsst_dt
      write(unit=*,fmt='(a,1x,i12)'   ) ' - Grid           :',ifm   
      write(unit=*,fmt='(a,1x,i12)'   ) ' - Polygon        :',ipy   
      write(unit=*,fmt='(a,1x,es12.5)') ' - Glon           :',lon
      write(unit=*,fmt='(a,1x,es12.5)') ' - Glat           :',lat
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      call fatal_error('Non-resolvable values','copy_met_2_lake','ed_lake_misc.f90')
   end if




   !----- Log of potential temperature. ---------------------------------------------------!
   lakemet%atm_lntheta  = log(lakemet%atm_theta)
   !---------------------------------------------------------------------------------------!



   !----- Pressure. -----------------------------------------------------------------------!
   lakemet%atm_prss  = (lakemet%atm_exner * cpdryi8) ** cpor8 * p008
   !---------------------------------------------------------------------------------------!


   !----- Air temperature. ----------------------------------------------------------------!
   lakemet%atm_tmp   = cpdryi8 * lakemet%atm_theta * lakemet%atm_exner
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Update properties that need to use therm_lib8.                                    !
   !---------------------------------------------------------------------------------------!
   lakemet%atm_theiv = thetaeiv8(lakemet%atm_theta,lakemet%atm_prss,lakemet%atm_tmp        &
                                ,lakemet%atm_rvap,lakemet%atm_rvap)
   lakemet%atm_rhos  = idealdenssh8(lakemet%atm_prss,lakemet%atm_tmp,lakemet%atm_shv)
   lakemet%atm_rhv   = rehuil8(lakemet%atm_prss,lakemet%atm_tmp,lakemet%atm_rvap,.false.)

   return
end subroutine copy_met_2_lakesite
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!     This subroutine will copy the variables from the integration buffer to the state     !
! site and patch.                                                                          !
!==========================================================================================!
subroutine initp2lakesite(csite,mzg,mzs,initp,wbudget_loss2atm,ebudget_loss2atm            &
                          ,co2budget_loss2atm,wbudget_loss2drainage,ebudget_loss2drainage  &
                          ,wbudget_loss2runoff,ebudget_loss2runoff)
   use lake_coms             , only : lakemet              ! ! intent(out)
   use consts_coms           , only : cliqvlme8            & ! intent(in)
                                    , tsupercool_liq8      & ! intent(in)
                                    , cliq                 & ! intent(in)
                                    , grav                 ! ! intent(in)
   use canopy_air_coms       , only : ubmin8               ! ! intent(in)
   use lake_coms             , only : lakesitetype         & ! structure
                                    , lakemet              & ! intent(in)
                                    , tiny_lakeoff         ! ! intent(in)
   use ed_misc_coms          , only : fast_diagnostics     & ! intent(in)
                                    , dtlsm                ! ! intent(in)
   use ed_state_vars         , only : sitetype             ! ! structure
   
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)    , target      :: csite
   integer           , intent(in)  :: mzg
   integer           , intent(in)  :: mzs
   type(lakesitetype), target      :: initp
   real              , intent(out) :: wbudget_loss2atm
   real              , intent(out) :: ebudget_loss2atm
   real              , intent(out) :: co2budget_loss2atm
   real              , intent(out) :: wbudget_loss2drainage
   real              , intent(out) :: ebudget_loss2drainage
   real              , intent(out) :: wbudget_loss2runoff
   real              , intent(out) :: ebudget_loss2runoff
   !----- Local variables. ----------------------------------------------------------------!
   integer                        :: k
   !----- External functions. -------------------------------------------------------------!
   real              , external   :: sngloff
   !----- Local constants -----------------------------------------------------------------!
   real             , parameter   :: z0fac_water  = 0.016/grav
   real             , parameter   :: z0_min_water = 0.0001
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Most variables require just a simple copy.  More comments will be made next to    ! 
   ! those in which this is not true.  All floating point variables are converted back     !
   ! to single precision.                                                                  !
   !---------------------------------------------------------------------------------------!
   csite%can_theiv       (1)  = sngloff(initp%can_theiv          ,tiny_lakeoff)
   csite%can_theta       (1)  = sngloff(initp%can_theta          ,tiny_lakeoff)
   csite%can_prss        (1)  = sngloff(initp%can_prss           ,tiny_lakeoff)
   csite%can_temp        (1)  = sngloff(initp%can_temp           ,tiny_lakeoff)
   csite%can_shv         (1)  = sngloff(initp%can_shv            ,tiny_lakeoff)
   csite%can_co2         (1)  = sngloff(initp%can_co2            ,tiny_lakeoff)
   csite%can_rhos        (1)  = sngloff(initp%can_rhos           ,tiny_lakeoff)
   csite%can_depth       (1)  = sngloff(initp%can_depth          ,tiny_lakeoff)
   csite%veg_displace    (1)  = 0.0
   csite%rough           (1)  = sngloff(initp%lake_rough         ,tiny_lakeoff)
   csite%snowfac         (1)  = 1.0-initp%lake_fliq
   csite%total_sfcw_depth(1)  = 0.0 
   csite%ggsoil          (1)  = sngloff(initp%gglake             ,tiny_lakeoff)
   csite%ggbare          (1)  = 0.0 
   csite%ggveg           (1)  = 0.0 
   csite%ggnet           (1)  = 0.0 
   
   csite%ustar           (1)  = sngloff(initp%ustar              ,tiny_lakeoff)
   csite%tstar           (1)  = sngloff(initp%tstar              ,tiny_lakeoff)
   csite%qstar           (1)  = sngloff(initp%qstar              ,tiny_lakeoff)
   csite%cstar           (1)  = sngloff(initp%cstar              ,tiny_lakeoff)
   
   csite%zeta            (1)  = sngloff(initp%zeta               ,tiny_lakeoff)
   csite%ribulk          (1)  = sngloff(initp%ribulk             ,tiny_lakeoff)
   
   csite%upwp            (1)  = sngloff(initp%upwp               ,tiny_lakeoff)
   csite%wpwp            (1)  = sngloff(initp%wpwp               ,tiny_lakeoff)
   csite%tpwp            (1)  = sngloff(initp%tpwp               ,tiny_lakeoff)
   csite%qpwp            (1)  = sngloff(initp%qpwp               ,tiny_lakeoff)
   csite%cpwp            (1)  = sngloff(initp%cpwp               ,tiny_lakeoff)

   !---------------------------------------------------------------------------------------!


   !----- Fluxes towards the free atmosphere. ---------------------------------------------!
   !ed_fluxf_g(ifm)%sflux_u(i,j,1) = sngloff(initp%avg_sflux_u ,tiny_lakeoff)
   !ed_fluxf_g(ifm)%sflux_v(i,j,1) = sngloff(initp%avg_sflux_v ,tiny_lakeoff)
   !ed_fluxf_g(ifm)%sflux_w(i,j,1) = sngloff(initp%avg_sflux_w ,tiny_lakeoff)
   !ed_fluxf_g(ifm)%sflux_t(i,j,1) = sngloff(initp%avg_sflux_t ,tiny_lakeoff)
   !ed_fluxf_g(ifm)%sflux_r(i,j,1) = sngloff(initp%avg_sflux_r ,tiny_lakeoff)
   !ed_fluxf_g(ifm)%sflux_c(i,j,1) = sngloff(initp%avg_sflux_c ,tiny_lakeoff)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    These variables are fast scale fluxes, and they may not be allocated, so just      !
   ! check this before copying.                                                            !
   !---------------------------------------------------------------------------------------!
   if (fast_diagnostics) then
      csite%avg_vapor_lc        (1) = 0.0
      csite%avg_vapor_wc        (1) = 0.0
      csite%avg_vapor_gc        (1) = sngloff(initp%avg_vapor_gc       ,tiny_lakeoff)
      csite%avg_wshed_vg        (1) = 0.0
      csite%avg_intercepted     (1) = 0.0
      csite%avg_throughfall     (1) = csite%avg_throughfall(1) + lakemet%pcpg * dtlsm
      csite%avg_vapor_ac        (1) = sngloff(initp%avg_vapor_ac       ,tiny_lakeoff)
      csite%avg_transp          (1) = 0.0
      csite%avg_evap            (1) = sngloff(initp%avg_vapor_ac           ,tiny_lakeoff)
      csite%avg_drainage        (1) = 0.0
      csite%avg_drainage_heat   (1) = 0.0
      csite%avg_rshort_gnd      (1) = sngloff(initp%avg_rshort_gnd     ,tiny_lakeoff)
      csite%avg_rlong_gnd       (1) = sngloff(initp%avg_rlong_gnd      ,tiny_lakeoff)
      csite%avg_sensible_lc     (1) = 0.0
      csite%avg_sensible_wc     (1) = 0.0
      csite%avg_qwshed_vg       (1) = 0.0
      csite%avg_qintercepted    (1) = 0.0
      csite%avg_qthroughfall    (1) = csite%avg_qthroughfall(1) + lakemet%qpcpg * dtlsm
      csite%avg_sensible_gc     (1) = sngloff(initp%avg_sensible_gc    ,tiny_lakeoff)
      csite%avg_sensible_ac     (1) = sngloff(initp%avg_sensible_ac    ,tiny_lakeoff)
      csite%avg_carbon_ac       (1) = 0.0
      csite%avg_carbon_st       (1) = 0.0
      csite%avg_ustar           (1) = sngloff(initp%ustar              ,tiny_lakeoff)
      csite%avg_tstar           (1) = sngloff(initp%tstar              ,tiny_lakeoff)
      csite%avg_qstar           (1) = sngloff(initp%qstar              ,tiny_lakeoff)
      csite%avg_cstar           (1) = sngloff(initp%cstar              ,tiny_lakeoff)
      do k = 1, mzg
         csite%avg_sensible_gg(k,1) = 0.0
         csite%avg_smoist_gg  (k,1) = 0.0
         csite%avg_transloss  (k,1) = 0.0
      end do

      !---------------------------------------------------------------------------------!
      !     These variables are integrated here, since they don't change with time.     !
      !---------------------------------------------------------------------------------!
      csite%avg_rlongup        (1) = csite%avg_rlongup        (1)                   &
                                     + csite%rlongup          (1) * dtlsm 
      csite%avg_albedo         (1) = csite%avg_albedo         (1)                   &
                                     + csite%albedo           (1) * dtlsm
      csite%avg_albedo_beam    (1) = csite%avg_albedo_beam    (1)                   &
                                     + csite%albedo_beam      (1) * dtlsm
      csite%avg_albedo_diffuse (1) = csite%avg_albedo_diffuse (1)                   &
                                     + csite%albedo_diffuse   (1) * dtlsm
      csite%avg_rlong_albedo   (1) = csite%avg_rlong_albedo   (1)                   &
                                     + csite%rlong_albedo     (1) * dtlsm
      !---------------------------------------------------------------------------------!
   end if
   !------------------------------------------------------------------------------------!


    co2budget_loss2atm    = sngloff(initp%co2budget_loss2atm   ,tiny_lakeoff)
    ebudget_loss2atm      = sngloff(initp%ebudget_loss2atm     ,tiny_lakeoff)
    ebudget_loss2drainage = sngloff(initp%ebudget_loss2drainage,tiny_lakeoff)
    ebudget_loss2runoff   = sngloff(initp%ebudget_loss2runoff  ,tiny_lakeoff)
    wbudget_loss2atm      = sngloff(initp%wbudget_loss2atm     ,tiny_lakeoff)
    wbudget_loss2drainage = sngloff(initp%wbudget_loss2drainage,tiny_lakeoff)
    wbudget_loss2runoff   = sngloff(initp%wbudget_loss2runoff  ,tiny_lakeoff)

   !------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   csite%nlev_sfcwater(1) = 0.
   do k=1, mzs
      csite%sfcwater_depth  (k,1)   = 0.
      csite%sfcwater_mass   (k,1)   = 0.
      csite%sfcwater_tempk  (k,1)   = 0.
      csite%sfcwater_fracliq(k,1)   = initp%lake_fliq
      csite%sfcwater_energy (k,1)   = 0.
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    "Soil" energy.  Because we can't store sea surface temperature, we store the       !
   ! internal energy. The "soil" energy will be updated in the below final_adjust_ebudget. !         
   !---------------------------------------------------------------------------------------!

   return
end subroutine initp2lakesite
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
! Because this is a simple lake modeland we don't have inputs of lake surface temperature, !
! we use an "energy balance closure" method to determine the lake/ground temperature.      !
!==========================================================================================!
subroutine final_adjust_ebudget(csite,lsl,qpcpg,ipa,ecurr_netrad,ecurr_loss2atm            &
                         ,ecurr_loss2drainage,ecurr_loss2runoff,old_can_theiv,old_can_shv  &
                         ,old_can_co2,old_can_rhos,old_can_temp)
   use ed_state_vars         , only : sitetype           ! ! structure   
   use ed_misc_coms          , only : dtlsm              ! ! intent(in)  
   use consts_coms           , only : cpdry              & ! intent(in)  
                                    , cliqvlme8          & ! intent(in)  
                                    , tsupercool_liq8    & ! intent(in)  
                                    , allivlme           & ! intent(in)  
                                    , cicevlme           & ! intent(in)  
                                    , wdns               ! ! intent(in)  
   use grid_coms             , only : nzg                & ! intent(in)  
                                    , nzs                ! ! intent(in)
   use lake_coms             , only : tiny_lakeoff       & ! intent(in)  
                                    , lakemet            ! ! intent(in)
   use soil_coms             , only : dslz               ! ! intent(in)  
   use therm_lib             , only : uextcm2tl          ! ! subroutine
   use rk4_coms              , only : rk4min_soil_temp   & ! intent(in)
                                    , rk4max_soil_temp   ! ! intent(in)
   use budget_utils          , only : compute_netrad     & ! function
                                    , ddens_dt_effect    ! ! function
                                    
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target        :: csite
   real                  , intent(in)    :: qpcpg
   real                  , intent(inout) :: ecurr_netrad
   real                  , intent(in)    :: ecurr_loss2atm
   real                  , intent(in)    :: ecurr_loss2drainage
   real                  , intent(in)    :: ecurr_loss2runoff
   integer               , intent(in)    :: lsl
   integer               , intent(in)    :: ipa
   real                  , intent(in)    :: old_can_theiv
   real                  , intent(in)    :: old_can_shv
   real                  , intent(in)    :: old_can_co2
   real                  , intent(in)    :: old_can_rhos
   real                  , intent(in)    :: old_can_temp
   !----- Local variables -----------------------------------------------------------------!
   real                                  :: ebudget_finalstorage
   real                                  :: ebudget_deltastorage
   real                                  :: ecurr_precipgain
   real                                  :: ecurr_denseffect
   real                                  :: ecurr_residual
   real                                  :: old_can_rhotemp
   real                                  :: curr_can_rhotemp
   real                                  :: old_can_lntheiv
   real                                  :: curr_can_lntheiv
   real                                  :: ecurr_soilstorage
   real                                  :: ecurr_casstorage
   real                                  :: depth
   real                                  :: lake_temp
   real                                  :: fracliq
   integer                               :: k
   !----- External functions. -------------------------------------------------------------!
   real                  , external      :: sngloff
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Compute gain in energy due to precipitation.                            !
   !---------------------------------------------------------------------------------------!
   ecurr_precipgain = qpcpg * dtlsm

   !---------------------------------------------------------------------------------------!
   !     Compute gain in energy due to radiation.                                          !
   !---------------------------------------------------------------------------------------!
   ecurr_netrad     = compute_netrad(csite,ipa) * dtlsm


   !---------------------------------------------------------------------------------------!
   !     For enthalpy, we must consider both density and temperature effects.              !
   !---------------------------------------------------------------------------------------!
   old_can_rhotemp  = old_can_rhos        * old_can_temp
   curr_can_rhotemp = csite%can_rhos(ipa) * csite%can_temp(ipa)
   old_can_lntheiv  = log(old_can_theiv)
   curr_can_lntheiv = log(csite%can_theiv(ipa))
   ecurr_denseffect    = ddens_dt_effect(old_can_rhotemp,curr_can_rhotemp                  &
                                        ,old_can_lntheiv,curr_can_lntheiv                  &
                                        ,csite%can_depth(ipa),cpdry)

   !----- Compute current cas storage term. -----------------------------------------------!
   ecurr_casstorage       = cpdry * csite%can_rhos(ipa) * csite%can_depth(ipa) *           &
                            csite%can_theiv(ipa)
                            
   !----- Determine current soil storage term according to energy budget closure. ---------!
   ebudget_finalstorage   = csite%ebudget_initialstorage(ipa) + ( ecurr_precipgain +       &
                            ecurr_netrad - ecurr_loss2atm - ecurr_loss2drainage -          &
                            ecurr_loss2runoff ) + ecurr_denseffect
   ecurr_soilstorage      = ebudget_finalstorage - ecurr_casstorage


   depth = 0.0
   do k = lsl, nzg
      depth = depth + dslz(k)
   end do

   ecurr_soilstorage      = ecurr_soilstorage / depth

   !----- Finding the average temperature and liquid fraction. ----------------------------!
   call uextcm2tl(ecurr_soilstorage,wdns,0.0,lake_temp,fracliq)
 
   !----- Sanity check --------------------------------------------------------------------!
   if ( lake_temp < sngl(rk4min_soil_temp) .or. lake_temp > sngl(rk4max_soil_temp) ) then
      write (unit=*,fmt='(80a)')         ('=',k=1,80)
      write (unit=*,fmt='(a)')           'LAKE_TEMP IS WRONG IN FINAL_ADJUST_EBUDGET' 
      write (unit=*,fmt='(80a)')         ('-',k=1,80)                                 
      write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',lakemet%lon           
      write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',lakemet%lat           
      write (unit=*,fmt='(a,1x,es12.4)') ' + LAKE_TEMP:    ',lake_temp           
      write (unit=*,fmt='(80a)') ('-',k=1,80)                                         
      call fatal_error('extreme lake temperature','final_adjust_ebudget'                   &
                      &,'ed_lake_driver.f90')                                             
   end if
  
   csite%soil_water  (nzg,ipa)    = 1.0
   csite%soil_fracliq(nzg,ipa)    = fracliq
   csite%soil_energy (nzg,ipa)    = sngloff( (1.0-fracliq)*cicevlme*dble(lake_temp)        &
                                           + fracliq *  cliqvlme8 * ( dble(lake_temp)      &
                                           - tsupercool_liq8 ) ,tiny_lakeoff)
  
   csite%soil_tempk  (nzg,ipa)    = sngloff(dble(lake_temp) , tiny_lakeoff)                        
   csite%ground_temp (ipa    )    = csite%soil_tempk  (nzg,ipa)
   
   do k=1, nzg-1
      csite%soil_water  (k,ipa)   = csite%soil_water  (nzg,ipa)
      csite%soil_fracliq(k,ipa)   = csite%soil_fracliq(nzg,ipa)
      csite%soil_energy (k,ipa)   = csite%soil_energy (nzg,ipa)
      csite%soil_tempk  (k,ipa)   = csite%soil_tempk  (nzg,ipa)
   end do

   do k=1, nzs-1
     csite%sfcwater_tempk(k,ipa)  = csite%soil_tempk  (nzg,ipa)
   end do

   return
end subroutine final_adjust_ebudget
!==========================================================================================!
!==========================================================================================!
