!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag
      integer*4 i, iStage
      real*8 fqs_fluct(3)

      real*8 aSDE,bq,bSDE,chi,denum,dW1,dW2,dW3,fq,Fs,gkern,
     >   sigF,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2,Z3
      real*8 TwoPi

!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,"(2x,E16.8,i5,16(1x,F13.8))") ppiclf_time,
     >                denum,upmean,vpmean,wpmean
         endif
      endif

      ! Lattenzi is valid only for incompressible flows,
      !   so here we use the compressible correction of Osnes
      !   Equations (9-12) of Osnes paper, with eqn (9) corrected
      !   Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigF = (fq + bq)*Fs

      chi = (1.0+2.50*rphip+4.51*(rphip**2)+4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      fqs_fluct = 0.0d0

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      ! Particle velocity fluctuation
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      tF_inv = (24.0*rphip*chi)/dp * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE = sigF*sqrt(2.0*tF_inv)

      call RANDOM_NUMBER(UnifRnd)

      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))
      Z3 = sqrt(-2.0d0*log(UnifRnd(5)))*cos(TwoPi*UnifRnd(6))

      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2
      dW3 = sqrt(fac)*Z3

      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE*dW1
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE*dW2
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE*dW3


      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE,sigF,                               ! 51-53
     >         fqs_fluct(1:3),                               ! 54-56
     >         Z1,Z2,Z3,dW1,dW2,dW3,                         ! 57-62
     >         ppiclf_np
         endif
         endif
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024 - T.L. Jackson
! Modified 3/6/24 - Balachandar
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag

      integer*4 i, iStage
      real*8 fqs_fluct(3)

      real*8 aSDE,bq,chi,denum,dW1,dW2,fq,Fs,gkern,
     >   sigD,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2
      real*8 TwoPi
      real*8 bSDE_CD, bSDE_CL,CD_frac,CD_prime
      real*8 sigT,sigCT
      real*8 sigmoid_cf, f_CF
      real*8 avec(3)
      real*8 bvec(3)
      real*8 cvec(3)
      real*8 dvec(3)
      real*8 cosrand,sinrand
      real*8 eunit(3)

!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',
     >          ppiclf_time,denum,upmean,vpmean,wpmean,
     >          abs(upmean),abs(vpmean),abs(wpmean),
     >          u2pmean,v2pmean,w2pmean
         endif
      endif

      ! Equations (9-12) of Osnes paper, with eqn (9) corrected
      ! Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigD = (fq + bq)*Fs

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      chi = (1.0 + 2.50*rphip + 4.51*(rphip**2) + 4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      tF_inv = (24.0*rphip*chi/dp) * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE_CD = sigD*sqrt(2.0*tF_inv)  ! Modified 3/6/24

! Fluctuating perpendicular force
! Compare CD' (units of force) against sigma_CD (units of force) 
!    to determine which of 5 bins to use for sigCT
!
! Added 3/6/24 
! Modified 3/14/24 
!
      ! 03/13/2025 - Thierry - if velocity is very small, don't impose fluctuations
      if(vmag > 1.d-8) then
        avec = [vx,vy,vz]/vmag

        CD_prime = ppiclf_rprop(PPICLF_R_FLUCTFX,i)*avec(1) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFY,i)*avec(2) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFZ,i)*avec(3)
        CD_frac  = CD_prime/sigD

      else
        avec     = [1.0, 0.0, 0.0]
        CD_prime = 0.0
        sigD     = 0.0
        CD_frac  = 0.0
      endif

      ! Thierry Daoud - Updated June 2, 2024
      sigmoid_cf = 1.0 / (1.0 + exp(-CD_frac))
      f_CF = 0.39356905*sigmoid_cf + 0.43758848
      sigT  = f_CF*sigD
      bSDE_CL = sigT*sqrt(2.0*tF_inv)

      if (ppiclf_debug==2) then
      if (i<=4) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',i,
     >        CD_prime,CD_frac,sigD,theta,bSDE_CD,bSDE_CL
         endif
      endif
      endif


! Calculate the three orthogonal unit vectors
! The first vector (avec) is vx/vmag, vy/vmag, and vz/vmag
! The second (bvec) is constructued by taking cross-product with eunit
!   Note: if avec is in dir. of e_x=(1,0,0), use e_y=(0,1,0) to get e_z
!       : if avec is in dir. of e_y=(0,1,0), use e_z=(0,0,1) to get e_x
!       : if avec is in dir. of e_z=(0,0,1), use e_x=(1,0,0) to get e_y
! The third (cvec) is cross product of the first two
! written 3/6/24
!
      eunit = [1,0,0]
      if (abs(avec(2))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,1,0]
      elseif (abs(avec(1))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,0,1]
      endif

      bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
      bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
      bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
      denum   = max(1.d-8,sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
      bvec    = bvec / denum

      cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
      cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
      cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
      denum   = max(1.d-8,sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
      cvec    = cvec / denum

      call RANDOM_NUMBER(UnifRnd)

      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))

      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2

! Random mixture of bvec and cvec - make sure the new one is a unit vector
! Added 3/6/24
!
      cosrand = cos(TwoPi*UnifRnd(5))
      sinrand = sin(TwoPi*UnifRnd(5)) 
      dvec(1) = bvec(1)*cosrand + cvec(1)*sinrand
      dvec(2) = bvec(2)*cosrand + cvec(2)*sinrand
      dvec(3) = bvec(3)*cosrand + cvec(3)*sinrand
      denum   = max(1.d-8,sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2))
      dvec    = dvec/denum

      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE_CD*dW1*avec(1) + bSDE_CL*dW2*dvec(1)
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE_CD*dW1*avec(2) + bSDE_CL*dW2*dvec(2)
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE_CD*dW1*avec(3) + bSDE_CL*dW2*dvec(3)


      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE_CD,bSDE_CL,                         ! 51-53
     >         sigD,sigT,                                    ! 54-55
     >         CD_prime,CD_frac,sigmoid_cf,f_CF,             ! 56-59
     >         eunit,avec,bvec,                              ! 60-68
     >         cvec,dvec,rpi,                                ! 69-75
     >         fqs_fluct(1:3),                               ! 76-78
     >         Z1,Z2,dW1,dW2,                                ! 79-82
     >         ppiclf_np,                                    ! 83
     >         ppiclf_y(PPICLF_JX:PPICLF_JZ,i)               ! 84-86
         endif
         endif
      endif


      return
      end
