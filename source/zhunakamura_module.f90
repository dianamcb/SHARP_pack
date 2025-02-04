module zhunakamura_module
!**********************************************************************
!     SHARP Pack module that contains all the subroutines for the
!     Zhu-Nakamura method. The subroutines are listed below.
!
!     - ComputeAdiaE: Searches the seam point between two potential
!                     energy surfaces: E_i and E_n.
!     - ZNParams: Computes the Zhu-Nakamura paramters.
!     - ZNHopping: Performs surface hopping according to the
!                  Zhu-Nakamura method.
!     - ZNJCorrection: corrects the total angular momentum of the
!                      system in case of non-vertical hops.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
  use global_module
  use models_module

  contains

  subroutine ComputeAdiaE(outPES, collE, vp, istate, inext)
!**********************************************************************
!     SHARP Pack subroutine that searches the seam point between
!     twho potential energy surfaces: E_i (current state) and E_n
!     (next state).
!
!     Input:
!     - Hel: diabatic hamiltonian.
!     - dHeldr: Derivative of diabatic hamiltonian.
!     - eva: Eigenvalues (adiabatic potentials or PESs).
!     - psi: Eigenvectors (adiabatic bases).
!     - istate: Current state.
!     - inext: Next state.
!
!     Output:
!     - collE: Collision energy.
!     - adE: Adiabatic energies.
!       + adE(:,1): Array of positions from -12 to 12.
!       + adE(:,istate): Adiabatic energy E_i.
!       + adE(:,inext): Adiabatic energy E_n.
!     - outPES: Array with values needed to compute ZN parameters.
!               It contains the values below.
!       + r0: Position at Min(Delta E). Minimum separation point.
!       + adEi_r0: E_i at r0.
!       + adEn_r0: E_n at r0.
!       + adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!       + ri_E0: Position at E_0 = E_i.
!       + rn_E0: Position at E_0 = E_n.
!       + adEi_rn_E0: E_i at rn_E0.
!       + adEn_ri_E0: E_n at ri_E0.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    implicit none

    real*8              :: hel(nstates,nstates)
    real*8              :: x(np), vp(np)
    real*8              :: dheldr(nstates,nstates,np)
    real*8              :: eva(nstates,2)
    real*8              :: psi(nstates,nstates,2)
    integer, intent(in) :: inext, istate

    real*8              :: deltaE(961)
    real*8              :: adE(961,nstates+1)
    integer             :: k, idx, idxMinDeltaE, idxMinEn

    real*8              :: r0, adEi_r0, adEn_r0, adE0
    real*8              :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8, intent(out) :: outPES(8)
    real*8, intent(out) :: collE

    do k = 1,961
      x(:) = -12+(k-1)*0.025

      call GetHel(x(:),hel(:,:),dheldr(:,:,:))
      call Diag(eva(:,1),psi(:,:,2),hel(:,:))

      adE(k,1) = x(1)
      ! Check: E_n >= E_i at r0 is a must.
      adE(k,:) = eva(:,1) !Adiabatic energies from E_1 to E_n.
    enddo

    deltaE(:) = abs(adE(:,inext)-adE(:,istate)) !We ensure DeltaE > 0

    idxMinDeltaE = MinLoc(deltaE, dim=1)
    r0 = adE(idxMinDeltaE,1) ! Position at Min(DeltaE)
    adEi_r0 = adE(idxMinDeltaE,istate) !E_i at r0
    adEn_r0 = adE(idxMinDeltaE,inext) !E_n at r0

    adE0 = 0.5*(adEi_r0 + adEn_r0) ! E0 = Mean(E_i,E_n) at r0

    idxMinEn = MinLoc(adE(:,inext), dim=1)
    if(idxMinEn .le. idxMinDeltaE)then
      idx = LocIdx(adE(idxMinEn:,istate), adE0)
      ri_E0 = adE(idx,1) ! Position at E_i = E0
      idx = LocIdx(adE(idxMinEn:,inext), adE0)
      rn_E0 = adE(idx,1) ! Position at E_n = E0
    endif
    if(idxMinEn .gt. idxMinDeltaE)then
      idx = LocIdx(adE(:idxMinEn,istate), adE0)
      ri_E0 = adE(idx,1) ! Position at E_i = E0
      idx = LocIdx(adE(:idxMinEn,inext), adE0)
      rn_E0 = adE(idx,1) ! Position at E_n = E0
    endif

    outPES(1:4) = (/ r0, adEi_r0, adEn_r0, adE0 /)
    outPES(5:8) = (/ ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0 /)

    collE = adE(idxMinDeltaE,istate) + 0.5*mp*vp(1)**2

    contains

    integer function LocIdx(arr, val)
      implicit none
      real*8, intent(in) :: arr(:), val
      real*8             :: boundR, boundL
      boundL = 1
      boundR = size(arr(:))
      do while(boundL .ne. boundR)
        locIdx = ceiling((boundL+boundR)/2.)
        if(arr(locIdx) .gt. val)then; boundR = locIdx-1
        else; boundL = locIdx
        endif
      enddo
      if(arr(locIdx) .eq. val)then; locIdx = boundL; endif
    end function LocIdx
  end subroutine ComputeAdiaE


  subroutine ZNParams(outPES, aSqr, bSqr)
!**********************************************************************
!     SHARP Pack subroutine that computes the Zhu-Nakamura paramters.
!
!     Input:
!     - mp: Mass of a particle.
!     - collE: Collision energy.
!     - outPES: Array of variables comming from ComputeAdiaE.
!               It contains the values below (output potential)
!       + r0: Position at Min(Delta E). Minimum separation point.
!       + adEi_r0: E_i at r0.
!       + adEn_r0: E_n at r0.
!       + adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!       + ri_E0: Position at E_i = E0.
!       + rn_E0: Position at E_n = E0.
!       + adEi_rn_E0: E_i at rn_E0.
!       + adEn_ri_E0: E_n at ri_E0.
!
!     Output:
!     - aSqr: a^2. Effective coupling constant.
!     - bSqr: b^2. Effective collision energy.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    implicit none

    real*8, intent(out) :: aSqr, bSqr

    real*8              :: dSqr

    real*8              :: collE
    real*8              :: r0, adEi_r0, adEn_r0, adE0
    real*8              :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8, intent(in)  :: outPES(8)

    r0 = outPES(1)
    adEi_r0 = outPES(2)
    adEn_r0 = outPES(3)
    adE0 = outPES(4)
    ri_E0 = outPES(5)
    rn_E0 = outPES(6)
    adEi_rn_E0 = outPES(7)
    adEn_ri_E0 = outPES(8)

    dSqr = (adEn_ri_E0-adE0)*(adE0-adEi_rn_E0)
    dSqr = dSqr/(adEn_r0-adEi_r0)**2

    aSqr = sqrt(dSqr-1.)*hbar**2
    aSqr = aSqr/(mp*(adEi_r0-adEn_r0)*(ri_E0-rn_E0)**2)
    bSqr = sqrt(dSqr-1.)*(collE-adE0)/(0.5*(adEi_r0-adEn_r0))
  end subroutine ZNParams


  subroutine ZNHopping(outPES,adE,collE,vp,aSqr,bSqr,dSqr,istate,inext)
!**********************************************************************
!     SHARP Pack subroutine that performs surface hopping according
!     to the Zhu-Nakamura method.
!
!     Note: This method was implemented for a single particle.
!           Next updates will include the implementation for multiple
!           particles.
!
!     Input:
!     - outPES: Array of variables comming from ComputeAdiaE.
!               It contains the values below.
!       + r0: Position at Min(Delta E). Minimum separation point.
!       + adEi_r0: E_i at r0.
!       + adEn_r0: E_n at r0.
!       + adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!       + ri_E0: Position at E_i = E0.
!       + rn_E0: Position at E_n = E0.
!       + adEi_rn_E0: E_i at rn_0.
!       + adEn_ri_E0: E_n at ri_0.
!     - aSqr: a^2. Effective coupling constant.
!     - bSqr: b^2. Effective collision energy.
!     - dSqr: d^2. No meaning.
!     - adE: Adiabatic potentials E_i and E_n.
!     - collE: Collision energy at current state (istate).
!     - inext: Index of surface energy to hop.
!     - istate: Current surface energy.
!
!     Output: inext, collE, aSqr
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
! Todo: Change FindLoc to LocIdx.
    implicit none

    real*8,  intent(in)    :: vp(np)
    real*8,  intent(in)    :: dSqr
    integer, intent(in)    :: istate

    real*8                 :: deltaE(961)
    real*8                 :: kappaUp(961,nstates+1)
    real*8                 :: aPar, bPar, g1, deltapsi, f1C, f2C
    real*8                 :: f1, f2, bx, dr, pZN, sigma0ZN, delta0ZN
    real*8                 :: gamma1, gamma2, prob(nstates)
    real*8                 :: sigmaZN, deltaZN
    real*8                 :: tmp, beta
    real*8                 :: acumulator, dice, ranf
    integer                :: idx1, idx2, i
    complex*16             :: phiS, ctmp

    real*8                 :: r0, adEi_r0, adEn_r0, adE0
    real*8                 :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8,  intent(inout) :: outPES(8)
    real*8,  intent(inout) :: adE(961,nstates+1)
    real*8,  intent(inout) :: aSqr, bSqr

    real*8,  intent(out)   :: collE
    integer, intent(out)   :: inext

    dr = abs(adE(2,1)-adE(1,1)) ! Deltar (constant value).

    ! Compute transition probabilities for all states.
    do i=1,nstates
      if(i .eq. istate)then; prob(i) = 0.d0; continue; endif

      call ComputeAdiaE(outPES,collE,vp,istate,i)
      call ZNParams(outPES,aSqr,bSqr)

      r0 = outPES(1)
      adEi_r0 = outPES(2)
      adEn_r0 = outPES(3)
      adE0 = outPES(4)
      ri_E0 = outPES(5)
      rn_E0 = outPES(6)
      adEi_rn_E0 = outPES(7)
      adEn_ri_E0 = outPES(8)

      aPar = sqrt(aSqr)
      bPar = sqrt(abs(bSqr))

      if(0.02 .gt. aSqr)then
        write(nrite_hopp,'("Adiabatic region. Unlikely non-adiabatic transition")')
        write(nrite_hopp,'("Skipping ZN hop from current step")')
        return
      endif
      if((0.02 .le. aSqr) .and. (aSqr .le. 1.0d4))then
        write(nrite_hopp,'("Non-adiabatic region. Likely non-ariabatic transition")')
        write(nrite_hopp,'("Performing ZN hop at current step")')
      endif
      if(aSqr .gt. 1.0d4)then
        write(nrite_hopp,'("Diabatic region. Unlikely non-adiabatic transition")')
        write(nrite_hopp,'("Skipping ZN hop from current step")')
        return
      endif

      ! Laundau-Zener case 1: E >= E_n(r0) ==> p_(ZN)
      if(collE .ge. adEn_r0)then
        prob(i)=exp(sqrt(2./(1.+sqrt(1.+bPar**(-4)*(0.7+0.4*aPar**2)))))
        prob(i)=real(prob(i)**(-pi/(4.*aPar*bPar)))
      endif

      ! Laundau-Zener case 2: E < E_n(r0) ==> p_(12)
      if(collE .lt. adEn_r0)then
        bx = bSqr-0.9553
        gamma1 = 0.9*sqrt(dSqr-1.)
        gamma2 = 7*sqrt(dSqr)/16.

        f1 =    sqrt(sqrt((bSqr+gamma1)**2+gamma2)-(bSqr+gamma1))
        f1 = f1+sqrt(sqrt((bSqr-gamma1)**2+gamma2)-(bSqr+gamma1))

        f2 =    sqrt(sqrt((bSqr+gamma1)**2+gamma2)+(bSqr+gamma1))
        f2 = f2+sqrt(sqrt((bSqr-gamma1)**2+gamma2)+(bSqr+gamma1))

        f1C = f1*(0.45*sqrt(dSqr))/(1.+1.5*exp(2.2*bx*abs(bx)**0.57))
        f2C = f2*(bSqr-0.16*bx/sqrt(1.+bSqr**2))

        ctmp = cmplx(f1C,f2C)
        ctmp = ctmp/(f2**2 + f1**2)
        ctmp = ctmp*sqrt(2.)*pi/(4.*aPar)
        sigma0ZN = real(ctmp)
        delta0ZN = aimag(ctmp)

        kappaUp(:,:) = sqrt(2*mp/hbar**2 * (collE-adE(:,:)))

        idx1 = FindLoc(adE(:,1), r0, dim=1)
        idx2 = FindLoc(adE(:,1), ri_E0, dim=1)
        sigmaZN = Integ(kappaUp(idx1:idx2,istate),dr) + sigma0ZN

        idx1 = FindLoc(adE(:,1), r0, dim=1)
        idx2 = FindLoc(adE(:,1), ri_E0, dim=1)
        deltaZN = -Integ(kappaUp(idx1:idx2,istate),dr)
        idx2 = FindLoc(adE(:,1), rn_E0, dim=1)
        deltaZN = deltaZN + Integ(kappaUp(idx1:idx2,i),dr) + delta0ZN

        g1 = 3*sigmaZN/(pi*deltaZN) * log(1.2+aPar**2) - 1./aPar**2

        deltapsi = 1. + 5*sqrt(aPar) * 10**(-sigmaZN)/(sqrt(aPar)+0.8)
        deltapsi = deltapsi*deltaZN

        ctmp = Cgamma(0.d0,deltapsi/pi)
        phiS = deltapsi/pi*log(deltapsi/pi)-atan2(aimag(ctmp),real(ctmp))
        phiS = phiS-deltapsi/pi-pi/4.0

        tmp = sigmaZN/pi
        beta = 2.0d0*pi*exp(-2.0d0*tmp)*tmp**(2.0d0*tmp)/(tmp*gamma(tmp)**2)
        pZN = 1.0d0 + beta*exp(2*deltaZN) - g1*sin(sigmaZN)**2
        pZN = 1.0d0/pZN

        prob(i) = real(4*pZN*(1-pZN)*sin(phiS+sigma0ZN)**2)
      endif
    enddo

    call random_number(ranf)
    dice = ranf
    prob(:) = prob(:)/nstates
    acumulator = 0.0
    do i=1,nstates ! Replaced inext with istate
      acumulator=acumulator+prob(i)
      if(acumulator .ge. dice)then
        inext=i ! We jumped (still have to check energy).
        exit
      endif
    enddo

    contains

    ! Complex gamma function
    complex*16 function Cgamma(zreal, zimag)
      real*8, intent(in) :: zreal, zimag
      complex*16         :: z
      complex*16         :: ctmp
      z = cmplx(zreal,zimag)
      ctmp = 1. + 1./(12.*z) + 1./(288.*z**2) - 139./(51840.*z**3)
      ctmp = ctmp - 571./(2488320.*z**4)
      cgamma = z**(z-0.5)*exp(-z)*sqrt(2.*pi) * ctmp
    end function Cgamma

    ! Composite Simpson integral
    ! Note: Arrays must be at least 9 elements long.
    real*8 function Integ(arr, width)
      real*8, intent(in) :: arr(:), width
      real*8             :: tmp
      integer            :: i, n
      n = size(arr(:))
      integ = 49*arr(n-3)+43*arr(n-2)+59*arr(n-1)+17*arr(n)
      do i=5,n-4
        integ = Integ + 48*arr(i)
      enddo
      integ = integ + 17*arr(1)+59*arr(2)+43*arr(3)+49*arr(4)
      integ = integ * width/48.
    end function Integ

    integer function LocIdx(arr, val)
      implicit none
      real*8, intent(in) :: arr(:), val
      real*8             :: boundR, boundL
      boundL = 1
      boundR = size(arr(:))
      do while(boundL .ne. boundR)
        locIdx = ceiling((boundL+boundR)/2.)
        if(arr(locIdx) .gt. val)then; boundR = locIdx-1
        else; boundL = locIdx
        endif
      enddo
      if(arr(locIdx) .eq. val)then; locIdx = boundL; endif
    end function LocIdx
  end subroutine ZNHopping


  subroutine ZNCorrection(inext,istate,collE,adE,outPES,rp,vp)
!**********************************************************************
!     SHARP Pack subroutine to correct velocity and coordinate after
!     successful hopping.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    use modelvar_module, only : nJump, nJumpFail, nFrust, nFrustR

    implicit none

    real*8                :: r0, adEi_r0, adEn_r0, adE0
    real*8                :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8, intent(in)    :: outPES(8)
    real*8, intent(in)    :: collE
    real*8, intent(in)    :: adE(961,nstates+1)
    real*8                :: deltaE(961)
    real*8                :: ctime, dt, avetim
    integer               :: ip, ibd
    integer               :: istate, inext

    real*8                :: adEn_rnext
    integer               :: upstate
    integer               :: idx, idxMinDeltaE, idxMinEp

    real*8, intent(inout) :: rp(np,nb), vp(np,nb)

    r0 = outPES(1)
    adEi_r0 = outPES(2)
    adEn_r0 = outPES(3)
    adE0 = outPES(4)
    ri_E0 = outPES(5)
    rn_E0 = outPES(6)
    adEi_rn_E0 = outPES(7)
    adEn_ri_E0 = outPES(8)

    ! Laundau-Zener case 1: E >= E_n(r0) ==> Classically allowed hop.
    if(collE .ge. adEn_r0)then
      if(ldtl)WRITE(nrite_hopp,'(" Classically accepted hop")')
      if((ctime*dt/41340.d0) > avetim)then
        nJump(istate,inext) = nJump(istate,inext) + 1
      endif

      ! If hopping (classically allowed) to greater FES.
      if(inext .ge. istate)then
        do ip=1,np
          do ibd=1,nb
            ! Rescaling velocity.
            ! Current and next velocities should have same direction.
            if(vp(ip,ibd) .ge. 0)then
              vp(ip,ibd)=sqrt(vp(ip,ibd)**2 - 2*(adEn_r0-adEi_r0)/mp)
            else
              vp(ip,ibd)=-sqrt(vp(ip,ibd)**2 - 2*(adEn_r0-adEi_r0)/mp)
            endif
          enddo
        enddo
      endif

      istate=inext
    endif

    ! Laundau-Zener case 2: E < E_n(r0) ==> Classically forbidden hop.
    if(collE .lt. adEn_r0)then
      if((ctime*dt/41340.d0)>avetim)then
        nfrust_hop = nfrust_hop + 1
        nFrust(istate,inext) = nFrust(istate,inext) + 1
      endif

      ! Never reverse velocity
      if((ctime*dt/41340.d0)>avetim)then
        nfrust_hop2 = 0
        nFrustR(istate,inext) = nFrustR(istate,inext) + 0
      endif

      if(ldtl)WRITE(nrite_hopp,'(" Classically forbidden hop")')
      if((ctime*dt/41340.d0)>avetim)then
        nJumpFail(istate,inext) = nJumpFail(istate,inext) + 1
      endif

      deltaE(:) = abs(adE(:,inext)-adE(:,istate))
      idxMinDeltaE = MinLoc(deltaE, dim=1) !index of r0

      ! If hopping (classically forbidden) to greater FES.
      if(inext .ge. istate)then
        idxMinEp = MinLoc(adE(:,inext), dim=1) !index of Min(Ep)
        do ip=1,np
          do ibd=1,nb
            adEn_rnext = adEi_r0 + 0.5*mp*vp(ip,ibd)**2
            if(idxMinEp .le. idxMinDeltaE)then
              idx = LocIdx(adE(idxMinEp:,inext),adEn_rnext)
            else
              idx = LocIdx(adE(:idxMinEp,inext),adEn_rnext)
            endif
            ! Rescaling velocity and adjusting position.
            rp(ip,ibd) = adE(idx,1)
            vp(ip,ibd) = 0
          enddo
        enddo
      endif

      ! If hopping (classically forbidden) to lower FES.
      if(inext .lt. istate)then
        idxMinEp = MinLoc(adE(:,istate), dim=1) !index of Min(Ep)
        do ip=1,np
          do ibd=1,nb
            adEn_rnext = adEi_r0 - 0.5*mp*vp(ip,ibd)**2
            if(idxMinEp .le. idxMinDeltaE)then
              idx = LocIdx(adE(idxMinEp:,istate),adEn_rnext)
            else
              idx = LocIdx(adE(:idxMinEp,inext),adEn_rnext)
            endif
            ! Rescaling velocity and adjusting position.
            rp(ip,ibd) = adE(idx,1)
            vp(ip,ibd) = -sqrt(2*(adEi_r0-adEn_rnext)/mp)
          enddo
        enddo
      endif

      istate = inext
    endif

    if(ldtl)WRITE(nrite_hopp,*)

    contains

    integer function LocIdx(arr, val)
      implicit none
      real*8, intent(in) :: arr(:), val
      real*8             :: boundR, boundL
      boundL = 1
      boundR = size(arr(:))
      do while(boundL .ne. boundR)
        locIdx = ceiling(0.5*(boundL+boundR))
        if(arr(locIdx) .gt. val)then; boundR = locIdx-1
        else; boundL = locIdx
        endif
      enddo
      if(arr(locIdx) .eq. val)then; locIdx = boundL; endif
    end function LocIdx
  end subroutine ZNCorrection
!**********************************************************************
end module zhunakamura_module
