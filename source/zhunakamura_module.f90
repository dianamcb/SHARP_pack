module zhunakamura_module
!**********************************************************************
!     SHARP Pack module that contains all the subroutines for the
!     Zhu-Nakamura method. The subroutines are listed below.
!
!     - ComputeAdiaE: Computes all the potential energy surfaces.
!     - ComputeSeamLine: Searches the seam point between two potential
!                        energy surfaces: E_i and E_n.
!     - ZNParams: Computes the Zhu-Nakamura paramters.
!     - ZNHopping: Performs surface hopping according to the
!                  Zhu-Nakamura method.
!     - ZNJCorrection: corrects the total angular momentum of the
!                      system in case of non-vertical hops.
!
!     Note:
!     - Currently, the method is implemented based on the
!       Landau-Zener curves type. Future updates aim to include the
!       non-adiabatic tunneling curves type.
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

  subroutine ComputeAdiaEAndCoupling(adE,d_ab)
!**********************************************************************
!     SHARP Pack subroutine that computes the potential energies
!     (PESs) as functions of coordinates and the nonadiabatic coupling
!     as a function of coordinates.
!
!     Input:
!     - x: Particle coordinate.
!     - Hel: diabatic hamiltonian.
!     - dHeldr: Derivative of diabatic hamiltonian.
!     - eva: Eigenvalues (adiabatic potentials or PESs).
!     - psi: Eigenvectors (adiabatic bases).
!
!     Output:
!     - adE: Adiabatic energies.
!       + adE(:,1): Array of positions from -12 to 12.
!       + adE(:,istate+1): Adiabatic energy E_i.
!       + adE(:,inext+1): Adiabatic energy E_n.
!     - d_ab: Nonadiabatic coupling vector.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    implicit none

    real*8              :: rc(np)
    real*8              :: hel(nstates,nstates)
    real*8              :: dheldr(nstates,nstates,np)
    real*8              :: eva(nstates,2)
    real*8              :: psi(nstates,nstates,2)

    integer             :: i,j,k,l,m,ip
    real*8              :: xlb, xub ! Boudaries of coordinate space
    real*8              :: dx

    real*8, intent(out) :: adE(mesh,nstates+1)
    real*8, intent(out) :: d_ab(mesh,nstates,nstates)

    xlb = -15; xub = 15
    dx = (xub-xlb)/mesh
    do k = 1,int(mesh)
      rc(:) = xlb + (k-1)*dx

      call GetHel(rc(:),hel(:,:),dheldr(:,:,:))
      call Diag(eva(:,1),psi(:,:,2),hel(:,:))

      if(k == 1)psi(:,:,1)=psi(:,:,2)

      do i = 1, nstates
        psi(:,i,2) = psi(:,i,2)*dot_product(psi(:,i,1),psi(:,i,2))/ &
                      abs(dot_product(psi(:,i,1),psi(:,i,2)))
      end do

      ! Compute the adiabatic coupling vector at each position
      d_ab(k,:,:)=0.d0
      do m=1,nstates
        do l=1,nstates
          if (m.ne.l) then
            do i=1,nstates
              do j=1,nstates
                do ip=1,np
                  d_ab(k,m,l) = d_ab(k,m,l)+psi(i,m,1)* &
                  dheldr(i,j,ip)*psi(j,l,1)/(eva(l,1)-eva(m,1))
                end do
              end do
            end do
          end if
        end do
      end do

      psi(:,:,1)=psi(:,:,2)

      adE(k,1) = rc(1)
      adE(k,2:) = eva(:,1) !Adiabatic energies from E_1 to E_n.
    end do
  end subroutine ComputeAdiaEAndCoupling

  subroutine ComputeSeamLine(adE,rc,vc,d_ab,inext, istate, &
                             collE,outPES,lb,ub,tranType)
!**********************************************************************
!     SHARP Pack subroutine that searches the seam point between
!     twho potential energy surfaces: E_i (current state) and E_n
!     (next state).
!
!     Input:
!     - istate: Current state.
!     - inext: Next state.
!     - adE: Adiabatic energies.
!       + adE(:,1): Array of positions from -12 to 12.
!       + adE(:,istate): Adiabatic energy E_i.
!       + adE(:,inext): Adiabatic energy E_n.
!     - xp: Particle coordinate.
!     - vp: Particle velocity.
!
!     Output:
!     - outPES: Array with values needed to compute ZN parameters.
!               It contains the values below.
!       + If LZ type:
!         * r0: Position at Min(Delta E). Minimum separation point.
!         * adEi_r0: E_i at r0.
!         * adEn_r0: E_n at r0.
!         * adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!         * ri_E0: Position at E_0 = E_i.
!         * rn_E0: Position at E_0 = E_n.
!         * adEi_rn_E0: E_i at rn_E0.
!         * adEn_ri_E0: E_n at ri_E0.
!       + If NT type:
!         * adEb: Minimum of E_n
!         * adEt: Minimum of E_i
!         * rb: Position at E_b
!         * rt: Position at E_t
!         * adEn_rbrt: E_n[(rb+rt)/2]
!         * adEi_rbrt: E_i[(rb+rt)/2]
!         * t1l: T_1^l. Left position at E_i = collE
!         * t2r: T_2^r. Right position at E_i = collE
!         * d2Endr2_rb: Second derivative of E_n at rb
!         * d2Eidr2_rt: Second derivative of E_i at rt
!     - lb: Lower bound of adiabatic energy (where LZ method applies).
!     - ub: Upper bound of adiabatic energy (where LZ method applies).
!     - tranType: Transition type. Landau-Zener (LZ) or
!                 nonadiabatic tunneling (NT).
!     - collE: Collision energy or translational energy along the
!              transition direction.
!
!     Note: We are assuming potential energies of the Landau-Zener
!           type: Both curves have same slope around crossing point.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
! Todo: Come with a way to define the boundaries lb and ub according
!       to the aCollE criteria.
    implicit none

    real*8,      intent(in)    :: rc(np), vc(np)
    real*8,      intent(in)    :: d_ab(mesh,nstates,nstates)
    real*8,      intent(in)    :: adE(mesh,nstates+1)
    integer,     intent(in)    :: inext

    real*8                     :: deltaE(mesh),dEn,dEi,dEn_k,dEn_l,dr
    real*8                     :: aCollE, bCollE, cCollE, gammaCollE
    real*8                     :: tmp
    integer                    :: ub, lb
    integer                    :: idx, idxMinDeltaE, idxMinEn, k

    real*8                     :: adEb,adEt,rb,rt,adEn_rbrt,adEi_rbrt
    real*8                     :: d2Endr2_rb, d2Eidr2_rt, t1l, t2r
    real*8                     :: r0,adEi_r0,adEn_r0,adE0
    real*8                     :: ri_E0, rn_E0,adEi_rn_E0,adEn_ri_E0
    real*8,      intent(out)   :: outPES(10)
    character*2, intent(out)   :: tranType

    real*8,      intent(inout) :: collE
    integer,     intent(inout) :: istate

    outPES(:) = 0.0d0
    tranType = 'NT'
    dr = adE(2,1)-adE(1,1) ! Constant value.
    lb = 1; ub = mesh

    ! Locate the index of particle position and its energy (collE).
    idx = LocIdx(adE(:,1),rc(1), lb,ub)

    ! Identify transition type.
    ! Curves with same slope are Landau-Zener type; curves with oposite
    ! slope are nonadiabatic tunneling type.
    dEn = adE(idx,inext+1)-adE(idx-1,inext+1)
    dEi = adE(idx,istate+1)-adE(idx-1,istate+1)
    if(dEn*dEi .gt. 0)tranType = 'LZ' ! Landau-Zener type.
    if(dEn*dEi .le. 0)tranType = 'NT' ! Nonadiabatic tunneling type.

    ! Locate the minimum energy surface gap in the local region
    deltaE(:) = abs(adE(:,inext+1)-adE(:,istate+1)) !We ensure DeltaE>0
    idxMinDeltaE = MinLoc(deltaE(:), lb,ub)

    ! Define an arbitrary local region where hop may be attempted
    aCollE = d_ab(idx,istate,inext)**2/(2.d0*mp)
    bCollE = vc(1)*d_ab(idx,istate,inext)
    cCollE = adE(idx,inext+1)-adE(idx,istate+1)
    if(aCollE .lt. 1.95d-4) istate = inext
    if(bCollE**2-4.d0*aCollE*cCollE .le. 0.d0) istate = inext

    ! Calculate the collision energy
    if(istate .ne. inext)then
      aCollE = d_ab(idx,istate,inext)**2/(2.d0*mp)
      bCollE = vc(1)*d_ab(idx,istate,inext)
      cCollE = adE(idx,inext+1)-adE(idx,istate+1)
      if(bCollE .ge. 0.d0)then
        gammaCollE = -bCollE+sqrt(bCollE**2-4.d0*aCollE*cCollE)
        gammaCollE = gammaCollE/(2.d0*aCollE)
      else
        gammaCollE = -bCollE-sqrt(bCollE**2-4.d0*aCollE*cCollE)
        gammaCollE = gammaCollE/(2.d0*aCollE)
      endif
      collE = gammaCollE*d_ab(idx,istate,inext)/mp
      collE = vc(1)-collE
      collE = 0.5*mp*collE**2
      collE = adE(idx,istate+1) + 2.0*collE
    end if

    ! Landau-Zener type
    if(tranType .eq. 'LZ')then
      r0 = adE(idxMinDeltaE,1) ! Position at Min(DeltaE)
      adEi_r0 = adE(idxMinDeltaE,istate+1) !E_i at r0
      adEn_r0 = adE(idxMinDeltaE,inext+1) !E_n at r0
      adE0 = 0.5*(adEi_r0 + adEn_r0) ! E0 = Mean(E_i,E_n) at r0

      idx = LocIdx(adE(:,istate+1), adE0, lb,ub)
      ri_E0 = adE(idx,1) ! Position at E_i = E0
      adEn_ri_E0 = adE(idx,inext+1)

      idx = LocIdx(adE(:,inext+1), adE0, lb,ub)
      rn_E0 = adE(idx,1) ! Position at E_n = E0
      adEi_rn_E0 = adE(idx,istate+1)

      outPES(1:4)  = (/ r0, adEi_r0, adEn_r0, adE0 /)
      outPES(5:8)  = (/ ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0 /)
      outPES(9:10) = (/ 0.0d00, 0.0d0 /)
    endif

    ! Nonadiabatic tunneling type.
    if(tranType .eq. 'NT')then
      idx = MinLoc(adE(:,inext+1), lb,ub)
      rb = adE(idx,1)
      adEb = adE(idx, inext+1)
      d2Endr2_rb = adE(idx-1,inext+1)+adE(idx+1,inext+1)
      d2Endr2_rb = (d2Endr2_rb-2.*adE(idx,inext+1))/(1.*dr)**2

      idx = MaxLoc(adE(:,istate+1), lb,ub)
      rt = adE(idx,1)
      adEt = adE(idx, istate+1)
      t1l = LocIdx(adE(:,istate+1), collE, lb,idx)
      t2r = LocIdx(adE(:,istate+1), collE, idx,ub)
      d2Eidr2_rt = adE(idx-1,istate+1)+adE(idx+1,istate+1)
      d2Eidr2_rt = (d2Eidr2_rt-2.*adE(idx,istate+1))/(1.*dr)**2

      idx = LocIdx(adE(:,1), 0.5*(rb+rt), lb,ub)
      adEi_rbrt = adE(idx, istate+1)
      adEn_rbrt = adE(idx, inext+1)

      outPES(1:6)  = (/ rb, adEb, rt, adEt, adEi_rbrt, adEn_rbrt /)
      outPES(7:10) = (/ d2Eidr2_rt, d2Endr2_rb, t1l, t2r /)
    endif

    contains

    ! Locate index of arr at max(arr(lb:ub)).
    integer function MaxLoc(arr,lb,ub)
      implicit none
      real*8,  intent(in)  :: arr(:)
      real*8               :: tol,val
      integer              :: i, lb, ub
      tol = 1.0d5
      val = nint(arr(lb)*tol)/tol
      maxLoc = lb
      do i=lb,ub
        if(arr(i) .ge. val)then
          maxLoc = i; val = nint(arr(i)*tol)/tol
        end if
      end do
    end function MaxLoc
    ! Locate index of arr at min(arr(lb:ub)).
    integer function MinLoc(arr,lb,ub)
      implicit none
      real*8,  intent(in)  :: arr(:)
      real*8               :: tol,val
      integer              :: i, lb, ub
      tol = 1.0d5
      val = nint(arr(lb)*tol)/tol
      minLoc = lb
      do i=lb,ub
        if(arr(i) .le. val)then
          minLoc = i; val = nint(arr(i)*tol)/tol
        end if
      end do
    end function MinLoc
    ! Locate index of arr where arr(idx) == val.
    integer function LocIdx(arr, val, lb, ub)
      implicit none
      real*8,  intent(in)  :: arr(:), val
      integer              :: i,lb,ub
      locIdx = lb
      do i=lb+1,ub-1
        if((arr(i-1) .lt. val) .and. (val .lt. arr(i+1)))then
          locIdx = i; exit
        end if
        if((arr(i-1) .gt. val) .and. (val .gt. arr(i+1)))then
          locIdx = i; exit
        end if
      end do
    end function LocIdx
  end subroutine ComputeSeamLine

  subroutine ZNParams(outPES,collE,tranType, aSqr,bSqr,dSqr)
!**********************************************************************
!     SHARP Pack subroutine that computes the Zhu-Nakamura paramters.
!
!     Input:
!     - mp: Mass of a particle.
!     - collE: Collision energy or translational energy along the
!              transition direction.
!     - outPES: Array of variables comming from ComputeSeamLine.
!               It contains the values below (output potential)
!       + If LZ type:
!         * r0: Position at Min(Delta E). Minimum separation point.
!         * adEi_r0: E_i at r0.
!         * adEn_r0: E_n at r0.
!         * adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!         * ri_E0: Position at E_0 = E_i.
!         * rn_E0: Position at E_0 = E_n.
!         * adEi_rn_E0: E_i at rn_E0.
!         * adEn_ri_E0: E_n at ri_E0.
!         * ri_collE: Position at E_i = collE
!         * rn_collE: Position at E_n = collE
!       + If NT type:
!         * adEb: Minimum of E_n
!         * adEt: Minimum of E_i
!         * rb: Position at E_b
!         * rt: Position at E_t
!         * adEi_rbrt: E_i[(rb+rt)/2]
!         * adEn_rbrt: E_n[(rb+rt)/2]
!         * t1l: T_1^l. Left position at E_i = collE
!         * t2r: T_2^r. Right position at E_i = collE
!         * d2Eidr2_rt: Second derivative of E_i at rt
!         * d2Endr2_rb: Second derivative of E_n at rb
!     - tranType: Transition type. Landau-Zener (LZ) or
!                  nonadiabatic tunneling (NT).
!
!     Output:
!     - aSqr: a^2. Effective coupling constant.
!     - bSqr: b^2. Effective collision energy.
!     - dSqr: d^2. No meaning. It plays the role of gammaSqr for
!                  the nonadiabatic tunneling transition type.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    implicit none

    real*8              :: rb, adEb, rt, adEt, adEi_rbrt, adEn_rbrt
    real*8              :: d2Endr2_rb, d2Eidr2_rt, t1l, t2r
    real*8              :: r0, adEi_r0, adEn_r0, adE0
    real*8              :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8              :: ri_collE, rn_collE
    real*8, intent(in)  :: outPES(10)
    real*8, intent(in)  :: collE
    character*2         :: tranType

    real*8              :: d2Edr2_rb, d2Edr2_rt

    real*8, intent(out) :: aSqr, bSqr, dSqr

    if(tranType .eq. 'LZ')then
      r0 = outPES(1)
      adEi_r0 = outPES(2)
      adEn_r0 = outPES(3)
      adE0 = outPES(4)
      ri_E0 = outPES(5)
      rn_E0 = outPES(6)
      adEi_rn_E0 = outPES(7)
      adEn_ri_E0 = outPES(8)
      ri_collE = outPES(9)
      rn_collE = outPES(10)

      dSqr = abs((adEn_ri_E0-adE0)*(adE0-adEi_rn_E0))
      dSqr = dSqr/(adEn_r0-adEi_r0)**2

      aSqr = sqrt(dSqr-1.)*hbar**2
      aSqr = abs(aSqr/(mp*(adEn_r0-adEi_r0)*(rn_E0-ri_E0)**2))
      bSqr = sqrt(dSqr-1.)*(collE-adE0)/(0.5*(adEn_r0-adEi_r0))
    endif

    if(tranType .eq. 'NT')then
      rb = outPES(1)
      adEb = outPES(2)
      rt = outPES(3)
      adEt = outPES(4)
      adEi_rbrt = outPES(5)
      adEn_rbrt = outPES(6)
      d2Eidr2_rt = outPES(7)
      d2Endr2_rb = outPES(8)
      t1l = outPES(9)
      t2r = outPES(10)

      dSqr = ((adEb-adEt)/(adEn_rbrt-adEi_rbrt))**2 ! gamma^2

      aSqr = ((1.-dSqr)*hbar**2/(mp*abs(adEb-adEt)*(rb-rt)**2))
      bSqr = (collE-0.5*(adEb+adEt))/(0.5*(adEb-adEt))

      if(rb .eq. rt)then
        aSqr = hbar**2/(4.*mp*(adEb-adEt))
        aSqr = (aSqr * (d2Endr2_rb-d2Eidr2_rt))
      endif
    endif

  end subroutine ZNParams


  subroutine ZNHopping(rc,vc,adE,d_ab ,istate, &
                       inext,collE,outPES,lb,ub,tranType)
!**********************************************************************
!     SHARP Pack subroutine that performs surface hopping according
!     to the Zhu-Nakamura method.
!
!     Important note:
!     - This method was implemented for a single particle. Next
!       updates will include the implementation for multiple particles.
!
!     Input:
!     - xp: Particle coordinate.
!     - vp: Particle velocity.
!     - adE: Adiabatic potentials E_i and E_n.
!     - istate: Current surface energy.
!     - d_ab: Nonadiabatic coupling vector.
!
!     Output:
!     - inext: Index of surface energy to hop.
!     - outPES: Array with values needed to compute ZN parameters.
!               It contains the values below.
!       + If LZ type:
!         * r0: Position at Min(Delta E). Minimum separation point.
!         * adEi_r0: E_i at r0.
!         * adEn_r0: E_n at r0.
!         * adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!         * ri_E0: Position at E_0 = E_i.
!         * rn_E0: Position at E_0 = E_n.
!         * adEi_rn_E0: E_i at rn_E0.
!         * adEn_ri_E0: E_n at ri_E0.
!         * ri_collE: Position at E_i = collE
!         * rn_collE: Position at E_n = collE
!       + If NT type:
!         * adEb: Minimum of E_n
!         * adEt: Minimum of E_i
!         * rb: Position at E_b
!         * rt: Position at E_t
!         * adEi_rbrt: E_i[(rb+rt)/2]
!         * adEn_rbrt: E_n[(rb+rt)/2]
!         * t1l: T_1^l. Left position at E_i = collE
!         * t2r: T_2^r. Right position at E_i = collE
!         * d2Eidr2_rt: Second derivative of E_i at rt
!         * d2Endr2_rb: Second derivative of E_n at rb
!     - lb: Lower bound of adiabatic energy (where LZ method applies).
!     - ub: Upper bound of adiabatic energy (where LZ method applies).
!     - tranType: Transition type. Landau-Zener (LZ) or
!                 nonadiabatic tunneling (NT).
!     - collE: Collision energy or translational energy along the
!              transition direction.
!
!     Notes:
!     - ComputeSeamLine computes collE, not ZNHopping. However,
!       ZNHopping retunrs collE for ZNCorection.
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    implicit none

    real*8                     :: rc(np)
    real*8,      intent(in)    :: vc(np),adE(mesh,nstates+1)
    real*8,      intent(in)    :: d_ab(nstates,nstates)

    real*8                     :: aSqr, bSqr, dSqr
    real*8                     :: aPar, deltapsi, bx, dr
    real*8                     :: deltaZN, sigma0ZN, delta0ZN
    real*8                     :: sigmaZN, g1, pZN, phiS
    real*8                     :: tmp, tmp2, beta
    real*8                     :: prob(nstates)
    real*8                     :: acumulator, dice, ranf
    real*8                     :: h4, h3, h2, wSqr
    real*8                     :: tmpArr(100,2), sigmac
    integer                    :: idx1, idx2, i, j, tmpu
    complex*16                 :: bPar, dPar, ctmp
    complex*16                 :: f1, f2, f1C, f2C, gamma1, gamma2
    complex*16                 :: kappaUp(mesh,nstates+1)

    real*8                     :: rb, adEb, rt, adEt, adEi_rbrt
    real*8                     :: t1l, t2r
    real*8                     :: adEn_rbrt, d2Endr2_rb, d2Eidr2_rt
    real*8                     :: r0, adEi_r0, adEn_r0, adE0
    real*8                     :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8                     :: ri_collE, rn_collE
    real*8,      intent(out)   :: outPES(10)
    real*8,      intent(out)   :: collE
    integer,     intent(out)   :: lb, ub !lower bound, upper bound.
    integer,     intent(out)   :: inext
    character*2, intent(out)   :: tranType

    integer,     intent(inout) :: istate

    prob(:) = 0.0d0
    dr = abs(adE(2,1)-adE(1,1)) ! Deltar (constant value).

    ! Compute transition probabilities for all states.
    do i=1,nstates
      ! By default, hopping prob. from state i to state i is zero.
      if(istate .eq. i) cycle

      ! For Nonadiabatic tunneling surfaces, hops to lower
      ! surfaces are not allowed. Hence prob(i) = 0.
      if(tranType .eq. 'NT')then
        if(istate .gt. i) cycle
      end if

      call ComputeSeamLine(adE,rc,vc,d_ab,i, istate, &
                           collE,outPES,lb,ub,tranType)

      if(istate .eq. i) cycle ! Far from the crossing point.

      call ZNParams(outPES,collE,tranType, aSqr,bSqr,dSqr)

      if(tranType .eq. 'LZ')then
        r0 = outPES(1)
        adEi_r0 = outPES(2)
        adEn_r0 = outPES(3)
        adE0 = outPES(4)
        ri_E0 = outPES(5)
        rn_E0 = outPES(6)
        adEi_rn_E0 = outPES(7)
        adEn_ri_E0 = outPES(8)
        ri_collE = outPES(9)
        rn_collE = outPES(10)
      end if
      if(tranType .eq. 'NT')then
        rb = outPES(1)
        adEb = outPES(2)
        rt = outPES(3)
        adEt = outPES(4)
        adEi_rbrt = outPES(5)
        adEn_rbrt = outPES(6)
        d2Eidr2_rt = outPES(7)
        d2Endr2_rb = outPES(8)
        t1l = outPES(9)
        t2r = outPES(10)
      end if

      aPar = sqrt(aSqr)
      bPar = sqrt(cmplx(bSqr,0.0d0))
      dPar = sqrt(cmplx(dSqr,0.0d0))

      kappaUp(:,1)  = adE(:,1)
      kappaUp(:,2:) = 2.*mp/hbar**2*cmplx(collE-adE(:,2:),0.0d0)
      kappaUp(:,2:) = sqrt(kappaUp(:,2:))

      if(tranType .eq. 'LZ')then
        ! Laundau-Zener case 1: E >= E_n(r0) ==> p_(ZN)
        if(collE .ge. adEn_r0)then
          prob(i)=-sqrt(2./(1.+sqrt(1.+(0.7+0.4*aSqr)/bSqr**2)))
          prob(i)=prob(i)*pi/(4.*aPar*abs(bPar))
          prob(i)=exp(prob(i))
        end if

        ! Laundau-Zener case 2: E < E_n(r0) ==> p_(12)
        if(collE .lt. adEn_r0)then
          bx = bSqr-0.9553
          gamma1 = 0.9*sqrt(cmplx(dSqr-1.0d0,0.0d0))
          gamma2 = 7.0*dPar/16.0

          f1=   sqrt(sqrt((bSqr+gamma1)**2+gamma2)-(bSqr+gamma1))
          f1=f1+sqrt(sqrt((bSqr-gamma1)**2+gamma2)-(bSqr-gamma1))

          f2=   sqrt(sqrt((bSqr+gamma1)**2+gamma2)+(bSqr+gamma1))
          f2=f2+sqrt(sqrt((bSqr-gamma1)**2+gamma2)+(bSqr-gamma1))

          f1C = f1*(0.45*dPar)/(1.+1.5*exp(2.2*bx*abs(bx)**0.57))
          f2C = f2*(bSqr-0.16*bx/sqrt(1.+bSqr**2))

          ctmp = f1C + cmplx(0.0d0,1.0d0)*f2C
          ctmp = ctmp/(f2**2 + f1**2)
          ctmp = ctmp*sqrt(2.0d0)*pi/(4.*aPar)
          sigma0ZN = real(ctmp)
          delta0ZN = aimag(ctmp)

          if(collE .gt. adEi_r0)then
            idx1 = LocIdx(adE(:,1), ri_collE, lb,ub)
            idx2 = LocIdx(adE(:,1), r0, lb,ub)
            sigmaZN = Integ(abs(kappaUp(idx1:idx2,istate+1)),dr)
            sigmaZN = sigmaZN + sigma0ZN

            idx1 = LocIdx(adE(:,1), r0, lb,ub)
            idx2 = LocIdx(adE(:,1), rn_collE, lb,ub)
            deltaZN = Integ(abs(kappaUp(idx1:idx2,i+1)),dr)
            deltaZN = deltaZN + delta0ZN
          end if

          if(collE .le. adEi_r0)then
            sigmaZN = sigma0ZN

            idx1 = LocIdx(adE(:,1), r0, lb,ub)
            idx2 = LocIdx(adE(:,1), ri_collE, lb,ub)
            deltaZN = -Integ(abs(kappaUp(idx1:idx2,istate+1)),dr)
            idx2 = LocIdx(adE(:,1), rn_collE, lb,ub)
            deltaZN = deltaZN+Integ(abs(kappaUp(idx1:idx2,i+1)),dr)
            deltaZN = deltaZN + delta0ZN
          end if

          g1=3.0*sigmaZN/(pi*deltaZN)*log(1.2+aSqr)-1.0/aSqr

          deltapsi=1.+5.*sqrt(aPar)*10**(-sigmaZN)/(sqrt(aPar)+0.8)
          deltapsi=deltapsi * deltaZN

          ctmp = Cgamma(0.d0,deltapsi/pi)
          phiS = atan2(aimag(ctmp),real(ctmp))
          phiS = deltapsi/pi * log(abs(deltapsi/pi)) - phiS
          phiS = phiS - deltapsi/pi - pi/4.0

          tmp = sigmaZN/pi
          beta = 2.0d0*pi*exp(-2.0d0*tmp)
          beta = beta*tmp**(2.0d0*tmp)/(tmp*gamma(tmp)**2)
          pZN=1.0d0+beta*exp(2.0d0*deltaZN)-g1*sin(sigmaZN)**2
          pZN = 1.0d0/pZN

          prob(i) = real(4*pZN*(1-pZN)*sin(phiS+sigma0ZN)**2)
        end if
      end if

      if(tranType .eq. 'NT')then
        ! Nonadiabatic tunneling case 1: E >= Eb ==> p_(ZN)
        if(collE .ge. adEb)then
          prob(i)=sqrt(1.-(0.72-0.62*aPar**(1.43))/bSqr**2)
          prob(i)=-sqrt(2./(1.+prob(i)))
          prob(i)=prob(i)*pi/(4.*aPar*abs(bPar))
          prob(i)=exp(prob(i))

        ! Nonadiabatic tunneling case 2: Et <= E < Eb ==> p_(12)
        elseif(collE .le. adEt)then
          h4 = 0.61*sqrt(2.+bSqr)
          h3 = (aPar-3.*bSqr)*sqrt(1.23+bSqr)/(aPar+3.)
          h2 = 1. + 0.38/aSqr * (1.+bSqr)**(1.2-0.4*bSqr)

          do j=1,99
            tmp = aSqr**(1./3.)
            tmp2 = aPar**(1./3.)
            tmpArr(j,1) = real(j-1,kind=8)/99.0d0
            tmpArr(j,2) = tmpArr(j,1)**3/(3.*(1.-tmpArr(j,1)**3))
            tmpArr(j,2) = tmpArr(j,2)-bSqr*tmpArr(j,1)/(tmp*(1.-tmpArr(j,1)))
            tmpArr(j,2) = tmpArr(j,2)-h3*(1.-tmpArr(j,1))/(tmp*h4*(1.-tmpArr(j,1))+tmp2*tmpArr(j,1))
            tmpArr(j,2) = cos(tmpArr(j,2))/(1.-tmpArr(j,1))**2
          end do

          tmp = tmpArr(2,1)-tmpArr(1,1)
          wSqr = Integ(tmpArr(:99,2),tmp)
          wSqr = wSqr*wSqr

          prob(i) = wSqr/(1.+wSqr)

        ! Nonadiabatic tunneling case 3: E < Et ==> p_(12)
        else
          sigmaZN=pi/(16.*aPar*abs(bPar))
          sigmaZN=sigmaZN*sqrt(6.+10.*sqrt(1.-1./bSqr**2))
          sigmaZN=sigmaZN/(1.+sqrt(1.-1./bSqr**2))

          idx1 = LocIdx(adE(:,1),t1l, lb,ub)
          idx2 = LocIdx(adE(:,1),t2r, lb,ub)
          deltaZN = Integ(abs(kappaUp(idx1:idx2,istate+1)),dr)

          sigmac=sigmaZN*(1.-0.32*10.**(-2./aSqr)*exp(-deltaZN))

          tmp = sigmac/pi
          beta = 2.0d0*pi*exp(-2.0d0*tmp)*tmp**(2.0d0*tmp)
          beta = beta/(tmp*gamma(tmp)**2)

          tmp = beta*exp(-2.*deltaZN)
          prob(i) = tmp/(tmp+(1.+(0.5*aPar/(aPar+1.))*tmp)**2)
        end if
      end if
    end do

    ! Attempting hop
    call random_number(ranf)
    dice = ranf
    prob(:) = prob(:)/nstates
    istate = inext
    acumulator = 0.0d0
    do i=1,nstates
      acumulator=acumulator+prob(i)
      if(acumulator .ge. dice)then
        inext=i ! We jumped (still have to check energy).
        if(ldtl)then
          write(nrite_hopp,*)
          write(nrite_hopp,'(" Trying to hop from ",I3," to ",I3)') &
                             istate,inext
        endif
        exit
      end if
    end do

    contains

    ! Complex gamma function
    complex*16 function Cgamma(zreal, zimag)
      real*8, intent(in) :: zreal, zimag
      complex*16         :: z
      complex*16         :: ctmp
      z = cmplx(zreal,zimag)
      ctmp=1.+1./(12.*z)+1./(288.*z**2)-139./(51840.*z**3)
      ctmp=ctmp - 571./(2488320.*z**4)
      cgamma = z**(z-0.5)*exp(-z)*sqrt(2.*pi) * ctmp
    end function Cgamma

    ! Composite Simpson integral
    ! Note: Arrays must be at least 9 elements long.
    real*8 function Integ(arr, width)
      real*8, intent(in) :: arr(:), width
      integer            :: i, n
      n = size(arr(:))
      integ = 49.*arr(n-3)+43.*arr(n-2)+59.*arr(n-1)+17.*arr(n)
      do i=5,n-4
        integ = integ + 48.*arr(i)
      end do
      integ = integ + 17.*arr(1)+59.*arr(2)+43.*arr(3)+49.*arr(4)
      integ = integ * width/48.
    end function Integ

    ! Locate index of arr where arr(idx) == val.
    integer function LocIdx(arr, val, lb,ub)
      implicit none
      real*8,  intent(in)  :: arr(:), val
      integer              :: i,lb,ub
      locIdx = lb
      do i=lb+1,ub-1
        if((arr(i-1) .gt. val) .and. (val .gt. arr(i+1)))then
          locIdx = i; exit
        end if
        if((arr(i+1) .gt. val) .and. (val .gt. arr(i-1)))then
          locIdx = i; exit
        end if
      end do
    end function LocIdx
  end subroutine ZNHopping

  subroutine ZNCorrect(istate,inext,adE,lb,ub,outPES,tranType,d_ab, &
                          collE,rp,vp,rc,vc, ctime)
!**********************************************************************
!     SHARP Pack subroutine that corrects velocity and coordinate
!     after successful hopping.
!
!     Input:
!     - istate: Current surface energy.
!     - inext: Index of surface energy to hop.
!     - rp: Particle coordinate.
!     - vp: Particle velocity.
!     - adE: Adiabatic potentials E_i and E_n.
!     - collE: Collision energy or translational energy along the
!              transition direction.
!     - lb: Lower bound of adiabatic energy (where LZ method applies).
!     - ub: Upper bound of adiabatic energy (where LZ method applies).
!     - outPES: Array with values needed to compute ZN parameters.
!               It contains the values below.
!       + If LZ type:
!         * r0: Position at Min(Delta E). Minimum separation point.
!         * adEi_r0: E_i at r0.
!         * adEn_r0: E_n at r0.
!         * adE0: E_0 = Mean(E_i, E_n) at r0. Crossing energy.
!         * ri_E0: Position at E_0 = E_i.
!         * rn_E0: Position at E_0 = E_n.
!         * adEi_rn_E0: E_i at rn_E0.
!         * adEn_ri_E0: E_n at ri_E0.
!       + If NT type:
!         * adEb: Minimum of E_n
!         * adEt: Minimum of E_i
!         * rb: Position at E_b
!         * rt: Position at E_t
!         * adEi_rbrt: E_i[(rb+rt)/2]
!         * adEn_rbrt: E_n[(rb+rt)/2]
!         * t1l: T_1^l. Left position at E_i = collE
!         * t2r: T_2^r. Right position at E_i = collE
!         * d2Eidr2_rt: Second derivative of E_i at rt
!         * d2Endr2_rb: Second derivative of E_n at rb
!     - tranType: Transition type. Landau-Zener (LZ) or
!                  nonadiabatic tunneling (NT).
!
!     Output:
!     - rp: Particle coordinate.
!     - vp: Particle velocity.
!
!
!     Authors    - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!     Copyright  - D.M. Castaneda-Bagatella & D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
    use modelvar_module, only : nJump, nJumpFail, nFrust, nFrustR

    implicit none

    real*8                  :: rb, adEb, rt, adEt, adEi_rbrt
    real*8                  :: t1l, t2r
    real*8                  :: adEn_rbrt, d2Endr2_rb, d2Eidr2_rt
    real*8                  :: r0, adEi_r0, adEn_r0, adE0
    real*8                  :: ri_E0, rn_E0, adEi_rn_E0, adEn_ri_E0
    real*8                  :: ri_collE, rn_collE
    real*8,  intent(in)     :: outPES(10)
    real*8,  intent(in)     :: adE(mesh,nstates+1)
    real*8,  intent(in)     :: d_ab(mesh,nstates,nstates)
    integer, intent(in)     :: lb, ub
    character*2, intent(in) :: tranType

    integer                 :: ctime
    real*8                  :: avetim
    integer                 :: ip, ibd
    integer                 :: istate, inext

    real*8                  :: adEn_rnext, adEn, adEi
    real*8                  :: aCollE, bCollE, cCollE, gammaCollE
    integer                 :: upstate
    integer                 :: idx, idxMinDeltaE, idxMinEp

    real*8, intent(inout)   :: rp(np,nb),vp(np,nb),rc(np),vc(np),collE

    avetim = 0.d0

    if(tranType .eq. 'LZ')then
      r0 = outPES(1)
      adEi_r0 = outPES(2)
      adEn_r0 = outPES(3)
      adE0 = outPES(4)
      ri_E0 = outPES(5)
      rn_E0 = outPES(6)
      adEi_rn_E0 = outPES(7)
      adEn_ri_E0 = outPES(8)
      ri_collE = outPES(9)
      rn_collE = outPES(10)

      ! Laundau-Zener case 1: E >= E_n(r0) ==> Classically allowed hop.
      if(collE .ge. adEn_r0)then
        if(ldtl)WRITE(nrite_hopp,'(" LZ: Enough Kinetic energy to hop, &
                                   accepted")')
        if((ctime*dt/41340.d0) > avetim)then
           nJump(istate,inext) = nJump(istate,inext) + 1
        endif

        ! If hopping (classically allowed) to higher FES.
        do ip=1,np
          do ibd=1,nb
            ! Rescaling velocity.
            ! Current and next velocities should have same direction.
            if(vp(ip,ibd) .ge. 0)then
              vp(ip,ibd)=sqrt(vp(ip,ibd)**2 - 2*(adEn_r0-adEi_r0)/mp)
            else
              vp(ip,ibd)=-sqrt(vp(ip,ibd)**2 - 2*(adEn_r0-adEi_r0)/mp)
            end if
          end do
        end do

        istate = inext
      end if

      ! Laundau-Zener case 2: E < E_n(r0) ==> Classically forbidden hop.
      if(collE .lt. adEn_r0)then
        if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop = nfrust_hop + 1
          nFrust(istate,inext) = nFrust(istate,inext) + 1
        end if

        ! Never reverse velocity
        if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = 0
          nFrustR(istate,inext) = nFrustR(istate,inext) + 0
        end if

        if(ldtl)WRITE(nrite_hopp,'(" LZ: Not enough Kinetic energy to hop. &
                                    Still, performing hop. &
                                    Classically forbidden hop.")')
        if((ctime*dt/41340.d0)>avetim)then
          nJumpFail(istate,inext) = nJumpFail(istate,inext) + 1
        end if

        ! If hopping (classically forbidden) to greater FES.
        if(inext .ge. istate)then
          do ip=1,np
            do ibd=1,nb
              adEn_rnext = adEi_r0 + 0.5*mp*vp(ip,ibd)**2
              idx = LocIdx(adE(:,inext+1),adEn_rnext, lb,ub)
              ! Rescaling velocity and adjusting position.
              rp(ip,ibd) = adE(idx,1)
              vp(ip,ibd) = 0.0d0
            end do
          end do
        end if

        ! If hopping (classically forbidden) to lower FES.
        if(inext .lt. istate)then
          idxMinEp = MinLoc(adE(:,istate+1), lb,ub) !index of Min(Ep)
          do ip=1,np
            do ibd=1,nb
              adEn_rnext = adEi_r0 - 0.5*mp*vp(ip,ibd)**2
              idx = LocIdx(adE(:,istate+1),adEn_rnext, lb,ub)
              ! Rescaling velocity and adjusting position.
              rp(ip,ibd) = adE(idx,1)
              vp(ip,ibd) = -sqrt(2*(adEi_r0-adEn_rnext)/mp)
            end do
          end do
        end if

        istate = inext
      end if
    end if

    if(tranType .eq. 'NT')then
      rb = outPES(1)
      adEb = outPES(2)
      rt = outPES(3)
      adEt = outPES(4)
      adEi_rbrt = outPES(5)
      adEn_rbrt = outPES(6)
      d2Eidr2_rt = outPES(7)
      d2Endr2_rb = outPES(8)
      t1l = outPES(9)
      t2r = outPES(10)

      ! Nonadiabatic tunneling case 1: E >= Eb ==> Accepted hop.
      if(collE .ge. adEb)then
      ! if(collE .ge. 0.5*(adEb+adEt))then
        if(ldtl)WRITE(nrite_hopp,'(" NT: Enough Kinetic energy to hop &
                                   (Eb <= E), accepted")')
        if((ctime*dt/41340.d0) > avetim)then
           nJump(istate,inext) = nJump(istate,inext) + 1
        endif

        ! Adjusting velocity of each bead.
        do ip=1,np
          do ibd=1,nb
            ! Rescaling velocity.
            ! Current and next velocities should have same direction.
            idx = LocIdx(adE(:,1),rp(ip,ibd), 1,mesh)
            adEn = adE(idx,inext+1)
            adEi = adE(idx,istate+1)
            if(vp(ip,ibd) .ge. 0)then
              vp(ip,ibd)= sqrt(vp(ip,ibd)**2 - 2*(adEn-adEi)/mp)
            else
              vp(ip,ibd)=-sqrt(vp(ip,ibd)**2 - 2*(adEn-adEi)/mp)
            end if
          end do
        end do

        ! Adjusting position and velocity of centroid.
        rc=0.d0; vc=0.d0
        do ibd=1,nb
          rc(:)=rc(:)+rp(:,ibd)
          vc(:)=vc(:)+vp(:,ibd)
        enddo
        rc(:)=rc(:)/real(nb)
        vc(:)=vc(:)/real(nb)

        ! Setting PES after attempted hop.
        istate = inext

      ! Nonadiabatic tunneling case 2: Et <= E < Eb ==> Rejected hop.
      ! Velocity is not rescaled as particle does not hop, so it
      ! stays at current state (same vp and same xp).
      elseif(collE .ge. adEt)then
        if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop = nfrust_hop + 1
          nFrust(istate,inext) = nFrust(istate,inext) + 1
        end if

        if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = 0
          nFrustR(istate,inext) = nFrustR(istate,inext) + 0
        end if

        if(ldtl)WRITE(nrite_hopp,'(" NT: Not enough kinetic energy &
                                    (Et <= E < Eb), rejected")')
        if((ctime*dt/41340.d0)>avetim)then
          nJumpFail(istate,inext) = nJumpFail(istate,inext) + 1
        end if

        ! Setting PES after attempted hop.
        inext = istate

      ! Nonadiabatic tunneling case 3: E < Et ==> Rejected hop.
      ! Velocity is not rescaled as particle does not hop, so it
      ! stays at current state (same vp and same xp).
      else
        if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop = nfrust_hop + 1
          nFrust(istate,inext) = nFrust(istate,inext) + 1
        end if

        if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = 0
          nFrustR(istate,inext) = nFrustR(istate,inext) + 0
        end if

        if(ldtl)WRITE(nrite_hopp,'(" NT: Not enough kinetic energy &
                                    (E < Et), rejected")')
        if((ctime*dt/41340.d0)>avetim)then
          nJumpFail(istate,inext) = nJumpFail(istate,inext) + 1
        end if

        inext = istate
      end if
    end if

    contains

    ! Locate index of arr at min(arr(lb:ub)).
    integer function MinLoc(arr,lb,ub)
      implicit none
      real*8,  intent(in)  :: arr(:)
      real*8               :: tol
      integer              :: val
      integer              :: i, lb, ub
      tol = 1.0d5
      val = nint(arr(lb)*tol)/tol
      minLoc = lb
      do i=lb,ub
        if(arr(i) .le. val)then
          minLoc = i; val = nint(arr(i)*tol)/tol
        end if
      end do
    end function MinLoc
    ! Locate index of arr where arr(idx) == val.
    integer function LocIdx(arr,val, lb,ub)
      implicit none
      real*8,  intent(in)  :: arr(:), val
      integer              :: i, lb, ub
      locIdx = lb
      do i=lb+1,ub-1
        if((arr(i-1) .gt. val) .and. (val .gt. arr(i+1)))then
          locIdx = i; exit
        end if
        if((arr(i+1) .gt. val) .and. (val .gt. arr(i-1)))then
          locIdx = i; exit
        end if
      end do
    end function LocIdx
  end subroutine ZNCorrect
!**********************************************************************
end module zhunakamura_module

