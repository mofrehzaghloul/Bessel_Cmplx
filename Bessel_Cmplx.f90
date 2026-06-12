!*******************************************************************************
! MODULE: Bessel_Cmplx
!
! PURPOSE
! -------
! This module evaluates Bessel and modified Bessel functions for complex
! arguments z and real order nu:
!
!   I_nu(z) : modified Bessel function of the first kind
!   J_nu(z) : Bessel function of the first kind
!   K_nu(z) : modified Bessel function of the second kind
!   Y_nu(z) : Bessel function of the second kind
!
! The implementation supports positive and negative real orders, including
! integer and non-integer values, subject to the usual branch conventions for
! complex powers, logarithms, and square roots.
!
! MAIN USER-CALLABLE ROUTINES
! ---------------------------
! The following routines are declared PUBLIC and are intended to be called from
! other programs:
!
!   call i_nu_of_z(nu, z, cbi, ierr)   ! returns I_nu(z)
!   call j_nu_of_z(nu, z, cbj, ierr)   ! returns J_nu(z)
!   call k_nu_of_z(nu, z, cbk, ierr)   ! returns K_nu(z)
!   call y_nu_of_z(nu, z, cby, ierr)   ! returns Y_nu(z)
!
! All other routines are PRIVATE helper routines used internally by the module.
!
! PRECISION MODEL
! ---------------
! The working floating-point kind is controlled by the module set_rk:
!
!   rk : working kind used by all computations and routine interfaces
!   dp : double precision kind, imported for consistency
!   qp : quadruple precision kind, imported for consistency
!
! The module assumes that constants such as pi, eps, cone, czero, one, two,
! log_rmax, log_rmin, etc. are supplied by parameters.f90.
!
! NUMERICAL METHODS
! -----------------
! The module uses a region-based strategy. Depending on the magnitudes of z and
! nu, the dispatcher selects among:
!
!   1. Power series for small and intermediate |z|;
!   2. Large-|z| asymptotic expansions;
!   3. Large-order uniform asymptotic expansions;
!   4. Backward recurrence in the order for I_nu;
!   5. Forward recurrence in the order for K_nu;
!   6. Reflection, connection, and analytic-continuation formulas.
!
! BRANCH AND CONTINUATION CONVENTIONS
! -----------------------------------
! The intrinsic complex log and sqrt functions are used on their principal
! branches. Analytic continuation into the left half-plane is handled using
! explicit phase factors. On the negative real axis, where a bank must be
! selected, the implementation adopts the upper-bank convention when Im(z)=0.
!
! ERROR FLAGS
! -----------
!   ierr =  0 : successful computation
!   ierr =  1 : overflow or no reliable computation
!   ierr = -1 : underflow; result is treated as zero
!   ierr = 99 : input error, for example K_nu(0)
!   ierr = -2 : internal dispatcher condition; should not occur through normal
!               use of the public routines
!
! AUTHOR
! ------
! Mofreh R. Zaghloul
! Department of Physics, UAE University
! Steven Johnson
! Deptartment of Mathematics, MIT
! June 2026
!***********************************************************************
    Module Bessel_Cmplx
      Use set_rk, Only: rk, dp, qp
      Use, Intrinsic :: ieee_arithmetic
      Use, Intrinsic :: ieee_exceptions

      Implicit None

      Private
       ! These are the only routines intended to be called directly by users of the
       ! module. The remaining routines are internal implementation details.
      Public :: i_nu_of_z, j_nu_of_z, k_nu_of_z, y_nu_of_z
    
      Include 'parameters.f90'    
     
      Real (rk), Parameter :: series_border = 2.20_rk + rk_by_qp*1.0_rk !Used mainly by 
     !K_nu to decide when a small-z series is appropriate.
     
      Real (rk), Parameter :: z_inf_border = 625.0_rk !Used by K_nu to separate 
     !intermediate-z and large-z regions.
      Real (rk), Parameter :: inf_nu_const = 25.0_rk + rk_by_qp*65.0_rk !Used to detect 
     ! the large-order region, where uniform asymptotic expansions are preferred.
      Complex (rk), Parameter :: j1_pi = j1*pi ! Frequently used constant i*pi.
    
    Contains

!==============================================================================
!                         1. Main user-callable routines
!==============================================================================

!*******************************************************************************
! SUBROUTINE: i_nu_of_z
!
! PURPOSE
! -------
! Compute the modified Bessel function of the first kind, I_nu(z), for real
! order nu and complex argument z.
!
! CALLING SEQUENCE
! ----------------
!   call i_nu_of_z(nu, z_in, cbi, ierr)
!
! INPUTS
! ------
!   nu   : real(rk)
!          Order of the modified Bessel function. The order may be positive,
!          negative, integer, or non-integer.
!
!   z_in : complex(rk)
!          Complex argument.
!
! OUTPUTS
! -------
!   cbi  : complex(rk)
!          Computed value I_nu(z_in).
!
!   ierr : integer
!          Status flag.
!
! METHOD
! ------
! This is the top-level dispatcher for I_nu(z).
!
! Positive orders are passed directly to i_abs_nu_of_z. Integer negative orders
! are also evaluated through the positive-order relation I_{-n}(z)=I_n(z).
!
! For negative non-integer orders, the routine applies the appropriate analytic
! continuation and connection formula. Depending on the region, it selects:
!
!   1. the small-z power series;
!   2. the negative-order uniform asymptotic expansion;
!   3. the large-|z| asymptotic expansion;
!   4. backward recurrence combined with K_nu through the connection formula.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - i_abs_nu_of_z: handles nonnegative orders and integer negative orders.
!   - bessel_series_core: evaluates the small-z I_nu(z) series for negative
!     non-integer order after continuation if needed.
!   - i_neg_nu_unif_sum: large-order uniform asymptotic treatment for
!     negative order.
!   - i_nu_inf_z: large-z asymptotic evaluation.
!   - i_nu_bk_recurr: backward recurrence for positive |nu|.
!   - k_nu_intrmed_z: supplies K_|nu|(z) for the I_{-nu} connection formula.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: pi, eps, cone, j1, one, two, zero, c1,
!     abs_z_brdr1, and related logarithmic/threshold constants.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, cos, floor, max, real, sin.
!*******************************************************************************

      Subroutine i_nu_of_z(nu, z_in, cbi, ierr)

        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbi
        Integer, Intent (Out) :: ierr
        Real (rk) :: xx, yy, abs_z, fnu, y_sign, abs_nu, small_z_border, arg
        Complex (rk) :: z, phase, cbk
        Real (rk), Parameter :: two_by_pi = two/pi
        Integer :: n

        n = floor(nu)
        fnu = nu - n
        arg = fnu*pi
        If (abs(nu-n)<eps) Then
          Call i_abs_nu_of_z(abs(nu), z_in, cbi, ierr)
          Return
        End If

        If (nu>zero) Then
          Call i_abs_nu_of_z(nu, z_in, cbi, ierr)
          Return
        End If

        If (nu<zero) Then
          abs_nu = abs(nu)
          xx = real(z_in, kind=rk)
          yy = aimag(z_in)

          If (yy==zero) Then
            y_sign = one ! convention for real axis: approach from above
          Else
            y_sign = yy/abs(yy)
          End If

          phase = cone
          z = z_in
          If (xx<zero) Then
            z = -z_in
            phase = ((-1)**n)*(cos(arg)+y_sign*j1*sin(arg))
          End If

          abs_z = abs(z_in)

          small_z_border = 324.0_rk + (8.0_rk)*abs_nu

          If (abs_z*abs_z<=small_z_border) Then
!            Call i_nu_series(nu, z, cbi, ierr)
            Call bessel_series_core(nu, z, +1, cbi, ierr)

            cbi = phase*cbi

          Else If (nu<=-(c1+abs_z)) Then
            Call i_neg_nu_unif_sum(abs_nu, z, cbi, ierr)
            cbi = phase*cbi

          Else If (abs_z*abs_z>max(small_z_border,abs_z_brdr1**2) .And. (two*abs_z)>=nu*nu) Then
            Call i_nu_inf_z(nu, z, cbi, ierr)
            cbi = phase*cbi


          Else

            Call i_nu_bk_recurr(abs_nu, z, cbi, ierr)
            Call k_nu_intrmed_z(abs_nu, z, cbk, ierr)
            n = floor(abs_nu)
            fnu = abs_nu - n
            cbi = phase*(cbi+two_by_pi*((-one)**n)*sin(fnu*pi)*cbk)


          End If
        End If
        Return

      End Subroutine i_nu_of_z




!*******************************************************************************
! SUBROUTINE: j_nu_of_z
!
! PURPOSE
! -------
! Compute the Bessel function of the first kind, J_nu(z), for real order nu
! and complex argument z.
!
! CALLING SEQUENCE
! ----------------
!   call j_nu_of_z(nu, z_in, cbj, ierr)
!
! INPUTS
! ------
!   nu   : real(rk)
!          Order of the Bessel function.
!
!   z_in : complex(rk)
!          Complex argument.
!
! OUTPUTS
! -------
!   cbj  : complex(rk)
!          Computed value J_nu(z_in).
!
!   ierr : integer
!          Status flag. The value is propagated from the selected numerical
!          route, usually i_nu_of_z or the shared power-series kernel.
!
! METHOD
! ------
! For nonnegative orders and sufficiently small |z|, J_nu is evaluated directly
! by the shared small-z series kernel using the alternating sign appropriate
! for J_nu.
!
! Outside that region, the routine uses the analytic relations
!
!   Im(z) >= 0 : J_nu(z) = exp(+i*pi*nu/2) I_nu(-i z)
!   Im(z) <  0 : J_nu(z) = exp(-i*pi*nu/2) I_nu(+i z)
!
! To reduce loss of significance in the phase factor, nu is split into an
! integer part and a small fractional part.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - bessel_series_core: evaluates the small-z J_nu(z) series for
!     nonnegative orders.
!   - i_nu_of_z: evaluates I_nu(+-i z), which is then mapped to J_nu(z)
!     using phase factors.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: pi, half_pi, cone, czero, j1, one, zero.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, cmplx, conjg, cos, modulo, nint, real, sin.
!*******************************************************************************
Subroutine j_nu_of_z(nu, z_in, cbj, ierr)
  Real (rk), Intent (In) :: nu
  Complex (rk), Intent (In) :: z_in
  Complex (rk), Intent (Out) :: cbj
  Integer, Intent (Out) :: ierr

  Complex (rk) :: z, zi, cbi, csgn, base, frac_phase, phase
  Real (rk) :: fnu, ang_small, small_z_border, abs_nu, y_sign
  Real (rk) :: abs_z2, xx, yy, arg, signpn
  Integer :: n, m4, ipar
  Logical :: need_conj, use_small_series

  ierr = 0
  cbj = czero
  xx = real(z_in, kind=rk)
  yy = aimag(z_in)
  abs_z2 = xx*xx + yy*yy
  abs_nu = abs(nu)

  ! -- stable integer/fractional split: n = nearest integer, fnu in [-1/2,1/2]
  n = nint(nu)
  fnu = nu - n   ! now |fnu| <= 0.5 (up to the rounding ties of nint)

  ! Decide whether to use the small-z J-series:
  ! Only for nonnegative orders; negative orders route via I_v mapping.
  small_z_border = 324.0_rk + (8.0_rk)*abs_nu
  use_small_series = (nu >= 0.0_rk .and. abs_z2 <= small_z_border)

  ! ----- Map to I_v on the appropriate half-plane
  If (yy >= 0.0_rk) Then
    zi = -j1*z_in
    y_sign = +one   ! corresponds to exp(+i*pi*nu/2)
  Else
    zi = +j1*z_in
    y_sign = -one   ! corresponds to exp(-i*pi*nu/2)
  End If

  If (use_small_series) Then
    ! ----- Region 1: Small z Series (safe for nu >= 0)
    phase = cone
    z = z_in

    If (nu < 0.0_rk .And. xx < 0.0_rk) Then
      z = -z_in
      arg = fnu * pi
      ! use parity to generate (-1)^n robustly:
      ipar = modulo(n,2)
      If (ipar == 0) Then
        signpn = +1.0_rk
      Else
        signpn = -1.0_rk
      End If
      phase = cmplx(signpn* cos(arg), signpn* y_sign * sin(arg), kind=rk)
    End If

    Call bessel_series_core(nu, z, -1, cbj, ierr)
    cbj = phase * cbj
    Return
  Else
    ! ----- Mapping via I_nu with stable split
    m4 = modulo(n, 4)

    Select Case (m4)
    Case (0)
      base = cmplx(+one, zero, kind=rk)
    Case (1)
      base = cmplx(zero, +one, kind=rk)
    Case (2)
      base = cmplx(-one, zero, kind=rk)
    Case (3)
      base = cmplx(zero, -one, kind=rk)
    End Select
    If (y_sign < 0.0_rk) base = conjg(base)

    ang_small = y_sign * fnu * half_pi
    frac_phase = cmplx(cos(ang_small), sin(ang_small), kind=rk)
    csgn = base * frac_phase

    ! Evaluate I_v; call directly (use conjugation trick in k_nu if required)
    Call i_nu_of_z(nu, zi, cbi, ierr)
    If (ierr /= 0) Then
      cbj = cbi
      Return
    End If

    cbj = cbi * csgn
  End If

  Return
End Subroutine j_nu_of_z


!*******************************************************************************
! SUBROUTINE: k_nu_of_z
!
! PURPOSE
! -------
! Compute the modified Bessel function of the second kind, K_nu(z), for real
! order nu and complex argument z.
!
! CALLING SEQUENCE
! ----------------
!   call k_nu_of_z(anu, z_in, cbk, ierr)
!
! INPUTS
! ------
!   anu  : real(rk)
!          Requested order. The routine uses the symmetry K_{-nu}(z)=K_nu(z),
!          so the computation is performed with abs(anu).
!
!   z_in : complex(rk)
!          Complex argument. The value z_in = 0 is singular and returns
!          ierr = 99.
!
! OUTPUTS
! -------
!   cbk  : complex(rk)
!          Computed value K_nu(z_in).
!
!   ierr : integer
!          Status flag.
!
! METHOD
! ------
! The routine first normalizes the order using K_{-nu}=K_nu. For arguments in
! the left half-plane, it maps the computation to -z and then applies the
! analytic-continuation formula
!
!   K_nu(z) = exp(-i*s*pi*nu) K_nu(-z) - i*s*pi I_nu(-z),
!
! where s = sign(Im(z)). If Im(z)=0, the upper-bank convention s=+1 is used.
!
! Region selection is based on |z| and |nu|:
!
!   1. small |z|          -> k_nu_small_z
!   2. large |nu|         -> unif_sum_core, K branch
!   3. intermediate |z|   -> k_nu_intrmed_z
!   4. large |z|          -> k_nu_inf_z
!   5. hybrid case        -> asymptotic seed plus forward recurrence
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - k_nu_small_z: small-z evaluation of K_nu(z).
!   - unif_sum_core: large-order uniform asymptotic evaluation of K_nu(z)
!     using kind_flag = -1.
!   - k_nu_intrmed_z: intermediate-z evaluation of K_nu(z).
!   - k_nu_inf_z: large-z asymptotic evaluation of K_nu(z).
!   - i_abs_nu_of_z: evaluates I_nu(z) needed by the left-half-plane
!     analytic-continuation formula.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - Module parameters defined in this file: series_border, z_inf_border,
!     inf_nu_const, and j1_pi.
!   - From parameters.f90: pi, cone, czero, j1, one, two, zero, and
!     ieee-compatible limits.
!   - From intrinsic IEEE support: ieee_value is used to return NaN for the
!     invalid input z = 0.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, cos, floor, real, sign, sin, sqrt.
!*******************************************************************************

      Subroutine k_nu_of_z(anu, z_in, cbk, ierr)

        Real (rk), Intent (In) :: anu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbk
        Integer, Intent (Out) :: ierr

        Complex (rk) :: k_nu, k_nu_min_1, two_by_z, z, cbi, phase
        Real (rk) :: abs_z, fnu
        Real (rk) :: xx, yy, nu, nutmp, nutmp_1, fnu_pi, y_sign
        Integer :: i, n, ntmp, ntmp_1, ierr_i
        Logical :: is_left

        If (z_in==czero) Then
          ierr = 99
          cbk = cmplx(ieee_value(1.0_rk,ieee_quiet_nan), ieee_value(1.0_rk,ieee_quiet_nan), kind=rk)
          Return
        End If

        nu = abs(anu) !: k_v(z)=k_-v(z) (even in v)

        xx = real(z_in, kind=rk)
        yy = aimag(z_in)
        abs_z = sqrt(xx*xx+yy*yy)

        z = z_in
        phase = cone
        is_left = (xx<zero)
        y_sign = one
        If (is_left) Then
          z = -z_in
          n = floor(nu)
          fnu = nu - n
          fnu_pi = fnu*pi
          y_sign = sign(one, yy)
          phase = (((-one)**n)*(cos(fnu_pi)-y_sign*j1*sin(fnu_pi)))
        End If

!:-- Region 1: Power Series & Forward Recurrence
!:
        If (abs_z<=series_border) Then
          Call k_nu_small_z(nu, z, cbk, ierr)

! :--Region 2: Asymptotic expansion for nu>nu_th (threshold) & |z|> series_border
! :
        Else If (nu>=(inf_nu_const+abs_z)) Then
!          Call k_nu_unif_sum(nu, z, cbk, ierr)
          Call unif_sum_core(nu, z, -1, cbk, ierr)
! :--Region 3: Intermediate region (Miller's Method & Forward Recurrence)
! :
        Else If (abs_z<=z_inf_border) Then
          Call k_nu_intrmed_z(nu, z, cbk, ierr)

! :--Region 4 
! :   Aymp. Expansion for |z|-->Inf with forward recurrence for nu^2>2|z| and nu<nu_th
        Else If (two*abs_z>=nu*nu) Then
          Call k_nu_inf_z(nu, z, cbk, ierr)

        Else

          If (abs(xx)>0.2_rk*abs(yy)) Then
!            Call k_nu_unif_sum(nu, z, cbk, ierr)
            Call unif_sum_core(nu, z, -1, cbk, ierr)
          Else
            n = floor(nu)
            fnu = nu - n
            ntmp = floor(sqrt(two*abs_z)) - 2
            ntmp_1 = ntmp - 1
            nutmp = ntmp + fnu
            nutmp_1 = ntmp_1 + fnu

            Call k_nu_inf_z(nutmp, z, k_nu, ierr)
            Call k_nu_inf_z(nutmp_1, z, k_nu_min_1, ierr)

            two_by_z = two/z
            Do i = ntmp, n - 1
              cbk = two_by_z*(fnu+i)*k_nu + k_nu_min_1
              k_nu_min_1 = k_nu
              k_nu = cbk
            End Do
          End If
        End If

        If (is_left) Then ! left half of the domain
          Call i_abs_nu_of_z(nu, z, cbi, ierr_i)
          cbk = (phase*cbk-y_sign*j1_pi*cbi)
        End If
!        If (yy==0.0_rk .And. xx>0.0_rk) cbk = cmplx(real(cbk,kind=rk), 0.0_rk, kind=rk)
        Return
      End Subroutine k_nu_of_z

!*******************************************************************************
! SUBROUTINE: y_nu_of_z
!
! PURPOSE
! -------
! Compute the Bessel function of the second kind, Y_nu(z), for real order nu
! and complex argument z.
!
! CALLING SEQUENCE
! ----------------
!   call y_nu_of_z(nu, z_in, cby, ierr)
!
! INPUTS
! ------
!   nu   : real(rk)
!          Order of the Bessel function.
!
!   z_in : complex(rk)
!          Complex argument.
!
! OUTPUTS
! -------
!   cby  : complex(rk)
!          Computed value Y_nu(z_in).
!
!   ierr : integer
!          Status flag. A nonzero value is propagated from the underlying
!          J_nu or K_nu computation.
!
! METHOD
! ------
! For real nu, the function satisfies conjugation symmetry:
!
!   Y_nu(conjg(z)) = conjg(Y_nu(z)).
!
! Therefore, points in the lower half-plane are evaluated by conjugating the
! input into the upper half-plane, computing the function there, and then
! conjugating the result back.
!
! Away from integer order, the routine uses the reflection identity
!
!   Y_nu(z) = [cos(pi*nu) J_nu(z) - J_{-nu}(z)] / sin(pi*nu).
!
! Near integer order, this formula is numerically ill-conditioned because
! sin(pi*nu) is small. In that case the routine computes H_nu^(1)(z) using
! K_nu(-i z), and then recovers
!
!   Y_nu(z) = [H_nu^(1)(z) - J_nu(z)] / i.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - j_nu_of_z: evaluates J_nu(z) and J_{-nu}(z) for the reflection formula.
!   - k_nu_of_z: evaluates K_nu(-i z) for the near-integer-order fallback
!     through the Hankel-function relation.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: pi, half_pi, sqrt_eps, cone, czero, j1, one, two.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, cmplx, conjg, cos, modulo, nint, sin.
!*******************************************************************************

      Subroutine y_nu_of_z(nu, z_in, cby, ierr) 
        Implicit None
        Real    (rk), Intent (In)  :: nu
        Complex (rk), Intent (In)  :: z_in
        Complex (rk), Intent (Out) :: cby
        Integer        , Intent (Out) :: ierr

! locals
        Integer :: n, ipar, m4, ierrjp, ierrjm, ierrk0
        Real    (rk) :: f, s_pi, c_pi, tau, phi_half
        Real    (rk) :: cf, sf
        Complex (rk) :: jp, jm, w, k0, h1
        Complex (rk) :: base, rot_half, e_mhalf
        Complex (rk) :: z              ! <--- work variable for z
        Logical :: need_conj
        Logical :: conj_input          ! <--- did we conjugate the input?

        ierr = 0
        cby  = cmplx(0.0_rk, 0.0_rk, rk)

!------------------------ handle lower half-plane via conjugation -------
! For real nu, Y_nu(conjg(z)) = conjg( Y_nu(z) )
        conj_input = .false.
        z = z_in
        if (aimag(z_in) < 0.0_rk) then
          conj_input = .true.
          z = conjg(z_in)
        end if

!------------------------- stable sin(pi*nu), cos(pi*nu) ----------------
! Split nu = n + f with f in [-1/2, 1/2]
        n = nint(nu)
        f = nu - n
        if (f >  0.5_rk) then
          n = n + 1
          f = f - 1.0_rk
        else if (f < -0.5_rk) then
          n = n - 1
          f = f + 1.0_rk
        end if

! sin(pi*nu) = (-1)^n sin(pi*f), cos(pi*nu) = (-1)^n cos(pi*f)
        ipar = modulo(n,2)          ! 0 or 1 even for negative n
        if (ipar == 0) then
          c_pi =  cos(pi*f)
          s_pi =  sin(pi*f)
        else
          c_pi = -cos(pi*f)
          s_pi = -sin(pi*f)
        end if

!----------------------- primary reflection identity --------------------
        tau = sqrt_eps              ! e.g. sqrt(epsilon(1.0_rk))
        if (abs(s_pi) >= tau) then
          ierrjp = 0
          ierrjm = 0

          call j_nu_of_z(nu,  z, jp, ierrjp)
          if (ierrjp /= 0) then
            ierr = ierrjp
            cby  = jp
            goto 999
          end if

          call j_nu_of_z(-nu, z, jm, ierrjm)
          if (ierrjm /= 0) then
            ierr = ierrjm
            cby  = jm
            goto 999
          end if

          cby = (c_pi*jp - jm) / s_pi
          goto 999
        end if

!===================== fallback: K-mapping via H^{(1)} ==================
! w = -i z
        w = -j1*z

! base = exp(i*pi*n/2) from exact table
        m4 = modulo(n,4)
        select case (m4)
        case (0)
          base = cone
        case (1)
          base = j1
        case (2)
          base = -cone
        case (3)
          base = -j1
        end select

        phi_half = half_pi * f
        cf       = cos(phi_half)
        sf       = sin(phi_half)
        rot_half = cmplx(cf, sf, rk)          ! exp(+i*pi*f/2)

! exp(-i*pi*nu/2) = conjg( base * rot_half ) = conjg( exp(+i*pi*nu/2) )
        e_mhalf = conjg(base * rot_half)

! Evaluate K_nu(w) robustly, via conjugation if needed
        ierrk0    = 0
        need_conj = (aimag(w) < 0.0_rk)

        if (need_conj) then
          call k_nu_of_z(abs(nu), conjg(w), k0, ierrk0)
          if (ierrk0 == 0) k0 = conjg(k0)
        else
          call k_nu_of_z(abs(nu), w, k0, ierrk0)
        end if

        if (ierrk0 /= 0) then
          ierr = ierrk0
          cby  = k0
          goto 999
        end if

! H^(1)_nu(z) = (2/(pi*i)) * exp(-i*pi*nu/2) * K_nu(w),
! 2/(pi*i) = -2*i/pi
        h1 = -j1 * (two/pi) * e_mhalf * k0

! J_nu(z) for final Y
        ierrjp = 0
        call j_nu_of_z(nu, z, jp, ierrjp)
        if (ierrjp /= 0) then
          ierr = ierrjp
          cby  = jp
          goto 999
        end if

! Y_nu(z) = (H^(1)_nu(z) - J_nu(z))/i = -i*(H^(1)_nu(z) - J_nu(z))
        cby = -j1 * (h1 - jp)

999     continue
! If we flipped z into the upper half-plane, flip the result back
        if (conj_input) cby = conjg(cby)

        return
      End Subroutine y_nu_of_z
      
      
      
!==============================================================================
!                     2. Private routines for I_nu(z)
!==============================================================================
      
      
!*******************************************************************************
! SUBROUTINE: i_abs_nu_of_z
!
! PURPOSE
! -------
! Compute I_nu(z) for nonnegative order nu.
!
! INPUTS
! ------
!   nu   : real(rk), assumed to be nonnegative by the caller
!   z_in : complex(rk)
!
! OUTPUTS
! -------
!   cbi  : complex(rk), computed I_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! This private helper is the main positive-order dispatcher. It selects among:
!
!   1. power series for small/intermediate |z|;
!   2. large-|z| asymptotic expansion;
!   3. uniform asymptotic expansion for large order;
!   4. backward recurrence in the remaining intermediate region.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - bessel_series_core: small- and intermediate-z I_nu(z) power series.
!   - i_nu_inf_z: large-z asymptotic expansion.
!   - unif_sum_core: large-order uniform asymptotic expansion for I_nu(z)
!     using kind_flag = +1.
!   - i_nu_bk_recurr: backward recurrence in the order.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: eps, cone, czero, one, two, zero, c1,
!     abs_z_brdr1, log_rmax, log_rmin, and rk_by_qp-dependent thresholds.
!   - From intrinsic IEEE support: ieee_value is used when constructing
!     infinity-valued overflow results.
!
! Fortran intrinsic procedures used:
!   - aimag, cmplx, max, real, sqrt.
!*******************************************************************************    
      Subroutine i_abs_nu_of_z(nu, z_in, cbi, ierr)

        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbi
        Integer, Intent (Out) :: ierr

        Real (rk) :: xx, yy, abs_z, small_z_brdr, abs_z2


! :..
! :.. Check input -->
! :.. Only (+) or zero real orders and (-) integers orders
        ierr = 0
        If (nu<zero .And. (nu-floor(nu))/=zero) ierr = 99
        If (ierr/=0) Return

! :..
!        abs_z = abs(z_in)
!        abs_z2 = abs_z*abs_z
        xx = real(z_in, kind=rk)
        yy = aimag(z_in)
        abs_z2 = xx*xx + yy*yy
        abs_z = sqrt(abs_z2)


! :-
! :--Region 1: Power Series small to intermediate |z|
! :
        small_z_brdr = 324.0_rk + 8.0_rk*nu

        If (abs_z2<=small_z_brdr) Then
          Call bessel_series_core(nu, z_in, +1, cbi, ierr)

          If (ierr==1) Then
            cbi = cmplx(ieee_value(1.0_rk,ieee_positive_inf), ieee_value(1.0_rk,ieee_positive_inf), kind=rk)
          Else If (ierr==-1) Then
            cbi = czero
          End If

          Return

! :-
! :--Region 2 (|z|^2> max(ABS_Z_BRDR1^2, small_z_brdr) &  2|z|>=v^2)
! :   Aymptotic Expansion for |z|-->Inf

        Else If (abs_z2>max(abs_z_brdr1*abs_z_brdr1,small_z_brdr) .And. (two*abs_z)>=nu*nu) Then

          Call i_nu_inf_z((nu), z_in, cbi, ierr)
          If (ierr==1) Then
            cbi = cmplx(ieee_value(1.0_rk,ieee_positive_inf), ieee_value(1.0_rk,ieee_positive_inf), kind=rk)
          Else If (ierr==-1) Then
            cbi = czero
          End If

          Return

! :--Region 3: Uniform Asymptotic Expansion (nu->+oo) (v>=(c1+|z|))
        Else If (nu>=(c1+abs_z)) Then
!   Call i_nu_unif_sum(nu, z_in, cbi, ierr)
          Call unif_sum_core(nu, z_in, +1, cbi, ierr)
          If (ierr==1) Then
            cbi = cmplx(ieee_value(1.0_rk,ieee_positive_inf), ieee_value(1.0_rk,ieee_positive_inf), kind=rk)
          Else If (ierr==-1) Then
            cbi = czero
          End If

          Return

!!:
!:-- Region 3 Extension: Uniform Asymptotic Expansion
        Else If ((abs_z>1.8_rk*abs_z_brdr1 .And. abs(xx)>0.577350269189626_rk*abs(yy))) Then
!  Call i_nu_unif_sum(nu, z_in, cbi, ierr)
          Call unif_sum_core(nu, z_in, +1, cbi, ierr)
          If (ierr==1) Then
            cbi = cmplx(ieee_value(1.0_rk,ieee_positive_inf), ieee_value(1.0_rk,ieee_positive_inf), kind=rk)
          Else If (ierr==-1) Then
            cbi = czero
          End If
          Return
        Else


! :-
! :--Region 4:  Backward recurrence
          Call i_nu_bk_recurr(nu, z_in, cbi, ierr)
          If (ierr==1) Then
            cbi = cmplx(ieee_value(1.0_rk,ieee_positive_inf), ieee_value(1.0_rk,ieee_positive_inf), kind=rk)
          Else If (ierr==-1) Then
            cbi = czero
          End If
     
          Return

        End If
      End Subroutine i_abs_nu_of_z



!*******************************************************************************
! SUBROUTINE: i_nu_bk_recurr
!
! PURPOSE
! -------
! Compute I_nu(z) by stable backward recurrence in the order.
!
! INPUTS
! ------
!   nu   : real(rk), target nonnegative order
!   z_in : complex(rk), argument
!
! OUTPUTS
! -------
!   cbi  : complex(rk), computed I_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine chooses a starting order above the target order, computes two
! anchor values at adjacent higher orders, and then recurs downward using
!
!   I_{v-1}(z) = I_{v+1}(z) + (2v/z) I_v(z).
!
! The anchor values are computed by the power series for moderate |z| and by
! the uniform asymptotic expansion for larger |z|.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - bessel_series_core: provides high-order anchor values when |z| is not
!     too large.
!   - unif_sum_core: provides high-order anchor values when |z| is large.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: c1, eps, j1, one, two, zero, rk_by_qp, and the
!     complex constants used for zero and unit values.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, floor, mod, real.
!*******************************************************************************
      Subroutine i_nu_bk_recurr(nu, z_in, cbi, ierr)


        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbi
        Integer, Intent (Out) :: ierr
        Complex (rk) :: i_nu_sum, i_nu_sum_2, two_by_z
        Real (rk) :: abs_z, nu_tmp, f_nu, factor, nu_tmp_pls_1
        Integer :: n, m, nn
        Real (rk), Parameter :: abs_zs2 = (42.0_rk+rk_by_qp*65.0_rk)

        n = nu
        f_nu = nu - n
        abs_z = abs(z_in)
        ierr = 0
! Determine NN (border for series & backward recurrence)
! series_border =ABS_Z_BRDR1! four*sqrt(c1+one)

        If (abs_z<=abs_zs2) Then ! Anchors with the small/intermediate-Z series
          nn = floor(0.125_rk*abs_z*abs_z-42.5_rk) + 3
          nu_tmp = real(nn, kind=rk) + f_nu
          nu_tmp_pls_1 = nu_tmp + one

          Call bessel_series_core(nu_tmp, z_in, +1, i_nu_sum, ierr)
          If (ierr/=0) Return
          Call bessel_series_core(nu_tmp_pls_1, z_in, +1, i_nu_sum_2, ierr)
          If (ierr/=0) Return

        Else ! Anchors with the Uniform Asymptotic Expansion
          nn = floor(c1+abs_z+20.0_rk+rk_by_qp*20.0_rk)
          nu_tmp = real(nn, kind=rk) + f_nu
          nu_tmp_pls_1 = nu_tmp + one

          Call unif_sum_core(nu_tmp, z_in, +1, i_nu_sum, ierr)
          If (ierr/=0) Return

          Call unif_sum_core(nu_tmp_pls_1, z_in, +1, i_nu_sum_2, ierr)
          If (ierr/=0) Return
        End If

! Backward Recurrence: I_{v-1} = I_{v+1} + (2v/z) I_v
        two_by_z = two/z_in
        factor = nu_tmp_pls_1 ! ONE + F_NU + NN
        Do m = 1, nn - n
          factor = factor - one
          cbi = i_nu_sum_2 + (factor*two_by_z)*i_nu_sum
          i_nu_sum_2 = i_nu_sum
          i_nu_sum = cbi
        End Do

! Handle special cases for real/imaginary outputs
        If (f_nu==zero) Then
          If (aimag(z_in)==zero) Then
            cbi = real(cbi, kind=rk)
          Else If (real(z_in,kind=rk)==zero) Then
            If (mod(nu,two)==zero) Then
              cbi = real(cbi, kind=rk)
            Else If (mod(nu,two)/=zero) Then
              cbi = j1*aimag(cbi)
            End If
          End If
        End If

        Return
      End Subroutine i_nu_bk_recurr


!*******************************************************************************
! SUBROUTINE: i_nu_inf_z
!
! PURPOSE
! -------
! Evaluate I_nu(z) using a large-|z| asymptotic expansion.
!
! INPUTS
! ------
!   nu   : real(rk), order
!   z_in : complex(rk), argument
!
! OUTPUTS
! -------
!   cbi  : complex(rk), computed asymptotic approximation to I_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine evaluates the dominant exp(+z) contribution and, where needed,
! the recessive Stokes-weighted exp(-z) contribution. Coefficients depend on
! nu^2 and are therefore valid for either sign of nu, with an additional
! connection contribution included for negative orders.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: lg_one_sqrt_2_pi, rtpi, pi, two_pi, eps,
!     sqrt_rmin, log_rmax, log_rmin, cone, czero, j1, one, two, half,
!     zero, and related constants.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, cmplx, cos, exp, floor, log, min, mod, real, sin.
!*******************************************************************************

      Subroutine i_nu_inf_z(nu, z_in, cbi, ierr)

        Complex (rk), Intent (In) :: z_in
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (Out) :: cbi
        Integer, Intent (Out) :: ierr

        Complex (rk) :: pre_sum1_fact, t_k, i_nu_sum1, i_nu_sum2, den, eight_z
        Complex (rk) :: exp_nu_hf_ipi, lg_inu, rtpi_z_in
        Real (rk) :: abs_t_k, ak, abs_eight_z, re_exp, im_exp, arg, atol, abs_z
        Real (rk) :: abs_den, two_nu, four_nu_sq, asign, numr, yy, fnu, two_fnu, abs_nu, y_sign
        Integer :: n, j, int_two_fnu

        ierr = 0
        cbi = czero
        lg_inu = (z_in+lg_one_sqrt_2_pi) - half*log(z_in)
        If (real(lg_inu,kind=rk)>log_rmax) Then
          ierr = 1 ! overflow occurs
          Return

        Else If (real(lg_inu,kind=rk)<log_rmin) Then
          cbi = czero
          ierr = -1
          Return
        End If

        abs_nu = abs(nu)

        abs_z = abs(z_in)
        two_nu = two*abs_nu

        four_nu_sq = zero
        If ((two_nu)>sqrt_rmin) four_nu_sq = two_nu*two_nu
        eight_z = z_in*8.0_rk

        yy = aimag(z_in)
        y_sign = one
        exp_nu_hf_ipi = czero
        If (yy/=zero) Then
          y_sign = yy/abs(yy)
! -----------------------------------------------------------------------
! CALCULATE EXP(PI*(NU+1/2)*I) TO MINIMIZE LOSSES OF
! SIGNIFICANCE WHEN FNU OR N IS LARGE
! -----------------------------------------------------------------------
          n = abs_nu ! integer part of nu
          fnu = abs_nu - n
          arg = fnu*pi ! the fractional part of nu  scaled by pi

! *********Handling Special Cases: Integer and Half-Integer Values
! EXP(PI*(NU+1/2)*I)= EXP((PI/2)*(2*NU+1)*I)
! = exp(N*PI) * exp((fnu+1/2)*pi)
! = exp(N*PI) * exp(2*fnu+1)*(pi/2)
! = exp(N*PI) * exp(i*pi/2)*exp(2*fnu*pi/2)
! = [-1 if N odd]*[(0+i)* cos(fnu*pi)+i*sin(fnu*pi)]

          two_fnu = two*fnu
          int_two_fnu = floor(two_fnu)
          If ((two_fnu-int_two_fnu)/=zero) Then
            re_exp = -sin(arg)
            im_exp = cos(arg)
          Else
            Select Case (mod(int_two_fnu,4))
            Case (0)
              re_exp = zero
              im_exp = one
            Case (1)
              re_exp = -one
              im_exp = zero
            Case (2)
              re_exp = zero
              im_exp = -one
            Case (3)
              re_exp = one
              im_exp = zero
            End Select
          End If

          exp_nu_hf_ipi = cmplx(re_exp, y_sign*im_exp, kind=rk)
          exp_nu_hf_ipi = ((-one)**n)*exp_nu_hf_ipi
        End If

        numr = four_nu_sq - one ! 4*nu^2-1  numerator of the 1st term in the series
! (4k^2-(2k+1)^2), k=0
        den = eight_z ! 8*z

! series|
        asign = one
        i_nu_sum1 = cone ! Sum((-1)^k a_k /z^k)
        i_nu_sum2 = cone ! Sum(a_k /z^k)
        t_k = cone
        ak = zero


        Do j = 1, 30 !
          asign = -asign
          t_k = t_k*numr/den
          i_nu_sum2 = i_nu_sum2 + t_k
          i_nu_sum1 = i_nu_sum1 + t_k*asign
          den = den + eight_z ! updating the denominator of each term

          ak = ak + 8.0_rk ! term used to update the numerator
          numr = numr - ak ! updating the numerator of each term

          If (abs(real(t_k,rk))<=eps*min(abs(real(i_nu_sum1,rk)),abs(real(i_nu_sum2,rk))) .And. &
          abs(aimag(t_k))<=eps*min(abs(aimag(i_nu_sum1)),abs(aimag(i_nu_sum2)))) Exit



        End Do

        rtpi_z_in = rtpi/z_in
        pre_sum1_fact = exp(lg_inu)

        If (real(z_in,kind=rk)>-half*log_rmin) Then
          cbi = pre_sum1_fact*i_nu_sum1
        Else

          cbi = pre_sum1_fact*i_nu_sum1 + exp_nu_hf_ipi*i_nu_sum2*(rtpi_z_in/pre_sum1_fact)

!!=== extra term for v<0 (connection component)
          If (nu<zero) Then
            cbi = cbi + two*((-one)**n)*sin(arg)*rtpi_z_in*i_nu_sum2/pre_sum1_fact
          End If
!!============================================
          Return

        End If
        Return
      End Subroutine i_nu_inf_z
      
      

!*******************************************************************************
! SUBROUTINE: i_neg_nu_unif_sum
!
! PURPOSE
! -------
! Evaluate I_nu(z) for negative order using uniform asymptotics and the
! negative-order connection formula.
!
! INPUTS
! ------
!   nu   : real(rk), negative order or magnitude depending on caller context
!   z_in : complex(rk), argument
!
! OUTPUTS
! -------
!   cbi  : complex(rk), computed I_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine forms positive-order I and K asymptotic components using abs(nu)
! and combines them according to the connection formula for negative order.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: c_uni_asy, lg_one_sqrt_2_pi, pi, eps, cone,
!     czero, one, two, half, zero, log_rmax, and log_rmin.
!
! Fortran intrinsic procedures used:
!   - abs, exp, log, min, mod, real, sin, sqrt.
!*******************************************************************************
      Subroutine i_neg_nu_unif_sum(nu, z_in, cbi, ierr)

! Input and output variables
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbi
        Integer, Intent (Out) :: ierr
!
! Local variables
        Complex (rk) :: i_nu_sum, k_nu_sum, i_nu_sum_2, increment
        Complex (rk) :: taw_old, p_sqr, log_in, z_sqr, p, taw_k, eta, increment2
        Complex (rk) :: log_in_k, z
        Integer :: kk, j, m, ksign
        Real (rk) :: abs_nu, xx, arg

! Initialize variables
        ierr = 0
        cbi = czero

        abs_nu = abs(nu)
        xx = real(z_in, kind=rk)

        z = z_in
        If (xx<zero) Then
          z = -z_in
        End If

! Compute core values
        z_sqr = (z/abs_nu)**2
        p = cone/sqrt(cone+z_sqr)
        p_sqr = p*p
        taw_k = p/abs_nu
        eta = cone/p + log(taw_k*z/(cone+p))

! Compute logarithm of i prefactor
        log_in = lg_one_sqrt_2_pi + abs_nu*eta + half*log(taw_k)


! Compute logarithm of k prefactor
!        log_in_k = lg_one_sqrt_2_pi - abs_nu*eta + half*log(taw_k)
        log_in_k = log_in - two*abs_nu*eta
! Overflow check
        If (real(log_in,kind=rk)>log_rmax) Then
          ierr = 1 ! Overflow detected
          Return
        Else If (real(log_in,kind=rk)<log_rmin) Then
          ierr = -1 !underflow detected
          Return
        End If

! Initialize summation
        i_nu_sum = cone ! Start with 1 in complex form
        k_nu_sum = cone
        taw_old = cone
        m = 1
        ksign = 1

! Summation loop
        Do kk = 1, 28
          taw_old = taw_old*taw_k
          i_nu_sum_2 = czero
          ksign = -ksign
! Inner sum
          Do j = 0, kk
            m = m + 1
            i_nu_sum_2 = i_nu_sum_2 + c_uni_asy(m)*(p_sqr**(kk-j))
          End Do
          increment = taw_old*i_nu_sum_2
          increment2 = ksign*taw_old*i_nu_sum_2

! Update main sum
          i_nu_sum = i_nu_sum + increment
          k_nu_sum = k_nu_sum + increment2
! Convergence check
          If (min(abs(increment),abs(increment2))<eps) Exit
        End Do


! Final result
        arg = mod(abs_nu*pi, two*pi)
        cbi = exp(log_in)*i_nu_sum + exp(log_in_k)*(two*sin(arg)*k_nu_sum)


        Return

      End Subroutine i_neg_nu_unif_sum
      
      
!===============================================================================
!                       3. Private routines for K_nu(z)
!===============================================================================


!*******************************************************************************
! SUBROUTINE: k_nu_small_z
!
! PURPOSE
! -------
! Compute K_nu(z) in the small-|z| region.
!
! INPUTS
! ------
!   nu   : real(rk), nonnegative order
!   z_in : complex(rk), argument
!
! OUTPUTS
! -------
!   cbk  : complex(rk), computed K_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine handles half-integer special cases explicitly. Otherwise it
! evaluates a small-z expansion and, when necessary, applies recurrence to
! reach the requested order. Range checks protect against overflow in the
! leading behavior near z=0.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: cc, pi, eps, eps2, lg_sqrt_half_pi,
!     sqrt_half_pi, abs_log_rmin, cone, czero, one, two, half, quarter,
!     three_halfs, zero, and related logarithmic constants.
!
! Fortran intrinsic procedures used:
!   - aimag, cos, exp, gamma, log, real, sin, sqrt.
!*******************************************************************************

      Subroutine k_nu_small_z(nu, z_in, cbk, ierr)
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbk
        Integer, Intent (Out) :: ierr
        Complex (rk) :: c_k, p_k, meo, cosh_meo, sinh_meo, f_k, q_k, lg_2_by_z, quarter_z_sqr, exp_zin, inv_z
        Complex (rk) :: k_nu, k_nu_pls_1, half_z_to_nu, half_z, two_by_z, incr, lg_half_z, inv_half_z_to_nu

        Real (rk) :: g1, g2, gam_1_min_nu_inv, gam_1_pls_nu_inv, lg_k_nu
        Real (rk) :: pi_nu, anu, sin_pi_nu, pi_nu_by_sin_pi_nu, kk_pls_nu_inv, kk_min_nu_inv
        Real (rk) :: re_incr, sum, nu_sqr, a_i, half_by_nu, abs_incr, quarter_abs_z_sqr, abs_z
        Real (rk) :: kk_inv, xx, yy, delta, abs_z2
        Integer :: kk, n, i
!      
        ierr = 0
        xx = real(z_in, kind=rk)
        yy = aimag(z_in)
        abs_z2 = (xx*xx+yy*yy)
        abs_z = sqrt(abs_z2)

        If (nu>eps) Then
          lg_k_nu = lg_sqrt_half_pi + (nu-half)*log(nu) - nu*(0.306852819440054690582767878541823432_rk+log(abs_z))

          If (lg_k_nu>log_rmax) Then
            ierr = 1 !overflow
            Return
          End If
        End If

        inv_z = one/z_in
        If (nu<=half) Then
          If (abs(nu-half)<eps) Then
            If (xx<abs_log_rmin) Then
              exp_zin = exp(-xx)*cmplx(cos(yy), -sin(yy), rk)
              cbk = exp_zin*sqrt(half_pi*inv_z)
            Else
              cbk = czero
            End If
            Return
          End If

          anu = nu

        Else
          n = nu
          delta = nu - n

          If (delta==zero) Then
            anu = zero
            n = n - 1
          Else If (delta>=half) Then
            anu = delta - one
          Else
            anu = delta
            n = n - 1
          End If

        End If

        nu_sqr = anu*anu
        c_k = cone
        half_z = half*z_in
        lg_half_z = log(half_z)
        lg_2_by_z = -lg_half_z

        quarter_z_sqr = half_z*half_z !quarter*z_in*z_in
        quarter_abs_z_sqr = quarter*abs_z2

        meo = czero
        half_z_to_nu = cone

        pi_nu_by_sin_pi_nu = one
        gam_1_pls_nu_inv = one
        gam_1_min_nu_inv = one

        If (anu/=zero) Then
          meo = anu*lg_2_by_z
          half_z_to_nu = exp(-meo) !half_z**anu
          pi_nu = pi*anu
          sin_pi_nu = sin(pi_nu)
          pi_nu_by_sin_pi_nu = pi_nu/sin_pi_nu
!-------------
! Use     GAMma(1-nu)*GAMma(1+nu)=PI*nu/SIN(PI*nu)
!------------- 
          gam_1_pls_nu_inv = one/gamma(one+anu)
          gam_1_min_nu_inv = one/(pi_nu_by_sin_pi_nu*gam_1_pls_nu_inv)

        End If
        inv_half_z_to_nu = one/half_z_to_nu

        If (nu_sqr<=eps) Then
! Limit as nu -> 0
!          g1 = -euler_cnst - euler_zeta_comb*nu_sqr 
          g1 = -0.5772156649015328606065120900824024_rk - 1.814958152161772578033120334712550_rk*nu_sqr !

          q_k = half*half_z_to_nu !gam_1_min_nu*half_z_to_nu
          p_k = half*inv_half_z_to_nu ! gam_1_pls_nu/half_z_to_nu
          f_k = (g1*(p_k+q_k)+lg_2_by_z) !  f_k = (g1*half*(inv_half_z_to_nu+half_z_to_nu)+lg_2_by_z) 

        Else If (nu_sqr<=1.0E-6_rk) Then
          a_i = one
          sum = cc(1)
          Do i = 2, 29
            a_i = a_i*nu_sqr
            re_incr = cc(i)*a_i
            sum = sum + re_incr
            If ((re_incr*re_incr)<eps2) Exit
          End Do
          g1 = -sum
          cosh_meo = half*(half_z_to_nu+inv_half_z_to_nu)
          sinh_meo = cosh_meo - half_z_to_nu

          g2 = (half/anu)*(gam_1_min_nu_inv+gam_1_pls_nu_inv)
          q_k = half*half_z_to_nu/gam_1_min_nu_inv
          f_k = pi_nu_by_sin_pi_nu*(g1*cosh_meo+g2*sinh_meo)
          p_k = half/(gam_1_pls_nu_inv*half_z_to_nu) ! gam_1_pls_nu/half_z_to_nu
        Else
          half_by_nu = half/anu
          g1 = half_by_nu*(gam_1_min_nu_inv-gam_1_pls_nu_inv)
          cosh_meo = half*(half_z_to_nu+inv_half_z_to_nu)
          sinh_meo = cosh_meo - half_z_to_nu !  

          g2 = half_by_nu*(gam_1_min_nu_inv+gam_1_pls_nu_inv)
          q_k = half*half_z_to_nu/gam_1_min_nu_inv
          f_k = pi_nu_by_sin_pi_nu*(g1*cosh_meo+g2*sinh_meo)
          p_k = half/(gam_1_pls_nu_inv*half_z_to_nu) ! gam_1_pls_nu/half_z_to_nu

        End If


        If (nu<half) Then

          cbk = c_k*f_k
          abs_incr = one

          Do kk = 1, 100
            kk_pls_nu_inv = one/(kk+anu)
            kk_min_nu_inv = one/(kk-anu)
            f_k = (kk*f_k+p_k+q_k)*kk_pls_nu_inv*kk_min_nu_inv !/((kk-anu)*(kk+anu))
            kk_inv = one/kk
            c_k = c_k*quarter_z_sqr*kk_inv
            incr = c_k*f_k
            abs_incr = abs_incr*quarter_abs_z_sqr*kk_inv
            cbk = cbk + incr
            If (abs_incr<=eps) Exit
            q_k = q_k*kk_pls_nu_inv !/(kk+anu)
            p_k = p_k*kk_min_nu_inv ! /(kk-anu)
          End Do

        Else If (nu>half .And. nu<=three_halfs) Then
          If (abs(nu-three_halfs)<eps) Then
            If (xx<abs_log_rmin) Then
              exp_zin = exp(-xx)*cmplx(cos(yy), -sin(yy), rk)
              cbk = exp_zin*sqrt(half_pi*inv_z)*(one+inv_z)
            Else
              cbk = czero
            End If
            Return
          End If


          k_nu = c_k*p_k
          two_by_z = two*inv_z
          abs_incr = one

          Do kk = 1, 100
            kk_pls_nu_inv = one/(kk+anu)
            kk_min_nu_inv = one/(kk-anu)
            f_k = (kk*f_k+p_k+q_k)*kk_pls_nu_inv*kk_min_nu_inv !/((kk-anu)*(kk+anu))
            kk_inv = one/kk
            q_k = q_k*kk_pls_nu_inv !/(kk+anu)
            p_k = p_k*kk_min_nu_inv ! /(kk-anu)
            c_k = c_k*quarter_z_sqr*kk_inv
            incr = c_k*(p_k-kk*f_k)
            abs_incr = abs_incr*quarter_abs_z_sqr*kk_inv
            k_nu = k_nu + incr
            If (abs_incr<=eps) Exit

          End Do
          cbk = two_by_z*k_nu

        Else
          two_by_z = two*inv_z
          k_nu = c_k*f_k
          k_nu_pls_1 = c_k*p_k
          abs_incr = one
          Do kk = 1, 100
            kk_pls_nu_inv = one/(kk+anu)
            kk_min_nu_inv = one/(kk-anu)
            f_k = (kk*f_k+p_k+q_k)*kk_pls_nu_inv*kk_min_nu_inv !/((kk-anu)*(kk+anu))
            kk_inv = one/kk
            c_k = c_k*quarter_z_sqr*kk_inv
            p_k = p_k*kk_min_nu_inv ! /(kk-anu)
            incr = c_k*(p_k-kk*f_k)
            k_nu = k_nu + c_k*f_k
            k_nu_pls_1 = k_nu_pls_1 + incr
            abs_incr = abs_incr*quarter_abs_z_sqr*kk_inv
            If (abs_incr<=eps) Exit
            q_k = q_k*kk_pls_nu_inv !/(kk+anu)
          End Do


          k_nu = k_nu
          k_nu_pls_1 = two_by_z*k_nu_pls_1
          Do i = 1, n
            cbk = two_by_z*(anu+i)*k_nu_pls_1 + k_nu
            k_nu = k_nu_pls_1
            k_nu_pls_1 = cbk

          End Do
        End If
        Return
      End Subroutine k_nu_small_z



!*******************************************************************************
! SUBROUTINE: k_nu_intrmed_z
!
! PURPOSE
! -------
! Compute K_nu(z) in the intermediate-|z| region.
!
! INPUTS
! ------
!   nu   : real(rk), nonnegative order
!   z_in : complex(rk), argument
!
! OUTPUTS
! -------
!   cbk  : complex(rk), computed K_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine uses a Miller-type process to obtain a stable base value near a
! fractional order and then applies forward recurrence to reach the requested
! order.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: lg_sqrt_half_pi, half_pi, eps, sqrt_half_pi,
!     abs_log_rmin, one, two, three_halfs, half, quarter, zero, and
!     rk_by_qp-dependent thresholds.
!
! Fortran intrinsic procedures used:
!   - abs, aimag, atan2, cos, exp, log, real, sqrt.
!*******************************************************************************
      Subroutine k_nu_intrmed_z(nu, z_in, cbk, ierr)
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbk
        Integer, Intent (Out) :: ierr

        Complex (rk) :: c_sum, two_by_z, k_n, k_n_min_1, k_n_pls_1
        Complex (rk) :: k_nu, k_nu_pls_1, b_n, z_inv, two_z

        Integer :: n, i, m, m_sqr, m_max
        Real (rk) :: qrt_nu_sqr, ln_k_nu, anu, a_n, abs_z, theta, a, b, bracket, xx, delta
        Real (rk), Parameter :: factor1 = one + rk_by_qp*0.4_rk
        Real (rk), Parameter :: border1 = (40.0_rk+rk_by_qp*120.0_rk)
        Real (rk), Parameter :: lg_c = 32.0_rk + rk_by_qp*24.0_rk
        Integer, Parameter :: intg_2 = 7 + rk_by_qp*9

        ierr = 0
        abs_z = abs(z_in)

        If (nu>half) Then
          ln_k_nu = lg_sqrt_half_pi + (nu-half)*log(nu) - nu*(0.306852819440054690582767878541823432_rk+log(abs_z))
          If (ln_k_nu>log_rmax) Then
            ierr = 1 !overflow
            Return
          End If
        End If

        z_inv = one/z_in
        xx = real(z_in, kind=rk)
        theta = half_pi
        If (xx/=zero) theta = abs(atan2(aimag(z_in),xx))
        m_max = intg_2
        If (abs_z<border1) Then
          a = three*theta/(one+abs_z)
          b = (14.70_rk*theta)/(28.0_rk+abs_z)
          bracket = (lg_c+abs_z*cos(a)/(one+0.008_rk*abs_z))/(cos(b))
          m_max = ((0.12125_rk/abs_z)*bracket*bracket+three_halfs)*factor1
        End If


        two_z = z_in + z_in
        If (nu<=half) Then
          If (nu==half) Then
            If (xx<abs_log_rmin) Then
              cbk = exp(-z_in)*sqrt(half_pi*z_inv)
            Else
              cbk = czero
            End If
            Return
          Else

            anu = nu
            k_n_pls_1 = czero
            k_n = cmplx(eps, zero, kind=rk)
            c_sum = k_n_pls_1 + k_n
            qrt_nu_sqr = quarter - anu*anu

            Do m = m_max, 1, -1
              m_sqr = m*m
              a_n = (m_sqr-m+qrt_nu_sqr)
              b_n = (m+m+two_z)/(m+one)
              k_n_min_1 = (-k_n_pls_1+b_n*k_n)*(m_sqr+m)/a_n
              c_sum = c_sum + k_n_min_1
              k_n_pls_1 = k_n
              k_n = k_n_min_1
            End Do

            If (xx<abs_log_rmin) Then
              cbk = sqrt(half_pi*z_inv)*exp(-z_in)*k_n_min_1/c_sum
            Else
              cbk = czero
            End If
          End If

        Else
          n = nu
          delta = nu - n

          If (delta==zero) Then
            anu = zero
            n = n - 1
          Else If (delta>=half) Then
            anu = delta - one
          Else
            anu = delta
            n = n - 1
          End If

          qrt_nu_sqr = quarter - anu*anu
          k_n_pls_1 = czero
          k_n = cmplx(eps, zero, kind=rk)
          c_sum = k_n_pls_1 + k_n

          Do m = m_max, 1, -1
            m_sqr = m*m
            a_n = (m_sqr-m+qrt_nu_sqr)
            b_n = (m+m+two_z)/(m+one)
            k_n_min_1 = (-k_n_pls_1+b_n*k_n)*(m_sqr+m)/a_n
            c_sum = c_sum + k_n_min_1
            k_n_pls_1 = k_n
            k_n = k_n_min_1
          End Do

          If (xx<abs_log_rmin) Then
            k_nu = sqrt(half_pi*z_inv)*exp(-z_in)*k_n_min_1/c_sum
          Else
            k_nu = czero
          End If
          k_nu_pls_1 = k_nu*(z_in+anu+half-k_n_pls_1/k_n)*z_inv

          two_by_z = two*z_inv
          cbk = k_nu_pls_1
          Do i = 1, n
            cbk = two_by_z*(anu+i)*k_nu_pls_1 + k_nu
            k_nu = k_nu_pls_1
            k_nu_pls_1 = cbk
          End Do

        End If

      End Subroutine k_nu_intrmed_z


      
!*******************************************************************************
! SUBROUTINE: k_nu_inf_z
!
! PURPOSE
! -------
! Evaluate K_nu(z) using a large-|z| asymptotic expansion.
!
! INPUTS
! ------
!   nu   : real(rk), order
!   z_in : complex(rk), argument
!
! OUTPUTS
! -------
!   cbk  : complex(rk), computed asymptotic approximation to K_nu(z_in)
!   ierr : integer, status flag
!
! METHOD
! ------
! The expansion has the form
!
!   K_nu(z) ~ sqrt(pi/(2z)) exp(-z) * sum_k a_k/z^k.
!
! A real-axis fast path is used for positive real z. The general complex path
! computes the phase explicitly and uses the principal square root.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: lg_sqrt_half_pi, sqrt_half_pi, eps, sqrt_rmin,
!     log_rmin, cone, czero, one, half, eight, zero, and rk_by_qp.
!
! Fortran intrinsic procedures used:
!   - aimag, cmplx, cos, exp, log, merge, real, sin, sqrt.
!*******************************************************************************
      Subroutine k_nu_inf_z(nu, z_in, cbk, ierr)
        Complex (rk), Intent (In) :: z_in
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (Out) :: cbk
        Integer, Intent (Out) :: ierr

! Locals
        Complex (rk) :: t_k, k_nu_sum2, den, eight_z, exp_z
        Complex (rk) :: s, pref ! for prefactor via sqrt(z)
        Real (rk) :: re_lg_k_nu
        Real (rk) :: abs_t_k, ak, abs_eight_z, abs_z
        Real (rk) :: abs_den, inv_abs_den, two_nu, four_nu_sq, numr
        Real (rk) :: xx, yy, atol
        Complex (rk) :: inv_den
        Integer, Parameter :: jmax = 50 + 50*rk_by_qp
        Integer :: j

        ierr = 0
        xx = real(z_in, kind=rk)
        yy = aimag(z_in)
        abs_z = sqrt(xx*xx+yy*yy)

        re_lg_k_nu = -xx + lg_sqrt_half_pi - half*log(abs_z)
        If (re_lg_k_nu<log_rmin) Then
          ierr = -1
          cbk = czero
          Return
        End If

        two_nu = two*nu
        four_nu_sq = merge(two_nu*two_nu, zero, two_nu>sqrt_rmin)

        eight_z = eight*z_in
        abs_eight_z = eight*abs_z
        numr = four_nu_sq - one ! 4*nu^2 - 1
        den = eight_z
        atol = eps*abs(four_nu_sq-one)/abs_eight_z
        atol = atol*atol
!  ! === Real-axis fast track 
        If (yy==0.0_rk .And. xx>0.0_rk) Then
          pref = sqrt_half_pi/sqrt(xx)
          k_nu_sum2 = cone
          t_k = cone
          ak = zero
          abs_t_k = one
          Do j = 1, jmax
            t_k = t_k*numr/den
            k_nu_sum2 = k_nu_sum2 + t_k
            den = den + eight_z
            ak = ak + eight
            numr = numr - ak
            If (real(t_k,rk)**2+aimag(t_k)**2<=atol) Exit
          End Do

          cbk = pref*exp(-xx)*k_nu_sum2
          Return
        End If

! === General complex case
        exp_z = exp(-xx)*cmplx(cos(yy), -sin(yy), rk)

! Prefactor via one complex sqrt 
        s = sqrt(z_in)
        pref = sqrt_half_pi/s

        k_nu_sum2 = cone
        t_k = cone
        ak = zero
        abs_t_k = one

        Do j = 1, jmax
          t_k = t_k*numr/den
          k_nu_sum2 = k_nu_sum2 + t_k
          den = den + eight_z
          ak = ak + eight
          numr = numr - ak
          If (real(t_k,rk)*real(t_k,rk)+aimag(t_k)*aimag(t_k)<=atol) Exit
        End Do

        cbk = pref*exp_z*k_nu_sum2

      End Subroutine k_nu_inf_z
      

         
!==============================================================================
!                         4. Shared numerical kernels
!==============================================================================


!*******************************************************************************
! SUBROUTINE: bessel_series_core
!
! PURPOSE
! -------
! Shared power-series kernel for I_nu(z) and J_nu(z).
!
! INPUTS
! ------
!   nu        : real(rk), order
!   z_in      : complex(rk), argument
!   kind_flag : integer
!               +1 selects the I_nu series (non-alternating series)
!               -1 selects the J_nu series  (alternating via (-1)^k)
!
! OUTPUTS
! -------
!   cval : complex(rk), computed value
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine evaluates the prefactor
!
!   (z/2)^nu / Gamma(nu+1)
!
! using logarithms on the principal branch, then sums the standard coefficient
! recurrence. The only difference between the I and J series is the sign of the
! recurrence factor (z/2)^2:
!
!   I_nu :  +(z/2)^2
!   J_nu :  -(z/2)^2
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk. This routine also contains a local Use set_rk
!     statement for explicit kind visibility.
!   - From parameters.f90: lg_half, pi, two_pi, eps2, cone, czero, j1,
!     one, quarter, zero, and log_rmin.
!
! Fortran intrinsic procedures used:
!   - aimag, cmplx, cos, exp, log, log_gamma, merge, mod, real, sign, sin.
!*******************************************************************************
      Subroutine bessel_series_core(nu, z_in, kind_flag, cval, ierr)
        Use set_rk
        Implicit None
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Integer, Intent (In) :: kind_flag ! +1 for I, -1 for J
        Complex (rk), Intent (Out) :: cval
        Integer, Intent (Out) :: ierr

        Complex (rk) :: t_k, halfz_sq, halfz_sq_eff, log_in
        Real (rk) :: nu_pls_one, re_lg_in, im_lg_in, angle_check
        Real (rk) :: halfz_sq2, t_k_abs, exp_tmp, den_rcp
        Real (rk) :: atol, eps_den_rcp
        Integer :: kk

        ierr = 0
        cval = czero

        nu_pls_one = nu + one

! ----- log-prefactor: (z/2)^nu / Gamma(nu+1) with principal branches
        If (nu==zero) Then
          log_in = czero
        Else If (nu_pls_one>=zero) Then
          log_in = nu*(lg_half+log(z_in)) - log_gamma(nu_pls_one)
        Else
          log_in = nu*(lg_half+log(z_in)) - (log_gamma(nu_pls_one)+j1*merge(zero,pi,sign(one,sin(pi*nu))<zero))
        End If

        re_lg_in = real(log_in, kind=rk)
        im_lg_in = aimag(log_in)

! underflow guards as in your original
        If (re_lg_in<log_rmin) Then
          ierr = -1
          Return
        Else If (im_lg_in<log_rmin) Then
          im_lg_in = zero
          re_lg_in = real(log_in, kind=rk)
        Else
          re_lg_in = real(log_in, kind=rk)
          im_lg_in = aimag(log_in)
        End If

        exp_tmp = exp(re_lg_in)
        angle_check = mod(im_lg_in, two_pi)

        halfz_sq = quarter*z_in*z_in
! ---------- single place where I vs J differs:
        halfz_sq_eff = merge(halfz_sq, -halfz_sq, kind_flag==+1) ! I: +, J: -

        halfz_sq2 = real(halfz_sq,rk)**2+aimag(halfz_sq)**2 ! magnitude same for +/-; keep your original
        den_rcp = one/nu_pls_one
        eps_den_rcp = eps2*den_rcp*den_rcp
        atol = halfz_sq2*eps_den_rcp
! series accumulation (identical for both)
        t_k = cone
        cval = cone

        If (halfz_sq2>=eps_den_rcp) Then
          Do kk = 1, 100
            den_rcp = one/(kk*(kk+nu))
            t_k = t_k*halfz_sq_eff*den_rcp
            cval = cval + t_k
            If (real(t_k,rk)*real(t_k,rk)+aimag(t_k)*aimag(t_k)<=atol) Exit
          End Do
        End If

        If (nu/=zero) Then
          cval = cval*cmplx(cos(angle_check), sin(angle_check), kind=rk)*exp_tmp
        End If
        
      End Subroutine bessel_series_core
      
      
      
!*******************************************************************************
! SUBROUTINE: unif_sum_core
!
! PURPOSE
! -------
! Shared large-order uniform asymptotic expansion for I_nu(z) and K_nu(z).
!
! INPUTS
! ------
!   nu        : real(rk), positive order used in the expansion
!   z_in      : complex(rk), argument
!   kind_flag : integer
!               +1 selects the I_nu prefactor and non-alternating series
!               -1 selects the K_nu prefactor and alternating series
!
! OUTPUTS
! -------
!   cval : complex(rk), computed value
!   ierr : integer, status flag
!
! METHOD
! ------
! The routine uses Debye-type variables p and eta and the coefficient table
! c_uni_asy. The same coefficient table is used for I and K; the exponential
! prefactor and signs are selected using kind_flag.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: c_uni_asy, lg_one_sqrt_2_pi, lg_sqrt_half_pi,
!     eps2, cone, czero, one, half, zero, log_rmax, and log_rmin.
!
! Fortran intrinsic procedures used:
!   - aimag, exp, log, real, sqrt.
!*******************************************************************************

      Subroutine unif_sum_core(nu, z_in, kind_flag, cval, ierr)
        Implicit None
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Integer, Intent (In) :: kind_flag ! +1: I,  -1: K
        Complex (rk), Intent (Out) :: cval
        Integer, Intent (Out) :: ierr

        Complex (rk) :: z_sqr, p, p_sqr, taw_k, eta, lg_pref
        Complex (rk) :: taw_old, sum_main, sum_k, increment
        Complex (rk) :: p_sqr_kk, p_sqr_j, inv_psqr, power
        Real (rk) :: nu_inv
        Integer :: kk, j, m, ksign

        ierr = 0
        cval = czero
        If (nu==0.0_rk) Then
          ierr = -2
          Return
        End If

        nu_inv = one/nu
        z_sqr = (z_in*nu_inv)
        z_sqr = z_sqr*z_sqr
        p = cone/sqrt(cone+z_sqr)
        p_sqr = p*p
        taw_k = p*nu_inv
        eta = cone/p + log(taw_k*z_in/(cone+p))

        If (kind_flag==+1) Then
          lg_pref = lg_one_sqrt_2_pi + nu*eta + half*log(taw_k) ! I_nu
        Else
          lg_pref = lg_sqrt_half_pi - nu*eta + half*log(taw_k) ! K_nu
        End If

        If (real(lg_pref,kind=rk)>log_rmax) Then
          ierr = 1
          cval = czero
          Return
        Else If (real(lg_pref,kind=rk)<log_rmin) Then
          ierr = -1
          cval = czero
          Return
        End If

        sum_main = cone
        taw_old = cone
        m = 1

        If (kind_flag==+1) Then
! -------- I_nu branch (non-alternating) --------
          Do kk = 1, 28
            taw_old = taw_old*taw_k
            sum_k = czero

            power = p_sqr**kk ! p_sqr^(kk), then down by /p_sqr each j
            Do j = 0, kk
              m = m + 1
              sum_k = sum_k + c_uni_asy(m)*power
              power = power/p_sqr
            End Do

            increment = taw_old*sum_k
            sum_main = sum_main + increment
            If (real(increment,kind=rk)**2+aimag(increment)**2<eps2) Exit
          End Do

        Else
! -------- K_nu branch (alternating, CORRECTED POWERS) --------
          ksign = +1
          p_sqr_kk = cone
          inv_psqr = cone/p_sqr

          Do kk = 1, 28
            taw_old = taw_old*taw_k
            sum_k = czero
            ksign = -ksign
            p_sqr_kk = p_sqr_kk*p_sqr ! now p_sqr_kk = p_sqr**kk

            p_sqr_j = cone ! so term uses p_sqr_kk * p_sqr_j = p_sqr**(kk-j)
            Do j = 0, kk
              m = m + 1
              sum_k = sum_k + c_uni_asy(m)*p_sqr_kk*p_sqr_j
              p_sqr_j = p_sqr_j*inv_psqr ! -> p_sqr^{-j}
            End Do

            increment = ksign*taw_old*sum_k
            sum_main = sum_main + increment
            If (real(increment, rk)**2+aimag(increment)**2<eps2) Exit
          End Do
        End If

        cval = exp(lg_pref)*sum_main
        
      End Subroutine unif_sum_core





!==============================================================================
!          5. Deprecated or validation routines, can be removed
!==============================================================================


!*******************************************************************************
! SUBROUTINE: k_nu_unif_sum
!
! PURPOSE
! -------
! Alternate large-order uniform asymptotic expansion for K_nu(z).
!
! NOTE
! ----
! This routine is retained as an alternate K_nu implementation. In the current
! dispatcher, the shared routine unif_sum_core is used for the main large-order
! K_nu path.
!
! DEPENDENCIES
! ------------
! Direct internal calls:
!   - None.
!
! Current use in this module:
!   - This routine is retained as a private K_nu uniform-asymptotic helper.
!     The active dispatcher currently calls the unified routine unif_sum_core
!     instead.
!
! Module-level constants and data used:
!   - From set_rk: rk.
!   - From parameters.f90: c_uni_asy, lg_sqrt_half_pi, eps2, cone, one,
!     half, zero, log_rmax, and log_rmin.
!
! Fortran intrinsic procedures used:
!   - aimag, exp, log, real, sqrt.
!*******************************************************************************

      Subroutine k_nu_unif_sum(nu, z_in, cbk, ierr)

! Input and output variables
        Real (rk), Intent (In) :: nu
        Complex (rk), Intent (In) :: z_in
        Complex (rk), Intent (Out) :: cbk
        Integer, Intent (Out) :: ierr
!
! Local variables
        Complex (rk) :: k_nu_sum, k_nu_sum_2, z_sqr, p, taw_k, eta, increment
        Complex (rk) :: taw_old, p_sqr, lg_k_nu, inv_psqr, p_sqr_kk, p_sqr_j
        Real (rk) :: nu_inv
        Integer :: kk, j, m, ksign

        ierr = 0
        nu_inv = one/nu

! Compute core values
        z_sqr = (z_in*nu_inv)
        z_sqr = z_sqr*z_sqr

        p = cone/sqrt(cone+z_sqr)
        p_sqr = p*p
        taw_k = p*nu_inv
        eta = cone/p + log(taw_k*z_in/(cone+p))

! Compute logarithm of prefactor
        lg_k_nu = lg_sqrt_half_pi - nu*eta + half*log(taw_k)

! Overflow check
        If (real(lg_k_nu,kind=rk)>log_rmax) Then
          ierr = 1 ! Overflow detected
          cbk = cmplx(ieee_value(1.0_rk,ieee_positive_inf), ieee_value(1.0_rk,ieee_positive_inf), kind=rk)
          Return
        Else If (real(lg_k_nu,kind=rk)<log_rmin) Then
          ierr = -1 !underflow detected
          cbk = czero
          Return
        End If

! Initialize summation
        k_nu_sum = cone ! Start with 1 in complex form
        taw_old = cone
        m = 1
        ksign = 1
        p_sqr_kk = cone
        inv_psqr = cone/p_sqr
! Summation loop
        Do kk = 1, 28
          taw_old = taw_old*taw_k
          k_nu_sum_2 = czero
          ksign = -ksign
          p_sqr_kk = p_sqr_kk*p_sqr

          p_sqr_j = p_sqr
          Do j = 0, kk
            m = m + 1
            p_sqr_j = p_sqr_j*inv_psqr
!            k_nu_sum_2 = k_nu_sum_2 + c_uni_asy(m)*(p_sqr**(kk-j))
            k_nu_sum_2 = k_nu_sum_2 + c_uni_asy(m)*p_sqr_kk*p_sqr_j !
          End Do

          increment = ksign*taw_old*k_nu_sum_2

! Update main sum
          k_nu_sum = k_nu_sum + increment

! Convergence check
          If (real(increment,kind=rk)**2+aimag(increment)**2<eps2) Exit
        End Do

        cbk = exp(lg_k_nu+log(k_nu_sum))
        Return

      End Subroutine  k_nu_unif_sum




    End Module Bessel_Cmplx

