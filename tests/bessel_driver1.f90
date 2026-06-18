!----------------------------------------------------------------------
    Program bessel_driver1
!----------------------------------------------------------------------
! Driver for verifying & timing I_nu(z), J_nu(z), K_nu(z), Y_nu(z)
! in double/quad precision.
!
! - Uses present implementation:
!       i_nu_of_z, j_nu_of_z, k_nu_of_z, y_nu_of_z
! - Optional compare against AMOS Algorithm 644:
!       cbesi, cbesj, cbesk, cbesy
! - Optional compare against Algorithm 912:
!       bessel1, bessel2, hankel1, hankel2
!
! Comparison options:
!       cmp_none : present algorithm only
!       cmp_644  : compare with Algorithm 644 only
!       cmp_912  : compare with Algorithm 912 only
!       cmp_both : compare with both 644 and 912
!----------------------------------------------------------------------

      Use set_rk
      Use, Intrinsic :: ieee_exceptions
      Use, Intrinsic :: ieee_arithmetic
      Use complex_bessel, Only: cbesi, cbesj, cbesk, cbesy
      Use Bessel_Cmplx, Only: i_nu_of_z, j_nu_of_z, k_nu_of_z, y_nu_of_z
      Use mod_zbes, Only: bessel1, bessel2, hankel1, hankel2

      Implicit None
      Include 'parameters.f90'

!=============================================
!==== CONTROL ================================
!=============================================
      Integer, Parameter :: case_i_dp_pos = 1
      Integer, Parameter :: case_i_dp_neg = 2
      Integer, Parameter :: case_i_qp_pos = 3
      Integer, Parameter :: case_i_qp_neg = 4

      Integer, Parameter :: case_j_dp_pos = 5
      Integer, Parameter :: case_j_dp_neg = 6
      Integer, Parameter :: case_j_qp_pos = 7
      Integer, Parameter :: case_j_qp_neg = 8

      Integer, Parameter :: case_k_dp_pos = 9
      Integer, Parameter :: case_k_dp_neg = 10
      Integer, Parameter :: case_k_qp_pos = 11
      Integer, Parameter :: case_k_qp_neg = 12

      Integer, Parameter :: case_y_dp_pos = 13
      Integer, Parameter :: case_y_dp_neg = 14
      Integer, Parameter :: case_y_qp_pos = 15
      Integer, Parameter :: case_y_qp_neg = 16

      Integer, Parameter :: run_case = 6

      Integer, Parameter :: cmp_none = 0
      Integer, Parameter :: cmp_644  = 1
      Integer, Parameter :: cmp_912  = 2
      Integer, Parameter :: cmp_both = 3

      Integer, Parameter :: comparison_mode = cmp_912

      Logical, Parameter :: do_per_point_timing = .True.

! Derived flags
      Logical, Parameter :: is_i_case = &
        run_case==case_i_dp_pos .Or. run_case==case_i_dp_neg .Or. &
        run_case==case_i_qp_pos .Or. run_case==case_i_qp_neg

      Logical, Parameter :: is_j_case = &
        run_case==case_j_dp_pos .Or. run_case==case_j_dp_neg .Or. &
        run_case==case_j_qp_pos .Or. run_case==case_j_qp_neg

      Logical, Parameter :: is_k_case = &
        run_case==case_k_dp_pos .Or. run_case==case_k_dp_neg .Or. &
        run_case==case_k_qp_pos .Or. run_case==case_k_qp_neg

      Logical, Parameter :: is_y_case = &
        run_case==case_y_dp_pos .Or. run_case==case_y_dp_neg .Or. &
        run_case==case_y_qp_pos .Or. run_case==case_y_qp_neg

      Logical, Parameter :: is_dp_case = &
        run_case==case_i_dp_pos .Or. run_case==case_i_dp_neg .Or. &
        run_case==case_j_dp_pos .Or. run_case==case_j_dp_neg .Or. &
        run_case==case_k_dp_pos .Or. run_case==case_k_dp_neg .Or. &
        run_case==case_y_dp_pos .Or. run_case==case_y_dp_neg

      Logical, Parameter :: is_qp_case = &
        run_case==case_i_qp_pos .Or. run_case==case_i_qp_neg .Or. &
        run_case==case_j_qp_pos .Or. run_case==case_j_qp_neg .Or. &
        run_case==case_k_qp_pos .Or. run_case==case_k_qp_neg .Or. &
        run_case==case_y_qp_pos .Or. run_case==case_y_qp_neg

      Logical, Parameter :: is_neg_case = &
        run_case==case_i_dp_neg .Or. run_case==case_i_qp_neg .Or. &
        run_case==case_j_dp_neg .Or. run_case==case_j_qp_neg .Or. &
        run_case==case_k_dp_neg .Or. run_case==case_k_qp_neg .Or. &
        run_case==case_y_dp_neg .Or. run_case==case_y_qp_neg

! Base toggles
      Logical, Parameter :: do_accuracy_base = .True.
      Logical, Parameter :: do_timing_base = .True.

      Logical, Parameter :: do_accuracy = do_accuracy_base
      Logical, Parameter :: do_timing = do_timing_base

! Algorithm 644: DP only; timing only for nonnegative nu
      Logical, Parameter :: do_compare644 = &
        (comparison_mode==cmp_644 .Or. comparison_mode==cmp_both) .And. is_dp_case

      Logical, Parameter :: do_time644 = do_compare644 .And. .Not. is_neg_case

! Algorithm 912: enabled for both DP and QP
      Logical, Parameter :: do_compare912 = &
        (comparison_mode==cmp_912 .Or. comparison_mode==cmp_both)

      Logical, Parameter :: do_time912 = do_compare912

! Optional region filter
      Logical, Parameter :: apply_region_filter = .False.

! Batch timing loops
      Integer, Parameter :: n_repeat = (qp/rk) + (1-rk_by_qp)*4  !9
      Integer, Parameter :: outer_timing_loops = (qp/rk) + (1-rk_by_qp)*4 !48

! Per-point timing controls
      Integer, Parameter :: n_pt_iters = 100 !1000
      Integer, Parameter :: pt_stride = 1

!==========================
!==== FILES ===============
!==========================
      Character (Len=*), Parameter :: f_i_dp_pos = 'data/small_dp_I_Maple.txt'
      Character (Len=*), Parameter :: f_i_dp_neg = 'data/small_dp_I_Maple_negnu.txt'
      Character (Len=*), Parameter :: f_i_qp_pos = 'data/small_qp_I_Maple.txt'
      Character (Len=*), Parameter :: f_i_qp_neg = 'data/small_qp_I_Maple_negnu.txt'

      Character (Len=*), Parameter :: f_j_dp_pos = 'data/small_dp_J_Maple.txt'
      Character (Len=*), Parameter :: f_j_dp_neg = 'data/small_dp_J_Maple_negnu.txt'
      Character (Len=*), Parameter :: f_j_qp_pos = 'data/small_qp_J_Maple.txt'
      Character (Len=*), Parameter :: f_j_qp_neg = 'data/small_qp_J_Maple_negnu.txt'

      Character (Len=*), Parameter :: f_k_dp_pos = 'data/small_dp_K_Maple.txt'
      Character (Len=*), Parameter :: f_k_dp_neg = 'data/small_dp_K_Maple_negnu.txt'
      Character (Len=*), Parameter :: f_k_qp_pos = 'data/small_qp_K_Maple.txt'
      Character (Len=*), Parameter :: f_k_qp_neg = 'data/small_qp_K_Maple_negnu.txt'

      Character (Len=*), Parameter :: f_y_dp_pos = 'data/small_dp_Y_Maple.txt'
      Character (Len=*), Parameter :: f_y_dp_neg = 'data/small_dp_Y_Maple_negnu.txt'
      Character (Len=*), Parameter :: f_y_qp_pos = 'data/small_qp_Y_Maple.txt'
      Character (Len=*), Parameter :: f_y_qp_neg = 'data/small_qp_Y_Maple_negnu.txt'

      Character (Len=*), Parameter :: acc_out_dp_i = 'dp_acc_cmp_I.txt'
      Character (Len=*), Parameter :: acc_out_qp_i = 'qp_acc_cmp_I.txt'
      Character (Len=*), Parameter :: acc_out_dp_j = 'dp_acc_cmp_J.txt'
      Character (Len=*), Parameter :: acc_out_qp_j = 'qp_acc_cmp_J.txt'
      Character (Len=*), Parameter :: acc_out_dp_k = 'dp_acc_cmp_K.txt'
      Character (Len=*), Parameter :: acc_out_qp_k = 'qp_acc_cmp_K.txt'
      Character (Len=*), Parameter :: acc_out_dp_y = 'dp_acc_cmp_Y.txt'
      Character (Len=*), Parameter :: acc_out_qp_y = 'qp_acc_cmp_Y.txt'

      Character (Len=*), Parameter :: time_out_i = 'timing_summary_I.txt'
      Character (Len=*), Parameter :: time_out_j = 'timing_summary_J.txt'
      Character (Len=*), Parameter :: time_out_k = 'timing_summary_K.txt'
      Character (Len=*), Parameter :: time_out_y = 'timing_summary_Y.txt'

      Character (Len=64) :: acc_file, time_file
      Character (Len=1) :: func_tag

!==========================
!==== DATA / STATE ========
!==========================
      Integer, Parameter :: nmax_in = 450000
      Integer :: nmin, nmax, jn, jkeep, total_read, total_kept
      Complex (rk), Allocatable :: z1(:)
      Real (rk), Allocatable :: fnu(:), xt(:), yt(:), ref_re(:), ref_im(:)

      Integer :: unit_in, unit_acc, unit_time, ios
      Integer :: nnz, ierr, ierr644, ierr912
      Complex (rk) :: cval, yloc(1), z_in, zans, val912

      Real (rk) :: err_re, err_im, err_re_max, err_im_max, nu

      Real (rk) :: err_re_644, err_im_644
      Real (rk) :: err_re_644_max, err_im_644_max
      Integer :: jn_err_re_max, jn_err_im_max
      Integer :: jn_err_re_644_max, jn_err_im_644_max

      Real (rk) :: err_re_912, err_im_912
      Real (rk) :: err_re_912_max, err_im_912_max
      Integer :: jn_err_re_912_max, jn_err_im_912_max

      Integer (selected_int_kind(18)) :: cs, cf, r, cm
      Real (rk) :: best_present, best_644, best_912, dt
      Real (rk) :: elapsed_pt, elapsed_644_pt, elapsed_912_pt
      Real (rk) :: elapsed_pt_sum, elapsed_644_sum, elapsed_912_sum
      Logical :: sampled_this_point
      Integer :: unit_ptime, s

!---------------------------
! Optional quick single-point check
!---------------------------
      Write (*,*) 'Optional quick single-point check'
      nu = -0.1000000000000000E+00_rk
      z_in = 1.4424172552410660E+01_rk + j1*1.2655308231907409E-06_rk

      Write (*,*) 'nu=',  nu, 'z=', z_in

      Call j_nu_of_z(nu,z_in,cval,ierr)
      Write (*,*) 'j_nu_of_z=', cval

      Call i_nu_of_z(nu,z_in,cval,ierr)
      Write (*,*) 'i_nu_of_z=', cval

      Call k_nu_of_z(nu,-j1*z_in,cval,ierr)
      Write (*,*) 'k_nu_of_z=', cval

      Call cbesj(z_in,nu,1,1,yloc(1),nnz,ierr644)
      Write (*,*) 'j_644=', yloc(1)

      Call cbesi(z_in,nu,1,1,yloc(1),nnz,ierr644)
      Write (*,*) 'i_644=', yloc(1)

      Call cbesk(-j1*z_in,nu,1,1,yloc(1),nnz,ierr644)
      Write (*,*) 'k_644=', yloc(1)

      Call bessel1(cmplx(nu,zero,rk),z_in,zans,ierr912)
      Write (*,*) 'j_912=', zans

      Call y_nu_of_z(nu,z_in,cval,ierr)
      Write (*,*) 'y_nu_of_z=', cval

      Call cbesy(z_in,nu,1,1,yloc(1),nnz,ierr644)
      Write (*,*) 'y_644=', yloc(1)

      Call bessel2(cmplx(nu,zero,rk),z_in,zans,ierr912)
      Write (*,*) 'y_912=', zans

!---------------------------
! Consistency: case vs. rk
!---------------------------
      If (is_dp_case .And. rk/=dp) Then
        Write (*,*) 'ERROR: case requires double precision, but rk /= dp.'
        Stop 1
      End If

      If (is_qp_case .And. rk/=qp) Then
        Write (*,*) 'ERROR: case requires quad precision, but rk /= qp.'
        Stop 1
      End If

!==========================
!==== RUNTIME GUARDS ======
!==========================
      Call ieee_set_halting_mode(ieee_underflow, .False.)
      Call ieee_set_halting_mode(ieee_invalid, .False.)

!==========================
!==== OPEN INPUT FILE =====
!==========================
      Select Case (run_case)

      Case (case_i_dp_pos)
        Open (Newunit=unit_in, File=f_i_dp_pos, Status='old', Action='read')
      Case (case_i_dp_neg)
        Open (Newunit=unit_in, File=f_i_dp_neg, Status='old', Action='read')
      Case (case_i_qp_pos)
        Open (Newunit=unit_in, File=f_i_qp_pos, Status='old', Action='read')
      Case (case_i_qp_neg)
        Open (Newunit=unit_in, File=f_i_qp_neg, Status='old', Action='read')

      Case (case_j_dp_pos)
        Open (Newunit=unit_in, File=f_j_dp_pos, Status='old', Action='read')
      Case (case_j_dp_neg)
        Open (Newunit=unit_in, File=f_j_dp_neg, Status='old', Action='read')
      Case (case_j_qp_pos)
        Open (Newunit=unit_in, File=f_j_qp_pos, Status='old', Action='read')
      Case (case_j_qp_neg)
        Open (Newunit=unit_in, File=f_j_qp_neg, Status='old', Action='read')

      Case (case_k_dp_pos)
        Open (Newunit=unit_in, File=f_k_dp_pos, Status='old', Action='read')
      Case (case_k_dp_neg)
        Open (Newunit=unit_in, File=f_k_dp_neg, Status='old', Action='read')
      Case (case_k_qp_pos)
        Open (Newunit=unit_in, File=f_k_qp_pos, Status='old', Action='read')
      Case (case_k_qp_neg)
        Open (Newunit=unit_in, File=f_k_qp_neg, Status='old', Action='read')

      Case (case_y_dp_pos)
        Open (Newunit=unit_in, File=f_y_dp_pos, Status='old', Action='read')
      Case (case_y_dp_neg)
        Open (Newunit=unit_in, File=f_y_dp_neg, Status='old', Action='read')
      Case (case_y_qp_pos)
        Open (Newunit=unit_in, File=f_y_qp_pos, Status='old', Action='read')
      Case (case_y_qp_neg)
        Open (Newunit=unit_in, File=f_y_qp_neg, Status='old', Action='read')

      Case Default
        Write (*,*) 'Invalid run_case.'
        Stop 1
      End Select

! Select output filenames + function tag
      If (is_i_case) Then
        acc_file = merge(acc_out_qp_i, acc_out_dp_i, rk==qp)
        time_file = time_out_i
        func_tag = 'I'
      Else If (is_j_case) Then
        acc_file = merge(acc_out_qp_j, acc_out_dp_j, rk==qp)
        time_file = time_out_j
        func_tag = 'J'
      Else If (is_k_case) Then
        acc_file = merge(acc_out_qp_k, acc_out_dp_k, rk==qp)
        time_file = time_out_k
        func_tag = 'K'
      Else
        acc_file = merge(acc_out_qp_y, acc_out_dp_y, rk==qp)
        time_file = time_out_y
        func_tag = 'Y'
      End If

!==========================
!==== READ & FILTER =======
!==========================
      Allocate (z1(nmax_in), fnu(nmax_in), xt(nmax_in), yt(nmax_in), &
                ref_re(nmax_in), ref_im(nmax_in))

      jkeep = 0
      total_read = 0
      total_kept = 0

      Do
        If (read_line(unit_in,fnu(jkeep+1),xt(jkeep+1),yt(jkeep+1), &
                      ref_re(jkeep+1),ref_im(jkeep+1))) Exit

        total_read = total_read + 1

        If (.Not. apply_region_filter .Or. &
            in_region(fnu(jkeep+1),xt(jkeep+1),yt(jkeep+1))) Then
          jkeep = jkeep + 1
          z1(jkeep) = cmplx(xt(jkeep),yt(jkeep),rk)
          total_kept = total_kept + 1
        End If

        If (jkeep==nmax_in) Exit
      End Do

      Close (unit_in)

      nmin = 1
      nmax = jkeep

      Write (*,'(A,I0,A,I0,A,I0)') 'Read: ', total_read, &
        '  kept: ', total_kept, '  dropped: ', total_read-total_kept

      If (nmax==0) Then
        Write (*,*) 'ERROR: No points loaded.'
        Stop 1
      End If

      Write (*,'(A,I0)') 'Loaded test points: ', nmax

!==========================
!==== ACCURACY PASS =======
!==========================
      If (do_accuracy) Then

        Open (Newunit=unit_acc, File=acc_file, Status='replace', Action='write')

        If (do_per_point_timing) Then
          Open (Newunit=unit_ptime, &
                File=trim(adjustl('per_point_times_'//func_tag))//'.txt', &
                Status='replace', Action='write')
          Call system_clock(cs,r,cm)
        End If

        err_re_max = 0.0_rk
        err_im_max = 0.0_rk
        err_re_644_max = 0.0_rk
        err_im_644_max = 0.0_rk
        err_re_912_max = 0.0_rk
        err_im_912_max = 0.0_rk

        jn_err_re_max = 1
        jn_err_im_max = 1
        jn_err_re_644_max = 1
        jn_err_im_644_max = 1
        jn_err_re_912_max = 1
        jn_err_im_912_max = 1

        Do jn = nmin, nmax

          sampled_this_point = do_per_point_timing .And. mod(jn,pt_stride)==0

          ierr = 0

          If (sampled_this_point) Then
            Call eval_one_present(fnu(jn),z1(jn),cval,ierr)

            elapsed_pt_sum = 0.0_rk
            Do s = 1, n_pt_iters
              Call system_clock(cs,r,cm)
              Call eval_one_present(fnu(jn),z1(jn),cval,ierr)
              Call system_clock(cf,r,cm)
              elapsed_pt_sum = elapsed_pt_sum + real(cf-cs,rk)/real(r,rk)
            End Do

            elapsed_pt = elapsed_pt_sum/real(n_pt_iters,rk)
          Else
            Call eval_one_present(fnu(jn),z1(jn),cval,ierr)
            elapsed_pt = 0.0_rk
          End If

          If (ierr==0 .And. ieee_is_finite(real(cval,rk)) .And. &
              ieee_is_finite(aimag(cval))) Then

            err_re = relerr(real(cval,rk),ref_re(jn))
            err_im = relerr(aimag(cval),ref_im(jn))

            If (err_re>err_re_max) Then
              err_re_max = err_re
              jn_err_re_max = jn
            End If

            If (err_im>err_im_max) Then
              err_im_max = err_im
              jn_err_im_max = jn
            End If

            yloc(1) = cmplx(zero,zero,rk)
            val912 = cmplx(zero,zero,rk)

            err_re_644 = 0.0_rk
            err_im_644 = 0.0_rk
            err_re_912 = 0.0_rk
            err_im_912 = 0.0_rk

! Algorithm 644
            If (do_compare644) Then
              ierr644 = 0

              If (sampled_this_point .And. do_time644) Then
                Call eval_one_644(fnu(jn),z1(jn),yloc(1),ierr644)

                elapsed_644_sum = 0.0_rk
                Do s = 1, n_pt_iters
                  Call system_clock(cs,r,cm)
                  Call eval_one_644(fnu(jn),z1(jn),yloc(1),ierr644)
                  Call system_clock(cf,r,cm)
                  elapsed_644_sum = elapsed_644_sum + real(cf-cs,rk)/real(r,rk)
                End Do

                elapsed_644_pt = elapsed_644_sum/real(n_pt_iters,rk)
              Else
                Call eval_one_644(fnu(jn),z1(jn),yloc(1),ierr644)
                elapsed_644_pt = 0.0_rk
              End If

              err_re_644 = relerr(real(yloc(1),rk),ref_re(jn))
              err_im_644 = relerr(aimag(yloc(1)),ref_im(jn))

              If (err_re_644>err_re_644_max) Then
                err_re_644_max = err_re_644
                jn_err_re_644_max = jn
              End If

              If (err_im_644>err_im_644_max) Then
                err_im_644_max = err_im_644
                jn_err_im_644_max = jn
              End If
            End If

! Algorithm 912
            If (do_compare912) Then
              ierr912 = 0

              If (sampled_this_point .And. do_time912) Then
                Call eval_one_912(fnu(jn),z1(jn),val912,ierr912)

                elapsed_912_sum = 0.0_rk
                Do s = 1, n_pt_iters
                  Call system_clock(cs,r,cm)
                  Call eval_one_912(fnu(jn),z1(jn),val912,ierr912)
                  Call system_clock(cf,r,cm)
                  elapsed_912_sum = elapsed_912_sum + real(cf-cs,rk)/real(r,rk)
                End Do

                elapsed_912_pt = elapsed_912_sum/real(n_pt_iters,rk)
              Else
                Call eval_one_912(fnu(jn),z1(jn),val912,ierr912)
                elapsed_912_pt = 0.0_rk
              End If

              err_re_912 = relerr(real(val912,rk),ref_re(jn))
              err_im_912 = relerr(aimag(val912),ref_im(jn))

              If (err_re_912>err_re_912_max) Then
                err_re_912_max = err_re_912
                jn_err_re_912_max = jn
              End If

              If (err_im_912>err_im_912_max) Then
                err_im_912_max = err_im_912
                jn_err_im_912_max = jn
              End If
            End If

! Accuracy output
            If (comparison_mode==cmp_none) Then

              Write (unit_acc,'(I9,2X,ES26.16E4,2X,2(ES26.16E4,2X), &
                   2(ES13.4E4,2X),4(ES26.16E4,2X))') &
                jn, fnu(jn), real(z1(jn),rk), aimag(z1(jn)), &
                err_re, err_im, &
                ref_re(jn), ref_im(jn), &
                real(cval,rk), aimag(cval)

            Else If (comparison_mode==cmp_644) Then

              Write (unit_acc,'(I9,2X,ES26.16E4,2X,2(ES26.16E4,2X), &
                   4(ES13.4E4,2X),6(ES26.16E4,2X))') &
                jn, fnu(jn), real(z1(jn),rk), aimag(z1(jn)), &
                err_re, err_im, err_re_644, err_im_644, &
                ref_re(jn), ref_im(jn), &
                real(cval,rk), aimag(cval), &
                real(yloc(1),rk), aimag(yloc(1))

            Else If (comparison_mode==cmp_912) Then

              Write (unit_acc,'(I9,2X,ES26.16E4,2X,2(ES26.16E4,2X), &
                   4(ES13.4E4,2X),6(ES26.16E4,2X))') &
                jn, fnu(jn), real(z1(jn),rk), aimag(z1(jn)), &
                err_re, err_im, err_re_912, err_im_912, &
                ref_re(jn), ref_im(jn), &
                real(cval,rk), aimag(cval), &
                real(val912,rk), aimag(val912)

            Else If (comparison_mode==cmp_both) Then

              Write (unit_acc,'(I9,2X,ES26.16E4,2X,2(ES26.16E4,2X), &
                   6(ES13.4E4,2X),8(ES26.16E4,2X))') &
                jn, fnu(jn), real(z1(jn),rk), aimag(z1(jn)), &
                err_re, err_im, err_re_644, err_im_644, &
                err_re_912, err_im_912, &
                ref_re(jn), ref_im(jn), &
                real(cval,rk), aimag(cval), &
                real(yloc(1),rk), aimag(yloc(1)), &
                real(val912,rk), aimag(val912)

            End If

            If (do_per_point_timing .And. sampled_this_point) Then
              Write (unit_ptime,'(I9,1X,ES26.16E4,1X,2(ES26.16E4,1X), &
                   3(ES13.4E4,1X))') &
                jn, fnu(jn), real(z1(jn),rk), aimag(z1(jn)), &
                elapsed_pt, elapsed_644_pt, elapsed_912_pt
            End If

          End If
        End Do

        If (do_per_point_timing) Close (unit_ptime)
        Close (unit_acc)

!==========================
!==== ACCURACY SUMMARY ====
!==========================
        Write (*,'(/A)') repeat('=',70)
        Write (*,'(A)') 'CROSS-COMPARISON OF MAXIMUM RELATIVE ERROR POINTS'
        Write (*,'(A)') repeat('=',70)

        Call report_errors_at_point( &
          'Point where Present algorithm has maximum RelErr(Re)', &
          jn_err_re_max)

        Call report_errors_at_point( &
          'Point where Present algorithm has maximum RelErr(Im)', &
          jn_err_im_max)

        If (do_compare644) Then
          Call report_errors_at_point( &
            'Point where Algorithm 644 has maximum RelErr(Re)', &
            jn_err_re_644_max)

          Call report_errors_at_point( &
            'Point where Algorithm 644 has maximum RelErr(Im)', &
            jn_err_im_644_max)
        End If

        If (do_compare912) Then
          Call report_errors_at_point( &
            'Point where Algorithm 912 has maximum RelErr(Re)', &
            jn_err_re_912_max)

          Call report_errors_at_point( &
            'Point where Algorithm 912 has maximum RelErr(Im)', &
            jn_err_im_912_max)
        End If

        Write (*,'(/A)') repeat('=',70)
        Write (*,'(A)') 'GLOBAL MAXIMUM RELATIVE ERRORS'
        Write (*,'(A)') repeat('=',70)

        Write (*,'(A,ES13.4E4,A,I0)') &
          'Present Max RelErr(Re) = ', err_re_max, ' at index ', jn_err_re_max
        Write (*,'(A,ES13.4E4,A,I0)') &
          'Present Max RelErr(Im) = ', err_im_max, ' at index ', jn_err_im_max

        If (do_compare644) Then
          Write (*,'(A,ES13.4E4,A,I0)') &
            'Alg. 644 Max RelErr(Re) = ', err_re_644_max, ' at index ', jn_err_re_644_max
          Write (*,'(A,ES13.4E4,A,I0)') &
            'Alg. 644 Max RelErr(Im) = ', err_im_644_max, ' at index ', jn_err_im_644_max
        End If

        If (do_compare912) Then
          Write (*,'(A,ES13.4E4,A,I0)') &
            'Alg. 912 Max RelErr(Re) = ', err_re_912_max, ' at index ', jn_err_re_912_max
          Write (*,'(A,ES13.4E4,A,I0)') &
            'Alg. 912 Max RelErr(Im) = ', err_im_912_max, ' at index ', jn_err_im_912_max
        End If

      End If

!==========================
!==== TIMING SWEEPS =======
!==========================
      If (do_timing) Then

        Open (Newunit=unit_time, File=time_file, Status='replace', Action='write')

        best_present = huge(1.0_rk)
        best_present = min(best_present,time_sweep_present())

        Write (*,'(/A,ES12.4)') 'Present min elapsed (s): ', best_present
        Write (unit_time,'(A,ES12.4)') 'present_min_elapsed = ', best_present

        If (do_time644) Then
          best_644 = huge(1.0_rk)
          best_644 = min(best_644,time_sweep_644())

          Write (*,'(A,ES12.4)') 'AMOS Algorithm 644 min elapsed (s): ', best_644
          Write (*,'(A,ES12.4)') 'ratio_present_over_644 = ', best_present/best_644

          Write (unit_time,'(A,ES12.4)') 'alg644_min_elapsed = ', best_644
          Write (unit_time,'(A,ES12.4)') 'ratio_present_over_644 = ', best_present/best_644
        End If

        If (do_time912) Then
          best_912 = huge(1.0_rk)
          best_912 = min(best_912,time_sweep_912())

          Write (*,'(A,ES12.4)') 'Algorithm 912 min elapsed (s): ', best_912
          Write (*,'(A,ES12.4)') 'ratio_present_over_912 = ', best_present/best_912

          Write (unit_time,'(A,ES12.4)') 'alg912_min_elapsed = ', best_912
          Write (unit_time,'(A,ES12.4)') 'ratio_present_over_912 = ', best_present/best_912
        End If

        Close (unit_time)
      End If

    Contains

!----------------------------------------------------------------------
      Logical Function read_line(u,f,xr,yi,rr,ri)
        Integer, Intent (In) :: u
        Real (rk), Intent (Out) :: f, xr, yi, rr, ri
        Integer :: ios

        Read (u,*,Iostat=ios) f, xr, yi, rr, ri
        read_line = ios/=0
      End Function

!----------------------------------------------------------------------
      Real (rk) Function relerr(val,ref)
        Real (rk), Intent (In) :: val, ref

        If (abs(ref)>rmin) Then
          relerr = abs(val-ref)/abs(ref)
        Else
          relerr = 0.0_rk
        End If
      End Function

!----------------------------------------------------------------------
      Subroutine report_errors_at_point(title,jidx)
        Character (Len=*), Intent (In) :: title
        Integer, Intent (In) :: jidx

        Complex (rk) :: vp, v644, v912
        Real (rk) :: ep_re, ep_im
        Real (rk) :: e644_re, e644_im
        Real (rk) :: e912_re, e912_im
        Integer :: iep, ie644, ie912

        vp = cmplx(zero,zero,rk)
        v644 = cmplx(zero,zero,rk)
        v912 = cmplx(zero,zero,rk)

        Call eval_one_present(fnu(jidx),z1(jidx),vp,iep)

        ep_re = relerr(real(vp,rk),ref_re(jidx))
        ep_im = relerr(aimag(vp),ref_im(jidx))

        Write (*,'(/A)') repeat('-',70)
        Write (*,'(A)') trim(title)
        Write (*,'(A,I0)') 'Index = ', jidx

        Write (*,'(A,ES26.16)') 'Order nu = ', fnu(jidx)
        Write (*,'(A,ES26.16,1X,ES26.16)') &
          'Argument z = (Re, Im) = ', real(z1(jidx),rk), aimag(z1(jidx))

        Write (*,'(A,ES26.16,1X,ES26.16)') &
          'Reference value = (Re, Im) = ', ref_re(jidx), ref_im(jidx)

        Write (*,'(/A)') 'Errors and values at this same point:'

        Write (*,'(A,2ES13.4E4)') &
          'Present  RelErr(Re,Im) = ', ep_re, ep_im
        Write (*,'(A,2ES26.16)') &
          'Present  value(Re,Im)  = ', real(vp,rk), aimag(vp)

        If (do_compare644) Then
          Call eval_one_644(fnu(jidx),z1(jidx),v644,ie644)

          e644_re = relerr(real(v644,rk),ref_re(jidx))
          e644_im = relerr(aimag(v644),ref_im(jidx))

          Write (*,'(A,2ES13.4E4)') &
            'Alg. 644 RelErr(Re,Im) = ', e644_re, e644_im
          Write (*,'(A,2ES26.16)') &
            'Alg. 644 value(Re,Im)  = ', real(v644,rk), aimag(v644)
        End If

        If (do_compare912) Then
          Call eval_one_912(fnu(jidx),z1(jidx),v912,ie912)

          e912_re = relerr(real(v912,rk),ref_re(jidx))
          e912_im = relerr(aimag(v912),ref_im(jidx))

          Write (*,'(A,2ES13.4E4)') &
            'Alg. 912 RelErr(Re,Im) = ', e912_re, e912_im
          Write (*,'(A,2ES26.16)') &
            'Alg. 912 value(Re,Im)  = ', real(v912,rk), aimag(v912)
        End If

        Write (*,'(A)') repeat('-',70)
      End Subroutine

!----------------------------------------------------------------------
      Logical Function in_region(f,xr,yi)
        Real (rk), Intent (In) :: f, xr, yi
        Real (rk) :: az2, nuabs
        Logical :: small_z_rgn, infit_z_rgn, infit_nu_rgn

        in_region = .False.

        nuabs = abs(f)
        az2 = xr*xr + yi*yi

        small_z_rgn = az2 < 324.0_rk + 8.0_rk*nuabs
        infit_nu_rgn = nuabs > c1 + abs(sqrt(az2))
        infit_z_rgn = az2 > max(324.0_rk + 8.0_rk*nuabs, abs_z_brdr1**2) .And. &
                      two*sqrt(az2) >= nuabs*nuabs

        If (small_z_rgn .Or. infit_z_rgn .Or. infit_nu_rgn) in_region = .True.
      End Function

!----------------------------------------------------------------------
      Subroutine eval_one_present(fn,z,val,ierr_local)
        Real (rk), Intent (In) :: fn
        Complex (rk), Intent (In) :: z
        Complex (rk), Intent (Out) :: val
        Integer, Intent (Out) :: ierr_local

        ierr_local = 0

        If (is_i_case) Then
          Call i_nu_of_z(fn,z,val,ierr_local)
        Else If (is_j_case) Then
          Call j_nu_of_z(fn,z,val,ierr_local)
        Else If (is_k_case) Then
          Call k_nu_of_z(fn,z,val,ierr_local)
        Else
          Call y_nu_of_z(fn,z,val,ierr_local)
        End If
      End Subroutine

!----------------------------------------------------------------------
      Subroutine eval_one_644(fn,z,val,ierr_local)
        Real (rk), Intent (In) :: fn
        Complex (rk), Intent (In) :: z
        Complex (rk), Intent (Out) :: val
        Integer, Intent (Out) :: ierr_local

        Integer :: nnz_l
        Complex (rk) :: tmp(1)

        ierr_local = 0

        If (is_i_case) Then
          Call cbesi(z,fn,1,1,tmp,nnz_l,ierr_local)
        Else If (is_j_case) Then
          Call cbesj(z,fn,1,1,tmp,nnz_l,ierr_local)
        Else If (is_k_case) Then
          Call cbesk(z,fn,1,1,tmp,nnz_l,ierr_local)
        Else
          Call cbesy(z,fn,1,1,tmp,nnz_l,ierr_local)
        End If

        val = tmp(1)
      End Subroutine

!----------------------------------------------------------------------
      Subroutine eval_one_912(fn,z,val,ierr_local)
        Real (rk), Intent (In) :: fn
        Complex (rk), Intent (In) :: z
        Complex (rk), Intent (Out) :: val
        Integer, Intent (Out) :: ierr_local

        Complex (rk) :: nu_c

        ierr_local = 0
        nu_c = cmplx(fn,zero,rk)

        If (is_j_case) Then
          Call bessel1(nu_c,z,val,ierr_local)

        Else If (is_y_case) Then
          Call bessel2(nu_c,z,val,ierr_local)

        Else If (is_i_case) Then
          Call bessel1(nu_c,j1*z,val,ierr_local)
          val = exp(-j1*pi*fn/two)*val

        Else If (is_k_case) Then
          Call hankel1(nu_c,j1*z,val,ierr_local)
          val = (pi*j1/two)*exp(j1*pi*fn/two)*val
        End If
      End Subroutine

!----------------------------------------------------------------------
      Real (rk) Function time_sweep_present()
        Integer :: k, i
        Real (rk) :: best

        best = huge(1.0_rk)

        Do k = 1, n_repeat
          Call system_clock(cs,r,cm)

          Do i = 1, outer_timing_loops
            Call eval_block_present()
          End Do

          Call system_clock(cf,r,cm)
          dt = real(cf-cs,rk)/real(r,rk)

          If (dt<best) best = dt
        End Do

        time_sweep_present = best
      End Function

!----------------------------------------------------------------------
      Subroutine eval_block_present()
        Integer :: j, ierr_local

        Do j = nmin, nmax
          ierr_local = 0
          Call eval_one_present(fnu(j),z1(j),cval,ierr_local)
        End Do
      End Subroutine

!----------------------------------------------------------------------
      Real (rk) Function time_sweep_644()
        Integer :: k, i
        Real (rk) :: best

        best = huge(1.0_rk)

        Do k = 1, n_repeat
          Call system_clock(cs,r,cm)

          Do i = 1, outer_timing_loops
            Call eval_block_644()
          End Do

          Call system_clock(cf,r,cm)
          dt = real(cf-cs,rk)/real(r,rk)

          If (dt<best) best = dt
        End Do

        time_sweep_644 = best
      End Function

!----------------------------------------------------------------------
      Subroutine eval_block_644()
        Integer :: j, ierr_local, nnz_l
        Complex (rk) :: ytmp(1)

        Do j = nmin, nmax
          ierr_local = 0

          If (is_i_case) Then
            Call cbesi(z1(j),fnu(j),1,1,ytmp,nnz_l,ierr_local)
          Else If (is_j_case) Then
            Call cbesj(z1(j),fnu(j),1,1,ytmp,nnz_l,ierr_local)
          Else If (is_k_case) Then
            Call cbesk(z1(j),fnu(j),1,1,ytmp,nnz_l,ierr_local)
          Else
            Call cbesy(z1(j),fnu(j),1,1,ytmp,nnz_l,ierr_local)
          End If
        End Do
      End Subroutine

!----------------------------------------------------------------------
      Real (rk) Function time_sweep_912()
        Integer :: k, i
        Real (rk) :: best

        best = huge(1.0_rk)

        Do k = 1, n_repeat
          Call system_clock(cs,r,cm)

          Do i = 1, outer_timing_loops
            Call eval_block_912()
          End Do

          Call system_clock(cf,r,cm)
          dt = real(cf-cs,rk)/real(r,rk)

          If (dt<best) best = dt
        End Do

        time_sweep_912 = best
      End Function

!----------------------------------------------------------------------
      Subroutine eval_block_912()
        Integer :: j, ierr_local
        Complex (rk) :: vtmp

        Do j = nmin, nmax
          ierr_local = 0
          Call eval_one_912(fnu(j),z1(j),vtmp,ierr_local)
        End Do
      End Subroutine

    End Program
