# Bessel_Cmplx

Fortran routines for computing Bessel and modified Bessel functions of complex argument and real order.

The package evaluates

- `I_nu(z)`: modified Bessel function of the first kind
- `J_nu(z)`: Bessel function of the first kind
- `K_nu(z)`: modified Bessel function of the second kind
- `Y_nu(z)`: Bessel function of the second kind

The order `nu` is real and may be positive or negative. The argument `z` is complex.

## Main routines

The main user-callable routines are provided by the module `Bessel_Cmplx`:

```fortran
Use Bessel_Cmplx, Only: i_nu_of_z, j_nu_of_z, k_nu_of_z, y_nu_of_z
```

| Routine                                | Computes |
|---                                          |---               |
| `i_nu_of_z(nu,z,cbi,ierr)`  | `I_nu(z)`   |
| `j_nu_of_z(nu,z,cbj,ierr)`  | `J_nu(z)`   |
| `k_nu_of_z(nu,z,cbk,ierr)` | `K_nu(z)` |
| `y_nu_of_z(nu,z,cby,ierr)` | `Y_nu(z)`  |

All routines use the same interface style:

Real(rk)     :: nu
Complex(rk)  :: z, value
Integer      :: ierr
Call i_nu_of_z(nu, z, value, ierr)



## Files

A typical source directory contains:

set_rk.f90          precision-control module
parameters.f90      constants and coefficient tables included by Bessel_Cmplx.f90
Bessel_Cmplx.f90    main Bessel-function module
bessel_driver1.f90  verification and timing driver
Makefile            optional build file

The verification driver can also compare against external/reference implementations. For that, the following files are used:

Cbessel_miller.f90
mod_zbes.f90

Large reference-data files are included in the repository for comparison and accuracy verification. See the section [Reference data](#reference-data).

## Precision

Precision is controlled in `set_rk.f90`:

Use
Integer, Parameter :: rk = dp              for double precision, or
and
Integer, Parameter :: rk = qp             for quadruple precision.

After changing `rk`, recompile the complete package.

## Compilation

### Option 1: direct compiler command

For a simple program that uses only the present package, compile `set_rk.f90`, `Bessel_Cmplx.f90`, and your driver program. The file `parameters.f90` is included by `Bessel_Cmplx.f90`, 
so it must be in the same directory or in the compiler include path.

Example:
gfortran -O3 set_rk.f90 Bessel_Cmplx.f90 bessel_driver.f90  -o bessel_driver


For the full verification/timing driver, including the comparison modules and external reference data, use for example:
gfortran -O3 set_rk.f90 Cbessel_miller.f90  mod_zbes.f90 Bessel_Cmplx.f90 bessel_driver1.f90  -o bessel_driver1 



### Option 2: using the provided Makefile

The repository also includes a `Makefile`. With GNU Fortran selected in the Makefile, build the driver by running:

make
or
make -f makefile

This creates the executable:
bessel_driver1


To clean object files, module files, and executables, run:

make clean

The supplied Makefile expects these source names:

set_rk.f90
Cbessel_miller.f90
mod_zbes.f90
Cbessel.f90
Bessel_Cmplx.f90
bessel_driver1.f90


If your driver file is named `bessel_driver.f90` instead of `bessel_driver1.f90`, either rename the file or edit the object list and dependency rule in the Makefile.

The Makefile includes commented compiler settings for NAG, Intel `ifort`, Intel `ifx`, and GNU Fortran. Uncomment the compiler block appropriate for your system.

Note for Linux/macOS users: if the `clean` rule uses the Windows command `del`, replace it with a Unix-style clean rule such as:

clean:
	$(RM) *.o *.mod $(TARGET)


## Error flags

| `ierr` | Meaning |
|---:|---|
| `0` | Successful computation |
| `1` | Overflow or no reliable computation |
| `-1` | Underflow; result treated as zero |
| `99` | Input error, for example `K_nu(0)` |
| `-2` | Internal dispatcher condition; not expected in normal user calls |

## Numerical methods

The module uses different methods in different regions of the order-argument plane, including:

- power series for small and intermediate `|z|`
- large-`|z|` asymptotic expansions
- large-order uniform asymptotic expansions
- backward recurrence for `I_nu(z)`
- forward recurrence for `K_nu(z)`
- connection and reflection formulas for negative orders and for `Y_nu(z)`

The complex logarithm and square root follow the principal branches used by the Fortran intrinsic functions. Analytic continuation across the left half-plane is handled internally using explicit phase factors.

## Verification driver

The driver `bessel_driver1.f90` is intended for accuracy checks and timing studies. It can compute the present implementation alone or compare it with optional reference implementations:

- AMOS Algorithm 644-style routines: `cbesi`, `cbesj`, `cbesk`, `cbesy`
- Algorithm 912-style routines: `bessel1`, `bessel2`, `hankel1`, `hankel2`

Inside the driver, select the case to run by changing:
Integer, Parameter :: run_case = 6


The driver includes cases for `I`, `J`, `K`, and `Y`, for double and quadruple precision, and for positive and negative orders.

Comparison mode is selected by:
Integer, Parameter :: comparison_mode = cmp_912


Available modes are:
cmp_none
cmp_644
cmp_912
cmp_both


## Reference data

The driver reads reference data files containing rows of the form:
nu   Re(z)   Im(z)   Re(reference)   Im(reference)

The reference data files are included in the repository arereduced from original reference data files by retaining one point of each 5 points. The 16 reduced reference data files are:
Reduced reference data files used by the verification driver
===========================================================

run_case   Function      Precision   Order sign     Reduced reference data file
--------   --------      ---------   ----------     ---------------------------
1          I_nu(z)       double      positive       small_dp_maple_ref_grid_128_pls_spcs.txt
2          I_nu(z)       double      negative       small_dp_I_neg_nu_dp.txt
3          I_nu(z)       quad        positive       small_qp_maple_ref_grid_128_pls_spcs.txt
4          I_nu(z)       quad        negative       small_qp_I_neg_nu.txt

5          J_nu(z)       double      positive       small_dp_maple_ref_grid_128_pls_spcsj.txt
6          J_nu(z)       double      negative       small_dp_maple_ref_grid_64j_neg.txt
7          J_nu(z)       quad        positive       small_qp_maple_ref_grid_128_pls_spcsj.txt
8          J_nu(z)       quad        negative       small_qp_maple_ref_grid_128j_neg.txt

9          K_nu(z)       double      positive       small_dp_maple_ref_grid_128_pls_spcsk.txt
10         K_nu(z)       double      negative       small_dp_k_neg_nu_dp.txt
11         K_nu(z)       quad        positive       small_qp_maple_ref_grid_128_pls_spcsk.txt
12         K_nu(z)       quad        negative       small_qp_maple_ref_grid_128k_neg.txt

13         Y_nu(z)       double      positive       small_dp_Maple_ref_grid_128_y.txt
14         Y_nu(z)       double      negative       small_dp_maple_ref_grid_64y_neg.txt
15         Y_nu(z)       quad        positive       small_qp_Maple_ref_grid_128_y.txt
16         Y_nu(z)       quad        negative       small_qp_maple_ref_grid_128y_neg.txt

The driver writes accuracy and timing output files such as:
dp_acc_cmp_I.txt
qp_acc_cmp_I.txt
dp_acc_cmp_J.txt
qp_acc_cmp_J.txt
dp_acc_cmp_K.txt
qp_acc_cmp_K.txt
dp_acc_cmp_Y.txt
qp_acc_cmp_Y.txt
timing_summary_I.txt
timing_summary_J.txt
timing_summary_K.txt
timing_summary_Y.txt
per_point_times_I.txt
per_point_times_J.txt
per_point_times_K.txt
per_point_times_Y.txt


## Authors and contributors

Mofreh R. Zaghloul  
Department of Physics, United Arab Emirates University

Steven G. Johnson  
Department of Mathematics, Massachusetts Institute of Technology


## Citation

If you use this package in published work, please cite the corresponding paper, software repository, or documentation once.
- Efficient Multi-Precision Computation of Bessel Functions for Real Orders and Complex Arguments with Fortran Implementation -- Part I: The Modified Bessel Function of the First Kind, I_nu(z)   https://arxiv.org/abs/2505.09770
- Efficient Multi-Precision Computation of Bessel Functions for Real Orders and Complex Arguments with Fortran Implementation -- Part II: The Modified Bessel Function of the Second Kind, K_nu(z)  https://arxiv.org/abs/2606.14839
- Efficient Multi-Precision Computation of Bessel Functions for Real Orders and Complex Arguments with a Fortran Implementation -- Part III: Regular Bessel Functions of the First and Second Kinds  J_nu(z) and Y_nu(z) 

