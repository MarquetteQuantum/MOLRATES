!--------------------------------------------------------------------
! This program computes rate coefficients for transitions in H2O + H2O system.
! For more details and to cite this work, please refer to:
! - Bikramaditya Mandal et al, 2024, to be published.
! Also see:
! - Bikramaditya Mandal and Dmitri Babikov, 2023, Astron Astrophys 671, A51.
! - Bikramaditya Mandal and Dmitri Babikov, 2023, Astron Astrophys 678, A51.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! This module "module_numerical_recipies" contains
! FORTRAN functions and subroutines from Numerical Recipes.
! Please do not change anything in this module!
!--------------------------------------------------------------------

module module_numerical_recipies
  use iso_fortran_env, only: real64
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Solves for a vector u of size N the tridiagonal linear set using a
! serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N,
! while a and c (off-diagonal elements) are size N − 1.
!-------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE tridag(a,b,c,r,u)
    REAL(real64), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(real64), DIMENSION(:), INTENT(OUT) :: u
    REAL(real64), DIMENSION(size(b)) :: gam ! One vector of workspace, gam is needed.
    INTEGER :: n,j
    REAL(real64) :: bet

    if (any(size(a) + 1 /= [size(b), size(c) + 1, size(r), size(u)])) then
      stop 'Size mismatch in tridag'
    end if
    n = size(a) + 1
    bet=b(1)
    if (bet == 0) then
      ! If this happens then you should rewrite your equations as a set of order N − 1, with u2 trivially eliminated.
      stop 'tridag_ser: Error at code stage 1'
    end if
    u(1)=r(1)/bet
    do j=2,n ! Decomposition and forward substitution.
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      if (bet == 0) then
        stop 'tridag_ser: Error at code stage 2'
      end if
      u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1 ! Backsubstitution.
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi), with x1 <
! x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
! function at points 1 and N, respectively, this routine returns an array y2 of length N
! that contains the second derivatives of the interpolating function at the tabulated points
! xi. If yp1 and/or ypn are equal to 10^30 or larger, the routine is signaled to set the
! corresponding boundary condition for a natural spline, with zero second derivative on that
! boundary.
!-------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE spline(x,y,yp1,ypn,y2)
    REAL(real64), DIMENSION(:), INTENT(IN) :: x,y
    REAL(real64), INTENT(IN) :: yp1,ypn
    REAL(real64), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER :: n
    REAL(real64), DIMENSION(size(x)) :: a,b,c,r

    if (size(x) /= size(y) .or. size(x) /= size(y2)) then
      stop 'Size of x, y and y2 has to be the same in spline'
    end if
    n = size(x)
    c(1:n-1)=x(2:n)-x(1:n-1) ! Set up the tridiagonal equations.
    r(1:n-1)=6*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2*(c(2:n-1)+a(2:n-1))
    b(1)=1
    b(n)=1
    if (yp1 > 0.99d30) then ! The lower boundary condition is set either to be "natural"
      r(1)=0
      c(1)=0
    else ! or else to have a specified first derivative.
      r(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      c(1)=0.5d0
    end if ! The upper boundary condition is set either to be "natural"
    if (ypn > 0.99d30) then
      r(n)=0
      a(n)=0
    else ! or else to have a specified first derivative.
      r(n)=(-3/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
      a(n)=0.5d0
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given an array xx(1:N), and given a value x, returns a value j such that x is between
! xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
! j = N is returned to indicate that x is out of range.
!-------------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION locate(xx,x)
    REAL(real64), DIMENSION(:), INTENT(IN) :: xx
    REAL(real64), INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd

    n=size(xx)
    ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise.
    jl=0 ! Initialize lower
    ju=n+1 ! and upper limits.
    do
      if (ju-jl <= 1) then
        exit
      end if
      jm=(ju+jl)/2 ! Compute a midpoint,
      if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm ! and replace either the lower limit
      else
        ju=jm ! or the upper limit, as appropriate.
      end if
    end do
    if (x == xx(1)) then ! set the output, being careful with the endpoints.
      locate=1
    else if (x == xx(n)) then
      locate=n-1
    else
      locate=jl
    end if
  END FUNCTION locate

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given the arrays xa and ya, which tabulate a function (with the xa_i’s in increasing or
! decreasing order), and given the array y2a, which is the output from spline, and
! given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
! and y2a are all of the same size.
!-------------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION splint(xa,ya,y2a,x)
    REAL(real64), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(real64), INTENT(IN) :: x
    REAL(real64) :: splint
    INTEGER :: khi,klo,n
    REAL(real64) :: a,b,h

    if (size(xa) /= size(ya) .or. size(xa) /= size(y2a)) then
      write(*,'(a,3(2x,i0))') 'Sizes of xa, ya and y2a have to be equal in splint', size(xa), size(ya), size(y2a)
          stop
    end if
    n=size(xa)
    ! We will find the right place in the table by means of locate’s bisection algorithm. This is
    ! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
    ! order, and closely spaced, one would do better to store previous values of klo and khi and
    ! test if they remain appropriate on the next call.
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1 ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0) then
      stop "Bad xa input in splint. The xa's must be distinct"
    end if
    a=(xa(khi)-x)/h ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
  END FUNCTION splint
        
end module

module module_engines
  implicit none
  integer, parameter :: n_e = 6, n_t = 882, n_qp = 22, n_qo = 21, n_q = 38                                ! n_e = #energy, n_t = #transitions, n_qp = #para_states, n_qo = #ortho_states, n_q = #quencher_states
  integer :: j1(n_t), ka1(n_t), kc1(n_t), j2(n_t), ka2(n_t), kc2(n_t)
  integer chp(n_q), jp(n_q), kap(n_q), kcp(n_q)
  integer cho(n_q), jo(n_q), kao(n_q), kco(n_q)
  real*8 :: T(n_e), sigma(n_e, n_t), ep(n_q), eo(n_q)
  real*8 :: sigma_para(n_e, n_qp, n_qp), sigma_ortho(n_e, n_qo, n_qo), sigma_new(n_e, n_t/2, 2)
  integer :: j_i(n_t/2), j_f(n_t/2), ka_i(n_t/2), ka_f(n_t/2), kc_i(n_t/2), kc_f(n_t/2)
  real*8 :: e_i(n_t/2), e_f(n_t/2)
  real*8 :: spln_der(n_e, n_t), new_x(n_e, n_t/2), new_y(n_e, n_t/2)
  integer, parameter :: n_e_dE = n_e-1
  real*8 :: spln_derr(n_e_dE, n_t), new_xx(n_e_dE, n_t/2), new_yy(n_e_dE, n_t/2)
  real*8 :: new_y1(n_e, n_t/2), new_yy1(n_e_dE, n_t/2), common_terms
  real*8, parameter :: Ktown = 0.6950345740d0, pi = dacos(-1.0d0), speed = 4.84997439836149d-13
  real*8 :: rp, ro
  
  contains
  
  subroutine generate_TACS(Temp_rot, ratio)
    implicit none
    integer funit, i, j, k, nf
    real*8 Temp_rot, kt, Q, wp(n_q), wo(n_q), rowp(n_q), rowo(n_q)
    real*8, dimension(6) :: U_val = [133.33d0, 200.d0, 266.67d0, 400.d0, 533.33d0, 707.96d0]
    real*8 tacs, sumd_sigma, ratio
    character (len = 100) :: fname
    logical :: initialize_TACS = .false.
  
    if(.not. initialize_TACS) then
      if(Temp_rot < 5.0d0 .or. Temp_rot > 1000.d0) then
        print*, "The required temperature is beyond capabilities of this program."
        print*, "Please indicated the temperature range of 5 <= T <= 1000 K."
        print*, "Program will stop now."
        stop
      end if

      rp = 1.0d0 / (ratio + 1.0d0)
      ro = ratio / (ratio + 1.0d0)
      print*, "para & ortho portions: ", rp, ro

!--------------------------------------------------------------------
! to read rotational quantum numbers of states and their energies
!--------------------------------------------------------------------
      open(newunit = funit, file = 'para_Water_Energies.dat', status = 'old')
      read(funit, *)                                                                                                ! Reading header
      do i = 1, n_q
        read(funit, '(i5, i11, i11, i11, f11.3)') chp(i), jp(i), kap(i), kcp(i), ep(i)
      end do
      close(funit)
      
      open(newunit = funit, file = 'ortho_Water_Energies.dat', status = 'old')
      read(funit, *)                                                                                                ! Reading header
      do i = 1, n_q
        read(funit, '(i5, i11, i11, i11, f11.3)') cho(i), jo(i), kao(i), kco(i), eo(i)
      end do
      close(funit)
  
      initialize_TACS = .true.
    end if

    kt = Temp_rot*Ktown

    Q = 0.d0
    do i = 1, n_q
      wp(i) = ( 2d0 * jp(i) + 1d0 ) * exp(-ep(i) / kt)
      Q = Q + wp(i)
    end do
    write(*, '(2(a, f0.2))') "Partition Function Computed for p-H2O for rotational Temp. ", Temp_rot, " K= ", Q
    
    do i = 1, n_q
      rowp(i) = wp(i) / Q
    end do
    
    Q = 0.d0
    do i = 1, n_q
      wo(i) = ( 2d0 * jo(i) + 1d0 ) * exp(-eo(i) / kt)
      Q = Q + wo(i)
    end do
    write(*, '(2(a, f0.2))') "Partition Function Computed for o-H2O for rotational Temp. ", Temp_rot, " K= ", Q
    
    do i = 1, n_q
      rowo(i) = wo(i) / Q
    end do

!--------------------------------------------------------------------
! making TACS here by product of summed quencher final state and averaged over initial quencher state
!--------------------------------------------------------------------
    open(1002, file = "MQCT_TACS.dat", action = 'write')
    do nf = 1, 6
      write(1002,'(a27,f16.8)') "Collision Energy (cm^-1) = ", U_val(nf)
      write(1002,'(2(3(a4,1x),10x),a19)') "J1", "KA1", "KC1", "J2", "KA2", "KC2", "TACS (A^2)"
      write(fname, '(a, i0)') "MQCT_data_", nf
      open(112, file = trim(fname), status = 'old', action = 'read')
      do i = 1, n_qp
        do k = 1, n_qp
          tacs = 0.d0
          read(112, *)
          read(112, *)
          do j = 1, n_q
            read(112, *) sumd_sigma
            tacs = tacs + sumd_sigma * rowp(j) * rp
          end do
          do j = 1, n_q
            read(112,*) sumd_sigma
            tacs = tacs + sumd_sigma * rowo(j) * ro
          end do
          if(i.ne.k) write(1002, '(2(3(i4,1x),10x),e19.12)') jp(i), kap(i), kcp(i), jp(k), kap(k), kcp(k), tacs
        end do
      end do
      do i = 1, n_qo
        do k = 1, n_qo
          tacs = 0.d0
          read(112, *)
          read(112, *)
          do j = 1, n_q
            read(112, *) sumd_sigma
            tacs = tacs + sumd_sigma * rowp(j) * rp
          end do
          do j = 1, n_q
            read(112,*) sumd_sigma
            tacs = tacs + sumd_sigma * rowo(j) * ro
          end do
          if(i.ne.k) write(1002, '(2(3(i4,1x),10x),e19.12)') jo(i), kao(i), kco(i), jo(k), kao(k), kco(k), tacs
        end do
      end do
      close(112)
    end do
    close(1002)

  end subroutine
  
  subroutine integrate(Temp_req, input_data_type, rji, rkai, rkci, rjf, rkaf, rkcf, rrate1, rrate2)
    use module_numerical_recipies
    implicit none
    integer funit, i, j, k, kount, method, int_cycle
    integer :: rji(n_t/2), rkai(n_t/2), rkci(n_t/2), rjf(n_t/2), rkaf(n_t/2), rkcf(n_t/2)
    real*8 :: rrate1(n_t/2), rrate2(n_t/2)
    real*8 Temp_rot, Temp_req, sigma_out(n_t/2), E_step, E_sum, U, dE, E_coll, gamma, g1, g2
    real*8 sig_tmp_LE, sig_tmp_HE, sig_tmp, int_max_LE, int_max_HE, int_max, erf1, erf2
    real*8 y1, yn, x1, x2, yy1, yy2, x3, yy3, c_t_factor, mid1, mid2
    real*8 :: left(n_e, n_t/2), right(n_e, n_t/2), e_prime, e_not_prime, int_min(n_t/2)
    real*8 :: a_lowE(n_t/2), b_lowE(n_t/2), a_highE(n_t/2), b_highE(n_t/2), c_highE(n_t/2), slope(n_t/2)
    real*8 :: f_u2
    character dummy1
    character (len = 100) :: fname
    character (len = 1) :: input_data_type
    logical :: initialization = .false., prn_int = .false., err = .false.
    logical :: chk_dE(n_t/2) = .false., p2_p3(n_t/2) = .false., pnt_rm(n_t/2) = .false.

!---------------------------------------------------------------------------------------------------------
! This piece is to compute and print numerical interpolation and extrapolation
! It is enabled/disabled by a logical keyword prn_int
! by default prn_int = .false.
! To print the data, please make sure prn_int = .true.
    prn_int = .false.
    initialization = .false.
!---------------------------------------------------------------------------------------------------------


    c_t_factor = speed * dsqrt(Temp_req) / (Ktown * Temp_req)**2d0

!--------------------------------------------------------------------
! to read input data from files to interpolate further
!--------------------------------------------------------------------
    if(.not. initialization) then
      if(Temp_req < 5.0d0 .or. Temp_req > 1000.d0) then
        print*, "The required temperature is beyond capabilities of this program."
        print*, "Please indicated the temperature range of 5 <= T <= 1000 K."
        print*, "Program will stop now."
        stop
      end if

      open(newunit = funit, file = "MQCT_TACS.dat", status = 'old')
        do j = 1, n_e
          read(funit, '(a27,f16.8)') dummy1, T(j)
          read(funit, *)
          do i = 1, n_t
            read(funit, '(2(3(i4,1x),10x),e19.12)') j1(i), ka1(i), kc1(i), j2(i), ka2(i), kc2(i), sigma(j, i)
          end do
        end do
      close(funit)
          
!--------------------------------------------------------------------
! to reshape input data to be used further
!--------------------------------------------------------------------
      do k = 1, n_e
        kount = 0
        do i = 1, n_qp
          do j = 1, n_qp
            if(i /= j) then
              kount = kount + 1
              sigma_para(k, i, j) = sigma(k, kount)
            end if
          end do
        end do
        do i = 1, n_qo
          do j = 1, n_qo
            if(i /= j) then
              kount = kount + 1
              sigma_ortho(k, i, j) = sigma(k, kount)
            end if
          end do
        end do
        if(kount /= n_t) then
          print*, "ERROR in reshaping data (a).", kount, n_t, k
          stop
        end if
      end do
          
      do k = 1, n_e
        kount = 0
        do i = 2, n_qp
          do j = 1, i-1
            kount = kount + 1
            sigma_new(k, kount, 1) = sigma_para(k, i, j)
            sigma_new(k, kount, 2) = sigma_para(k, j, i)
            if (k == 1) then
              j_i(kount) = jp(i)
              j_f(kount) = jp(j)
              e_i(kount) = ep(i)
              e_f(kount) = ep(j)
              ka_i(kount) = kap(i)
              ka_f(kount) = kap(j)
              kc_i(kount) = kcp(i)
              kc_f(kount) = kcp(j)
            end if
          end do
        end do
        do i = 2, n_qo
          do j = 1, i-1
            kount = kount + 1
            sigma_new(k, kount, 1) = sigma_ortho(k, i, j)
            sigma_new(k, kount, 2) = sigma_ortho(k, j, i)
            if (k == 1) then
              j_i(kount) = jo(i)
              j_f(kount) = jo(j)
              e_i(kount) = eo(i)
              e_f(kount) = eo(j)
              ka_i(kount) = kao(i)
              ka_f(kount) = kao(j)
              kc_i(kount) = kco(i)
              kc_f(kount) = kco(j)
            end if
          end do
        end do
        if(kount /= n_t/2) then
          print*, "ERROR in reshaping data (b).", kount, n_t/2, k
          stop
        end if
      end do

      chk_dE = .false.
      open(1, file = "Check_Low_Energy.dat", action = 'read')
      do i = 1, n_t/2
        read(1, *) j
        if(j == 1) chk_dE(i) = .true.
      end do
      close(1)

      pnt_rm = .false.
      slope = 0.d0
      open(1, file = "Check_High_Energy.dat", action = 'read')
      do i = 1, n_t/2
        read(1, *) slope(i)
        if(slope(i) /= 0.d0) pnt_rm(i) = .true.
      end do
      close(1)

      if(input_data_type == 'q') then
        print*, 'Using QUENCHING as input'
      else if(input_data_type == 'e') then
        print*, 'Using EXCITATION as input'
      else if(input_data_type == 'a') then
        print*, 'Using AVERAGE of quenching and excitation as input'
      else
        print*, "No Input data-type or wrong type selected."
        print*, "Input data type are indicated as:"
        print*, "'q' for quenching,"
        print*, "'e' for excitation, and"
        print*, "'a' for average of quenching and excitation."
        print*, "Program will stop now."
        stop
      end if

      initialization = .true.
    end if

!--------------------------------------------------------------------
! this piece below is to construct different methods
! as indicated by case number
!--------------------------------------------------------------------
    method = 3
    select case (method)
    case (1)
      ! print*, "Following Method 1:"
            
    case (2)
      ! print*, "Following Method 2:"
      do k = 1, n_e
        ! U = 4d0*Ktown*T(k)/pi
        U = T(k)
        do i = 1, n_t/2
          dE = abs(e_f(i) - e_i(i))
          if(k == 1) int_min(i) = dE
          e_prime = U + dE/2d0 + dE**2/16d0/U
          left(k, i) = (2d0*j_i(i) + 1d0)* (e_prime - dE)* sigma_new(k, i, 1)* exp(-e_prime/ (Ktown* Temp_req))
          right(k, i) = (2d0*j_f(i) + 1d0)* e_prime* sigma_new(k, i, 2)* exp(-e_prime/ (Ktown* Temp_req))
          new_x(k, i) = e_prime
          new_y(k, i) = (left(k, i) + right(k, i))* 0.50d0
        end do
      end do
                
      ! open(1, file = "M2_data.dat")
      do i = 1, n_t/2
        do k = 1, n_e
          ! write(1,*) new_x(k, i), new_y(k, i)
        end do
      end do
      ! close(1)
                
    case (3)
      ! print*, "Following Method 3:"
      do k = 1, n_e
        ! U = 4d0*Ktown*T(k)/pi
        U = T(k)
        do i = 1, n_t/2
          dE = abs(e_f(i) - e_i(i))
          if(k == 1) int_min(i) = dE/4d0
          e_prime = U + dE/2d0 + dE**2/16d0/U
          e_not_prime = U - dE/2d0 + dE**2/16d0/U
          gamma = (dE/(4d0*U))**2
          g1 = 1.d0 + gamma
          g2 = 1.d0 - gamma
          common_terms = U* exp(-U/(Ktown* Temp_req)* g1)* g2
          left(k, i) =  (2d0*j_i(i) + 1d0)* sigma_new(k, i, 1)* common_terms
          right(k, i) = (2d0*j_f(i) + 1d0)* sigma_new(k, i, 2)* common_terms
          new_x(k, i) = U
          if(input_data_type == 'q') then
            new_y(k, i) = left(k, i)
          else if(input_data_type == 'e') then
            new_y(k, i) = right(k, i)
          else if(input_data_type == 'a') then
            new_y(k, i) = (left(k, i) + right(k, i))* 0.50d0
          end if
          new_y1(k, i) = new_y(k, i) / common_terms
        end do
      end do

      if(prn_int) then
        open(2, file = "M3_data.dat")
        do i = 1, n_t/2
          do k = 1, n_e
          write(2,*) new_x(k, i), new_y(k, i)
          end do
        end do
        close(2)
      end if
                
    case (4)
      ! print*, "Following Method 4:"
                
    case default
      print*, "No method selected.", method
      stop
    end select

!--------------------------------------------------------------------
! This is to filter out the first point if dE/4 > U
!--------------------------------------------------------------------
    do i = 1, n_t/2
      if(chk_dE(i)) then
        do k = 1, n_e_dE
          new_xx(k, i)  = new_x(k+1, i)
          new_yy(k, i)  = new_y(k+1, i)
          new_yy1(k, i) = new_y1(k+1, i)
          new_x (k, i)  = new_x(k+1, i)
          new_y (k, i)  = new_y(k+1, i)
          new_y1 (k, i) = new_y1(k+1, i)
        end do
      end if
    end do

!--------------------------------------------------------------------
! this is to compute the derivative to be used in the
! splint program later for the new constructed data
!--------------------------------------------------------------------             
    do i = 1, n_t/2
      if(.not. chk_dE(i)) then
        y1 = (new_y1 (2, i) - new_y1 (1, i)) / (new_x (2, i) - new_x (1, i))
        yn = (new_y1 (n_e - 1, i) - new_y1 (n_e, i)) / (new_x (n_e - 1, i) - new_x (n_e, i))
        call spline(new_x(:, i), new_y1 (:, i), y1, yn, spln_der(:,i))
      else if(chk_dE(i)) then
        y1 = (new_yy1 (2, i) - new_yy1 (1, i)) / (new_xx (2, i) - new_xx (1, i))
        yn = (new_yy1 (n_e_dE - 1, i) - new_yy1 (n_e_dE, i)) / (new_xx (n_e_dE - 1, i) - new_xx (n_e_dE, i))
        call spline(new_xx(:, i), new_yy1 (:, i), y1, yn, spln_derr(:,i))
      else
        print*, "Error in Spline calling.", i
        stop
      end if
    end do

!--------------------------------------------------------------------
! this is to compute the extrapolation at the low energy regime
!--------------------------------------------------------------------
    do i = 1, n_t/2
      if(.not. chk_dE(i)) then
        x1  = new_x(1, i)
        x2  = new_x(2, i)
        yy1 = new_y(1, i)
        yy2 = new_y(2, i)
      else if(chk_dE(i)) then
        x1  = new_xx(1, i)
        x2  = new_xx(2, i)
        yy1 = new_yy(1, i)
        yy2 = new_yy(2, i)
      else
        print*, "Error in low energy extrapolation.", i
        stop
      end if
      b_lowE(i) = -(x1-x2)/ log(yy1/ yy2* (x2-int_min(i))/ (x1-int_min(i)))
      a_lowE(i) = (x1-int_min(i))/ yy1* exp(-(x1-int_min(i))/b_lowE(i))
    end do

!--------------------------------------------------------------------
! this is to check % of adjusting last point at the high energy regime
!--------------------------------------------------------------------
    ! p2_p3 = .false.
    ! do i = 1, n_t/2
    ! kount = 0
      ! if(.not. chk_dE(i)) then
        ! x1  = new_x(n_e-2, i)
        ! x2  = new_x(n_e-1, i)
        ! x3  = new_x(n_e, i)
        ! yy1 = new_y(n_e-2, i)
        ! yy2 = new_y(n_e-1, i)
        ! yy3 = new_y(n_e, i)
        ! b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                     ! (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                     ! ((x3**2 - x2**2) * log(yy3 / yy1)))
        ! c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                     ! (((x3 - x2) * log(yy3 / yy1)) - &
                     ! ((x3 - x1) * log(yy3 / yy2)))
        ! a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
        ! if(c_highE(i) < 0.d0 .or. c_highE(i) > 1.d6) pnt_rm(i) = .true.
        ! if(pnt_rm(i)) then
          ! do while (.not. p2_p3(i))
            ! kount = kount + 1
            ! x1  = new_x(n_e-2, i)
            ! x2  = new_x(n_e-1, i)
            ! x3  = new_x(n_e, i)
            ! yy1 = new_y(n_e-2, i)
            ! yy2 = new_y(n_e-1, i)
            ! yy3 = new_y(n_e, i) * (1.d0 - kount / 10d0)
            ! b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                         ! (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                         ! ((x3**2 - x2**2) * log(yy3 / yy1)))
            ! c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                         ! (((x3 - x2) * log(yy3 / yy1)) - &
                         ! ((x3 - x1) * log(yy3 / yy2)))
            ! a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
            ! if(c_highE(i) > 0.d0 .and. c_highE(i) < 1.d6) then
              ! p2_p3(i) = .true.
              ! if(slope(i) == 0.d0 .or. slope(i) < 0.00009d0) then
                ! print*, i, slope(i), 1.d0 - kount / 10d0
                ! slope(i) = 1.d0 - kount / 10d0
              ! end if
              ! if(slope(i) > 1.d0 - kount / 10d0) then
                ! print*, i, slope(i), 1.d0 - kount / 10d0
                ! slope(i) = 1.d0 - kount / 10d0
              ! end if
            ! end if
          ! end do
        ! end if
      ! else if(chk_dE(i)) then
        ! x1  = new_xx(n_e_dE-2, i)
        ! x2  = new_xx(n_e_dE-1, i)
        ! x3  = new_xx(n_e_dE, i)
        ! yy1 = new_yy(n_e_dE-2, i)
        ! yy2 = new_yy(n_e_dE-1, i)
        ! yy3 = new_yy(n_e_dE, i)
        ! b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                     ! (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                     ! ((x3**2 - x2**2) * log(yy3 / yy1)))
        ! c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                     ! (((x3 - x2) * log(yy3 / yy1)) - &
                     ! ((x3 - x1) * log(yy3 / yy2)))
        ! a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
        ! if(c_highE(i) < 0.d0 .or. c_highE(i) > 1.d6) pnt_rm(i) = .true.
        ! if(pnt_rm(i)) then
          ! do while (.not. p2_p3(i))
            ! kount = kount + 1
            ! x1  = new_xx(n_e_dE-2, i)
            ! x2  = new_xx(n_e_dE-1, i)
            ! x3  = new_xx(n_e_dE, i)
            ! yy1 = new_yy(n_e_dE-2, i)
            ! yy2 = new_yy(n_e_dE-1, i)
            ! yy3 = new_yy(n_e_dE, i) * (1.d0 - kount / 10d0)
            ! b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                         ! (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                         ! ((x3**2 - x2**2) * log(yy3 / yy1)))
            ! c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                         ! (((x3 - x2) * log(yy3 / yy1)) - &
                         ! ((x3 - x1) * log(yy3 / yy2)))
            ! a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
            ! if(c_highE(i) > 0.d0 .and. c_highE(i) < 1.d6) then
              ! p2_p3(i) = .true.
              ! if(slope(i) == 0.d0 .or. slope(i) < 0.00009d0) then
                ! print*, i, slope(i), 1.d0 - kount / 10d0
                ! slope(i) = 1.d0 - kount / 10d0
              ! end if
              ! if(slope(i) > 1.d0 - kount / 10d0) then
                ! print*, i, slope(i), 1.d0 - kount / 10d0
                ! slope(i) = 1.d0 - kount / 10d0
              ! end if
            ! end if
              ! end do
            ! end if
      ! else
        ! print*, "Error in high energy point removal.", i
        ! stop
      ! end if
    ! ! print*, i, chk_dE(i), pnt_rm(i), p2_p3(i), kount, slope(i), &
    ! ! a_highE(i), b_highE(i), c_highE(i)
    ! end do
    ! p2_p3 = .false.
    ! open(1, file = "Check_High_Energy.dat", action = 'write')
    ! do i = 1, n_t/2
      ! if(pnt_rm(i)) then
        ! if(slope(i) == 0.d0 .or. slope(i) < 0.00090d0) stop "Error computing shift to last point."
        ! write(1, *) slope(i)
      ! else
        ! write(1, *) "0.d0"
      ! end if
    ! end do
    ! close(1)

!--------------------------------------------------------------------
! this is to compute the extrapolation at the high energy regime
!--------------------------------------------------------------------
    err = .false.
    do i = 1, n_t/2
      if(.not. chk_dE(i)) then
        if(.not. pnt_rm(i)) then
          x1  = new_x(n_e-2, i)
          x2  = new_x(n_e-1, i)
          x3  = new_x(n_e, i)
          yy1 = new_y(n_e-2, i)
          yy2 = new_y(n_e-1, i)
          yy3 = new_y(n_e, i)
          b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                       ((x3**2 - x2**2) * log(yy3 / yy1)))
          c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3 - x2) * log(yy3 / yy1)) - &
                       ((x3 - x1) * log(yy3 / yy2)))
          a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
        else if(pnt_rm(i)) then
          x1  = new_x(n_e-2, i)
          x2  = new_x(n_e-1, i)
          x3  = new_x(n_e, i)
          yy1 = new_y(n_e-2, i)
          yy2 = new_y(n_e-1, i)
          yy3 = new_y(n_e, i) * slope(i)
          b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                       ((x3**2 - x2**2) * log(yy3 / yy1)))
          c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3 - x2) * log(yy3 / yy1)) - &
                       ((x3 - x1) * log(yy3 / yy2)))
          a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
        end if
        ! if(c_highE(i) < 0.d0) then
          ! p2_p3(i) = .true.
          ! x1  = new_x(n_e-2, i)
          ! x2  = new_x(n_e-3, i)
          ! yy1 = new_y(n_e-2, i)
          ! yy2 = new_y(n_e-3, i)
          ! a_highE(i) = (yy2**(x1/ (x1-x2))) / (yy1**(x2/ (x1-x2)))
          ! b_highE(i) = -(x1-x2)/ log(yy1/ yy2)
          ! c_highE(i) = 0.d0
        ! end if
      else if(chk_dE(i)) then
        if(.not. pnt_rm(i)) then
          x1  = new_xx(n_e_dE-2, i)
          x2  = new_xx(n_e_dE-1, i)
          x3  = new_xx(n_e_dE, i)
          yy1 = new_yy(n_e_dE-2, i)
          yy2 = new_yy(n_e_dE-1, i)
          yy3 = new_yy(n_e_dE, i)
          b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                       ((x3**2 - x2**2) * log(yy3 / yy1)))
          c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3 - x2) * log(yy3 / yy1)) - &
                       ((x3 - x1) * log(yy3 / yy2)))
          a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
        else if(pnt_rm(i)) then
          x1  = new_xx(n_e_dE-2, i)
          x2  = new_xx(n_e_dE-1, i)
          x3  = new_xx(n_e_dE, i)
          yy1 = new_yy(n_e_dE-2, i)
          yy2 = new_yy(n_e_dE-1, i)
          yy3 = new_yy(n_e_dE, i) * slope(i)
          b_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3**2 - x1**2) * log(yy3 / yy2)) - &
                       ((x3**2 - x2**2) * log(yy3 / yy1)))
          c_highE(i) = ((x2 - x1) * (x3 - x2) * (x3 - x1)) / &
                       (((x3 - x2) * log(yy3 / yy1)) - &
                       ((x3 - x1) * log(yy3 / yy2)))
          a_highE(i) = yy3 * exp(x3 / b_highE(i)) * exp(x3**2 / c_highE(i))
        end if
        ! if(c_highE(i) < 0.d0) then
          ! p2_p3(i) = .true.
          ! x1  = new_xx(n_e_dE-1, i)
          ! x2  = new_xx(n_e_dE-2, i)
          ! yy1 = new_yy(n_e_dE-1, i)
          ! yy2 = new_yy(n_e_dE-2, i)
          ! a_highE(i) = (yy2**(x1/ (x1-x2))) / (yy1**(x2/ (x1-x2)))
          ! b_highE(i) = -(x1-x2)/ log(yy1/ yy2)
          ! c_highE(i) = 0.d0
        ! end if

      else
        print*, "Error in high energy extrapolation.", i
        stop
      end if

      if(c_highE(i) < 0.d0) then
        print*, "Cannot proceed as C^2 < 0. C^2 = ", c_highE(i), " , i = ", i
        err = .true.
        !stop
      end if
      if(c_highE(i) > 1.d6) then
        print*, "Cannot proceed as C^2 > 1.E6 C^2 = ", c_highE(i), " , i = ", i
        err = .true.
        !stop
      end if
      ! if(i == 127) print*, i, pnt_rm(i), chk_dE(i), a_highE(i), b_highE(i), c_highE(i)
      ! if(1/c_highE(i) < ctol) print*, i, a_highE(i), b_highE(i), c_highE(i)
    end do
    if(err) stop
      

!---------------------------------------------------------------------------------------------------------
! This piece is to compute and print numerical interpolation and extrapolation
! It is enabled/disabled by a logical keyword prn_int
! by default prn_int = .false.
! To print the data, please make sure prn_int = .true.
!---------------------------------------------------------------------------------------------------------
    if(prn_int) then
      open(3, file = "interpolation_exterpolation.out")
      int_max = 100000.0d0
      int_cycle = 1000
      do i = 1, n_t/2

!---------------------------------------------------------------------------------------------------------
! deciding values of different range of integration
!---------------------------------------------------------------------------------------------------------
        if(.not. chk_dE(i)) then
          mid1 = new_x(2, i)
          if(.not. pnt_rm(i)) then
            mid2 = new_x(n_e-1, i)
          else
            mid2 = new_x(n_e-2, i)
          end if
        else if(chk_dE(i)) then
          mid1 = new_xx(2, i)
          if(.not. pnt_rm(i)) then
            mid2 = new_xx(n_e_dE-1, i)
          else
            mid2 = new_xx(n_e_dE-2, i)
          end if
        end if

        E_step = (int_max - int_min(i))/ int_cycle
        E_sum = 0.0d0
        do j = 1, int_cycle
          E_coll = int_min(i) + (j-1)*E_step
          if(E_coll < mid1) sig_tmp = (E_coll - int_min(i))/a_lowE(i)*exp(- ((E_coll - int_min(i))/ b_lowE(i)))
          if(E_coll == int_min(i)) sig_tmp = 1.d-100
          if(E_coll >= mid1 .and. E_coll < mid2) then
            U = E_coll
            gamma = (int_min(i) / U)**2
            g1 = 1.d0 + gamma
            g2 = 1.d0 - gamma
            common_terms = U* exp(-U/(Ktown* Temp_req)* g1)* g2
            if(.not. chk_dE(i)) then
              sig_tmp = splint(new_x(:, i), new_y1(:, i), spln_der(:,i), E_coll) * common_terms
            else if(chk_dE(i)) then
              sig_tmp = splint(new_xx(:, i), new_yy1(:, i), spln_derr(:,i), E_coll) * common_terms
            end if
          else if(E_coll >= mid2) then
            !  if(p2_p3(i)) then
              !  sig_tmp = a_highE(i)* exp(-(E_coll/ b_highE(i)))
            !  else
              sig_tmp = a_highE(i)* exp( - (E_coll / b_highE(i)) - (E_coll**2 / c_highE(i)))
            !  end if
          end if

          if(sig_tmp < 0.d0) write(*,*) j_i(i), j_f(i), E_coll, sig_tmp
          if(sig_tmp /= sig_tmp) write(*,*) j_i(i), j_f(i), E_coll, sig_tmp
          E_sum = E_sum + sig_tmp*E_step
          write(3,*) E_coll, sig_tmp
        end do

      end do
      close(3)
    end if

!--------------------------------------------------------------------
! this is the main piece to compute integral following a traditional equations.
! here the range of intergration is also decided
!--------------------------------------------------------------------

    do i = 1, n_t/2

!--------------------------------------------------------------------
! This is computing the integral to the left end.
!--------------------------------------------------------------------
      int_max_LE = new_x(2, i)
      ! sig_tmp_LE = b_lowE(i) / a_lowE(i) * (b_lowE(i) - (b_lowE(i) - int_min(i) + int_max_LE) * &
      !              exp((int_min(i) - int_max_LE) / b_lowE(i)))
      f_u2 = (int_max_LE - int_min(i)) / a_lowE(i) * exp(- ((int_max_LE - int_min(i)) / b_lowE(i)))
      sig_tmp_LE = b_lowE(i)**2 / a_lowE(i) * (1.d0 - f_u2 * &
                   a_lowE(i) / (b_lowE(i) * (int_max_LE - int_min(i))) * &
                  (b_lowE(i) + int_max_LE - int_min(i))) 

!--------------------------------------------------------------------
! This is to check the analytical integral with numerical integration.
!--------------------------------------------------------------------
!          int_cycle = 20000
!          E_step = (int_max_LE - int_min(i))/ (int_cycle - 1)
!          E_sum = 0.0d0
!          do j = 1, int_cycle
!            E_coll = int_min(i) + (j-1)*E_step
!            E_sum = E_sum + (E_coll - int_min(i))/a_lowE(i)*exp(- ((E_coll - int_min(i))/ b_lowE(i))) * E_step
!          end do
!          print*, i, E_sum, sig_tmp_LE, abs(E_sum - sig_tmp_LE)/E_sum * 100
!          cycle

!---------------------------------------------------------------------------------------------------------
! deciding values of different range of integration
!---------------------------------------------------------------------------------------------------------
      if(.not. chk_dE(i)) then
        if(.not. pnt_rm(i)) then
          int_max_HE = new_x(n_e-1, i)
        else
          int_max_HE = new_x(n_e-2, i)
        end if
      else if(chk_dE(i)) then
        if(.not. pnt_rm(i)) then
          int_max_HE = new_xx(n_e_dE-1, i)
        else
          int_max_HE = new_xx(n_e_dE-2, i)
        end if
      end if

!--------------------------------------------------------------------
! This is computing the numerical integration for the spline piece.
!--------------------------------------------------------------------
      int_cycle = 200000
      E_step = (int_max_HE - int_max_LE)/ (int_cycle - 1)
      E_sum = 0.0d0
      do j = 1, int_cycle
        E_coll = int_max_LE + (j-1)*E_step
        U = E_coll
        gamma = (int_min(i) / U)**2
        g1 = 1.d0 + gamma
        g2 = 1.d0 - gamma
        common_terms = U* exp(-U/(Ktown* Temp_req)* g1)* g2
        if(.not. chk_dE(i)) then
          sig_tmp = splint(new_x(:, i), new_y1(:, i), spln_der(:,i), E_coll) * common_terms
        else if(chk_dE(i)) then
          sig_tmp = splint(new_xx(:, i), new_yy1(:, i), spln_derr(:,i), E_coll) * common_terms
        end if
        E_sum = E_sum + sig_tmp*E_step
      end do

!--------------------------------------------------------------------
! This is computing the integral to the right end.
!--------------------------------------------------------------------
      int_max = 50000.d0
      !  if(p2_p3(i)) then
        !  erf1 = exp(-int_max / b_highE(i))
        !  erf2 = exp(-int_max_HE / b_highE(i))
        !  sig_tmp_HE = a_highE(i) * b_highE(i) * (erf2 - erf1)
      !  else
        ! erf1 = erf((2d0 * b_highE(i) * int_max + c_highE(i)) / (2d0 * b_highE(i) * dsqrt(c_highE(i))))
        ! erf2 = erf((2d0 * b_highE(i) * int_max_HE + c_highE(i)) / (2d0 * b_highE(i) * dsqrt(c_highE(i))))
        f_u2 = a_highE(i)* exp( - (int_max_HE / b_highE(i)) - (int_max_HE**2 / c_highE(i)))
        erf1 = (2.d0 * b_highE(i) * int_max_HE + c_highE(i)) / (2.d0 * b_highE(i) * dsqrt(c_highE(i)))
        erf2 = bk_erf(erf1)
        sig_tmp_HE = 0.5d0 * dsqrt(pi * c_highE(i)) * erf2 * f_u2
      !  end if

!--------------------------------------------------------------------
! This is to check the analytical integral with numerical integration.
!--------------------------------------------------------------------
!      E_step = (int_max - int_max_HE)/ (int_cycle - 1)
!      sig_tmp = 0.d0
!      do j = 1, int_cycle
!        E_coll = int_max_HE + (j-1)*E_step
!        if(p2_p3(i)) then
!          sig_tmp = sig_tmp + a_highE(i)* exp(-(E_coll/ b_highE(i))) * E_step
!        else
!          sig_tmp = sig_tmp + a_highE(i)* exp( - (E_coll / b_highE(i)) - (E_coll**2 / c_highE(i))) * E_step
!        end if
!      end do
!      print*, i, p2_p3(i), sig_tmp, sig_tmp_HE, abs(sig_tmp - sig_tmp_HE)/sig_tmp*100, a_highE(i), b_highE(i), c_highE(i)
!      cycle

      E_sum = E_sum + sig_tmp_LE + sig_tmp_HE
      sigma_out(i) = E_sum

      ! print*, E_sum, sig_tmp_HE, abs(E_sum - sig_tmp_HE)/E_sum*100.d0
      ! write(*,'(i4,3x,3(i0),3x,3(i0),3x,2(e19.12,3x))')i, j_i(i), ka_i(i), kc_i(i), j_f(i), ka_f(i), kc_f(i), &
                 ! E_sum* exp(+(abs(e_f(i) - e_i(i)))/ (2d0* Ktown* Temp_req))/ (2d0* j_i(i) + 1) * c_t_factor, &
                 ! E_sum* exp(-(abs(e_f(i) - e_i(i)))/ (2d0* Ktown* Temp_req))/ (2d0* j_f(i) + 1) * c_t_factor
    end do

!--------------------------------------------------------------------
! This is to print the computed rates.
!--------------------------------------------------------------------
!    write(fname, '(a, 2(f0.5, a))') "Rates_for_Trot_", Temp_rot, "K_Tkin_", Temp_req, "K.out"
!    open(1, file = trim(fname))
!    write(1, '(a, f0.5)') "Temperature(K) = ", Temp_req
!    write(1, *) "Unit of Rates is cm^3/s"
!    write(1, '(2(3(a4,1x),2x),2(5x,a19))') "J1", "KA1", "KC1", "J2", "KA2", "KC2", "Quenching", "Excitation"
    do i = 1, n_t/2
!      write(1, '(2(3(i4,1x),2x),2(5x,e19.12))') j_i(i), ka_i(i), kc_i(i), j_f(i), ka_f(i), kc_f(i), &
!      sigma_out(i)* exp(+(abs(e_f(i) - e_i(i)))/ (2d0* Ktown* Temp_req))/ (2d0* j_i(i) + 1) * c_t_factor, &
!      sigma_out(i)* exp(-(abs(e_f(i) - e_i(i)))/ (2d0* Ktown* Temp_req))/ (2d0* j_f(i) + 1) * c_t_factor
  
      rji(i) = j_i(i)
      rkai(i) = ka_i(i)
      rkci(i) = kc_i(i)
      rjf(i) = j_f(i)
      rkaf(i) = ka_f(i)
      rkcf(i) = kc_f(i)
      rrate1(i) = sigma_out(i)* exp(+(abs(e_f(i) - e_i(i)))/ (2d0* Ktown* Temp_req))/ (2d0* j_i(i) + 1) * c_t_factor
      rrate2(i) = sigma_out(i)* exp(-(abs(e_f(i) - e_i(i)))/ (2d0* Ktown* Temp_req))/ (2d0* j_f(i) + 1) * c_t_factor
    end do
!    close(1)

    return
  end subroutine

  function bk_erf(x)
    implicit none
    real*8 x, bk_erf, t
    real*8, parameter :: p = 0.3275911d0
    real*8, parameter :: a1 = 0.254829592d0, a2 = -0.284496736d0, a3 = 1.421413741d0, a4 = -1.453152027, a5 = 1.061405429d0

    t = 1.d0 / (1.d0 + p * x)
    bk_erf = a1 * t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5
  end function
  
end module

