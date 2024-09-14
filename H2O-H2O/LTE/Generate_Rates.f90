!--------------------------------------------------------------------
! This program computes rate coefficients for transitions in H2O + H2O system.
! It is designed for LTE case when kinetic temperature is equal to the rotational temperature, i.e., Temp_kin = Temp_rot.
! For more details and to cite this work, please refer to:
! - Bikramaditya Mandal et al, 2024, Astronomy & Astrophysics, 688, A208.
! - Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
! - Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! This is the main driver program.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! User should change the value of input variable Temp to the desired one.
! E.g., the current value is set to 300 K.
! The unit of input temterature is K and the unit of output rate is cm^3/s.
!--------------------------------------------------------------------

program Generate_Rates
  use module_engines
  implicit none
  integer i, j, k
  character (len = 1) :: input_data_type
  character (len = 100) :: label

! The total number of transitions (dont change):
  integer, parameter :: n = 441

! Rotational and kinetic temperatures (in unit of Kelvin)
! are set equal to each other to obtain LTE:
  real*8 :: Temp

! Initial and final quantum numbers of the target H2O:
  integer :: ji(n), kai(n), kci(n), jf(n), kaf(n), kcf(n)

! Quenching and excitation rate coefficients:
  real*8 :: rateq(n), ratee(n)

! Value of ortho-to-para ratio:
  real*8, parameter :: ortho_to_para_ratio = 3.0d0

! Which MQCT cross sections are employed to predict rate coefficients:
! 'e' for excitaion
! 'q' for quenching
! 'a' for average (recomended)
  input_data_type = 'a' 
 
! Generate thermally averaged cross sections:
  Temp = 300.d0
  call generate_TACS(Temp, ortho_to_para_ratio)

! Compute rate coefficients:
  call integrate(Temp, input_data_type, ji, kai, kci, jf, kaf, kcf, rateq, ratee)

! Simple output:
  write(*, *) "Temp_rot (K) = Temp_kin (K) = ", Temp
  write(*, '(2(3(a4,1x),2x),2(5x,a19))') "J1", "KA1", "KC1", "J2", "KA2", "KC2", "k_Quench (cm^3/s)", "k_Excite (cm^3/s)"
  do i = 1, n
    write(*, '(2(3(i4,1x),2x),2(5x,e19.12))') &
          ji(i), kai(i), kci(i), jf(i), kaf(i), kcf(i), rateq(i), ratee(i)
  end do

end program
