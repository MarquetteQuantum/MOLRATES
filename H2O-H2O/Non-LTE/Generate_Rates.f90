!--------------------------------------------------------------------
! This program computes rate coefficients for transitions in H2O + H2O system.
! For more details and to cite this work, please refer to:
! - Bikramaditya Mandal et al, Rotational state-to-state transition rate coefficients for H$_{2}$O + H$_{2}$O collisions at nonequilibrium conditions, Astronomy & Astrophysics, 2024.
! - Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
! - Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! This is the main driver program.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! User should change the values of input variables Temp_rot and Temp_kin to the 
! desired ones. E.g., the current values are set to 250 and 300 K, respectively.
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

! Rotational and kinetic temperatures (in unit of Kelvin):
  real*8 :: Temp_rot, Temp_kin 

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
  Temp_rot = 250.d0
  call generate_TACS(Temp_rot, ortho_to_para_ratio)

! Compute rate coefficients:
  Temp_kin = 350.d0
  call integrate(Temp_kin, input_data_type, ji, kai, kci, jf, kaf, kcf, rateq, ratee)

! Simple output:
  write(*, *) "Temp_rot (K), Temp_kin (K) = ", Temp_rot, Temp_kin
  write(*, '(2(3(a4,1x),2x),2(5x,a19))') "J1", "KA1", "KC1", "J2", "KA2", "KC2", "k_Quench (cm^3/s)", "k_Excite (cm^3/s)"
  do i = 1, n
    write(*, '(2(3(i4,1x),2x),2(5x,e19.12))') &
          ji(i), kai(i), kci(i), jf(i), kaf(i), kcf(i), rateq(i), ratee(i)
  end do

end program
