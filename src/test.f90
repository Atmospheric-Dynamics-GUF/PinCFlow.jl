module types

  implicit none

  type ice_rayType2
    real :: wwp ! wPrime: vertical vertical GW fluctuation
    real :: epp ! expPrime: Exner pressure GW fluctuation
    real :: thp ! thetaPrime: potential temperature GW fluctuation
    real :: si

    !evolved ice fields
    real :: Ni, Qi, Qv
    !auxillary arrays required for RK integration
    real :: qNi, qQi, qQv
    !tendency ice fields
    real :: tNi, tQi, tQv

  end type ice_rayType2

end module types

module operators
  use types
  implicit none
  interface operator(+)
    module procedure scpsc
  end interface operator(+)

  contains

  function scpsc(cl1, cl2) result(sum)

    implicit none

    type(ice_rayType2), INTENT(IN) :: cl1, cl2
    type(ice_rayType2) :: sum

    sum%wwp = cl1%wwp + cl2%wwp
    sum%epp = cl1%epp + cl2%epp
    sum%thp = cl1%thp + cl2%thp
    sum%si = cl1%si + cl2%si

  end function scpsc

  subroutine test_one
    implicit none
    type(ice_rayType2) :: ice_ray, sum

    ice_ray%wwp = 1
    ice_ray%si = 2

    sum = ice_ray + ice_ray
    print *, sum%wwp, sum%si
  end subroutine test_one
end module operators

program test
  use types
  use operators
  implicit none

  call test_one
end program test