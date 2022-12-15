module tracer_module

  use type_module
  use atmosphere_module

  implicit none

  private

  public :: setup_tracer

contains

  subroutine setup_tracer (var)

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    real :: tracer0 ! initial max. tracer ampl.
    real :: xt0 , yt0 , zt0
    real :: sigtx, sigty, sigtz

    integer :: ii , jj , kk

    select case (tracerSetup)

    case ( "gaussian_tracer" )

       tracer0 = 1.0
       xt0 = nx/2.*dx
       yt0 = ny/2.*dy
       zt0 = nz/2.*dz
       sigtx = 0.05
       sigty = 0.06
       sigtz = 0.1

       do ii = 1,nx
          do jj = 1,ny
             do kk = 1,nz
                var(ii,jj,kk,iVart) = tracer0* &
                     & exp(-(x(ii)-xt0)**2./(2.*sigtx**2)) * &
                     & exp(-(y(jj)-yt0)**2./(2.*sigty**2)) * &
                     & exp(-(z(kk)-zt0)**2./(2.*sigtz**2))
             end do
          end do
       end do
       

    case ( "layer_tracer" )

       tracer0 = 1.0
       zt0 = nz/2.*dz
       
       do kk = 1,nz
          var(:,:,kk,iVart) = tracer0*exp(-(z(kk)-zt0)**2./(2.*2.**2.))
       end do
       

    case default
       stop "tracer.f90: unkown tracer distribution setup"
    end select

  end subroutine setup_tracer
  
end module tracer_module

