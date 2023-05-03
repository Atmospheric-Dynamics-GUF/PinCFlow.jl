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
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz)::rho

    real :: tracer0 ! initial max. tracer ampl.
    real :: xt0 , yt0 , zt0
    real :: sigtx, sigty, sigtz
    real :: k0t , m0t
    real :: alpha
    real :: zcart
    real :: dz_tr

    integer :: ii , jj , kk

    select case (tracerSetup)


    case ( "thin_layer" )

       tracer0 = 1.0

       var(:,:,:,iVart) = 0.
       var(:,:,nz/2,iVart) = tracer0
       

    case ( "gaussian_tracer_2D" )

       tracer0 = 1.0
       xt0 = 0.5*nx*dx!0.75*nx*dx
       zt0 = nz/2.*dz
       sigtx = 0.1*nx*dx
       sigtz = 0.05*nz*dz

       do ii = 1,nx
          do jj = 1,ny
             do kk = 1,nz
                if (topography) then
                   zcart = jac(ii,jj,kk) * z(kk) + topography_surface(ii,jj)
                else
                   zcart = z(kk)
                end if
                var(ii,jj,kk,iVart) = tracer0* &
                     & exp(-(x(ii)-xt0)**2./(2.*sigtx**2)) * &
                     & exp(-(zcart-zt0)**2./(2.*sigtz**2))
             end do
          end do
       end do

    case ( "gaussian_tracer_3D" )

       tracer0 = 1.0
       xt0 = 0.75*nx*dx
       yt0 = 0.75*ny*dy
       zt0 = nz/2.*dz
       sigtx = 0.1*nx*dx
       sigty = 0.1*ny*dy
       sigtz = 0.05*nz*dz

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

    case ( "sin_cos_tracer" )

       tracer0 = 1.0
       
       k0t = pi/(mountainWidth_dim / lRef)
       m0t = pi/(lz_dim(1) / lRef)
       
       do ii = 1,nx
          do jj = 1,ny
             do kk = 1,nz
                if (topography) then
                   zcart = jac(ii,jj,kk) * z(kk) + topography_surface(ii,jj)
                else
                   zcart = z(kk)
                end if
                var(ii,jj,kk,iVart) = tracer0/2.* &
                     &(1+ cos(k0t*x(ii))) * sin(m0t*zcart)**2.
             end do
          end do
       end do

    case( "increase_in_z_tracer" )

       
       do kk = 1,nz
          do jj = 1,ny
             do ii = 1,nx
                if (fluctuationMode) then
                   if (topography) then
                      rho(ii,jj,kk) = (var(ii,jj,kk,1) + rhoStratTFC(ii,jj,kk))
                   else
                      rho(ii,jj,kk) = var(ii,jj,kk,1) + rhoStrat(kk)
                   end if
                else
                   rho(ii,jj,kk) = var(ii,jj,kk,1)
                end if
             end do
          end do
       end do

       dz_tr = 0.1 * dz

       print*, "dz_tr = ", dz_tr
       
       if(topography) then
          do kk = 1, nz
             do jj = 1, ny
                do ii = 1, nx
                   
                   var(ii,jj,kk, iVart) = rho(ii,jj,kk) * dz_tr * heightTFC(ii, jj, kk)

                end do
             end do
          end do
       else
          do kk = 1, nz
             var(:,:,kk, iVart) = rho(:,:,kk) * dz_tr * (z(kk) -z(1))
          end do
       end if

       

    case( "like_density")

       do kk = 1,nz
          do jj = 1,ny
             do ii = 1,nx
                if (fluctuationMode) then
                   if (topography) then
                      var(ii,jj,kk,iVart) = var(ii,jj,kk,1) + rhoStratTFC(ii,jj,kk)
                   else
                      var(ii,jj,kk,iVart) = var(ii,jj,kk,1) + rhoStrat(kk)
                   end if
                else
                   var(ii,jj,kk,iVart) = var(ii,jj,kk,1)
                end if
             end do
          end do
       end do
  
    case ("test_tracer")

       call random_number(var(:,:,:,1))
       var(:,:,:,iVart) = var(:,:,:,1)

       


    case default
       stop "tracer.f90: unkown tracer distribution setup"
    end select


  end subroutine setup_tracer

end module tracer_module

