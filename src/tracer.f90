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
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz)::rho, tracerprime

    real :: tracer0 ! initial max. tracer ampl.
    real :: xt0 , yt0 , zt0
    real :: sigtx, sigty, sigtz
    real :: k0t , m0t
    real :: alpha
    real :: zcart
    real :: dz_tr

    integer :: ii , jj , kk

    select case (tracerSetup)
   
    case( "increase_in_z_tracer" )

       if (testCase == 'wavePacket') then
          tracerprime = var(:, :, :, iVart)
          !tracerprime = - tracerprime
       end if
       
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

       alpha = alphaTracer!/ (lz(1) - lz(0))

       if(topography) then
          do kk = 1, nz
             do jj = 1, ny
                do ii = 1, nx
                   if (testCase == 'wavePacket') then
                      var(ii, jj, kk, iVart) = rho(ii,jj,kk)*(tracerprime(ii, jj, kk) + alpha * heightTFC(ii, jj, kk))
                   else
                      var(ii,jj,kk, iVart)   = rho(ii,jj,kk)*alpha * heightTFC(ii, jj, kk)
                   end if
                   initialtracer(ii, jj, kk) = alpha * heightTFC(ii, jj, kk) ! rho(ii,jj,kk)*
                end do
             end do
          end do
       else
          do kk = 1, nz
             if (testCase == 'wavePacket') then
                var(:, :, kk, iVart) = rho(:, :, kk) * (tracerprime(:, :, kk) + alpha * (z(kk)-z(1)))
             else
                var(:,:,kk, iVart) = rho(:,:,kk) * alpha * (z(kk) -z(1))
             end if
             initialtracer(:, :, kk) = alpha * (z(kk) -z(1)) ! rho(:,:,kk) * 
          end do
       end if

    case( "quadratic_increase" )

         if (testCase == 'wavePacket') then
            tracerprime = var(:, :, :, iVart)
            !tracerprime = - tracerprime
         end if
         
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
  
         alpha = alphaTracer!/ (lz(1) - lz(0))
  
         if(topography) then
            do kk = 1, nz
               do jj = 1, ny
                  do ii = 1, nx
                     if (testCase == 'wavePacket') then
                        var(ii, jj, kk, iVart) = rho(ii,jj,kk)*(tracerprime(ii, jj, kk) + alpha * heightTFC(ii, jj, kk))**2.0
                     else
                        var(ii,jj,kk, iVart)   = rho(ii,jj,kk)*alpha * heightTFC(ii, jj, kk)**2.0
                     end if
                     initialtracer(ii,jj,kk)  = alpha * heightTFC(ii, jj, kk)**2.0
                  end do
               end do
            end do
         else
            do kk = 1, nz
               if (testCase == 'wavePacket') then
                  var(:, :, kk, iVart) = rho(:, :, kk) * (tracerprime(:, :, kk) + alpha * (z(kk)-z(1))**2.0)
               else
                  var(:,:,kk, iVart) = rho(:,:,kk) * alpha * (z(kk) -z(1))**2.0
               end if
               initialtracer(:,:,kk) = alpha * (z(kk) -z(1))**2.0
            end do
         end if

 
    
    case default
       stop "tracer.f90: unkown tracer distribution setup"
    end select

    !initialtracer = var(:, :, :, iVart)


  end subroutine setup_tracer

end module tracer_module

