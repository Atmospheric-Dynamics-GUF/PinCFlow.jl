module tracer_module

  use type_module
  use atmosphere_module

  implicit none

  private

  public :: setup_tracer

  contains

   subroutine setup_tracer (var)
      ! setup initial tracer mixing ratio distribution
      ! rho*chi = var(:, :, :, iVarT)
      ! 
      ! initial large-scale tracer mixing ratio saved in 
      ! initialtracer(:, :, :) or 
      ! initialtracerrho(:, :, :) (multiplied by rho)

      real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
            & intent(inout) :: var
      real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz)::rho, tracerprime

      integer :: ii , jj , kk

      select case (tracerSetup)
      
         case( "increase_in_z_tracer" )
            ! chi = alphaTracer * z
            
            ! wavepacket initial tracer chi = alphaTracer*z + tracerprime
            if (testCase == 'wavePacket') then
               tracerprime = var(:, :, :, iVarT)
            end if
            
            ! determine total density
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

            if(topography) then
               do kk = 1, nz
                  do jj = 1, ny
                     do ii = 1, nx
                        if (testCase == 'wavePacket') then
                           var(ii, jj, kk, iVarT) = rho(ii,jj,kk) & 
                              *(tracerprime(ii, jj, kk) & 
                                 + alphaTracer * heightTFC(ii, jj, kk))
                        else
                           var(ii, jj, kk, iVarT) = rho(ii,jj,kk) &
                              * alphaTracer * heightTFC(ii, jj, kk)
                        end if
                        initialtracer(ii, jj, kk) = &
                              alphaTracer * heightTFC(ii, jj, kk)
                        initialtracerrho(ii, jj, kk) = &
                              alphaTracer * heightTFC(ii, jj, kk) &
                              * rho(ii, jj, kk)
                     end do
                  end do
               end do
            else
               do kk = 1, nz
                  if (testCase == 'wavePacket') then
                     var(:, :, kk, iVarT) = rho(:, :, kk) &
                        * (tracerprime(:, :, kk) + alphaTracer * (z(kk)-z(1)))
                  else
                     var(:, :, kk, iVarT) = rho(:, :, kk) &
                        * alphaTracer * (z(kk) -z(1))
                  end if
                  initialtracer(:, :, kk) = alphaTracer * (z(kk) -z(1))
                  initialtracerrho(:, :, kk) = alphaTracer * (z(kk) -z(1)) &
                                                * rho(:, :, kk)
               end do
            end if

         case default
            stop "tracer.f90: unkown tracer distribution setup"
      end select


   end subroutine setup_tracer

end module tracer_module

