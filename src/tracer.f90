module tracer_module

  use type_module
  use atmosphere_module

  implicit none

  private

  public :: setup_tracer

  contains

  subroutine setup_tracer(var)
    ! setup initial tracer mixing ratio distribution
    ! rho*chi = var(:, :, :, iVarT)
    !
    ! initial large-scale tracer mixing ratio saved in
    ! initialtracer(:, :, :) or
    ! initialtracerrho(:, :, :) (multiplied by rho)

    type(var_type), intent(inout) :: var
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: rho, &
        &tracerprime

    integer :: ii, jj, kk

    select case(tracerSetup)

    case("alpha_z")
      ! chi = alphaTracer * z

      ! wavepacket initial tracer chi = alphaTracer*z + tracerprime
      if(testCase == 'wavePacket') then
        tracerprime = var%chi(:, :, :)
      end if

      ! determine total density
      do kk = 0, nz + 1
        do jj = 0, ny + 1
          do ii = 0, nx + 1
            if(topography) then
              rho(ii, jj, kk) = (var%rho(ii, jj, kk) + rhoStratTFC(ii, jj, kk))
            else
              rho(ii, jj, kk) = var%rho(ii, jj, kk) + rhoStrat(kk)
            end if
          end do
        end do
      end do

      if(topography) then
        do kk = 0, nz + 1
          do jj = 0, ny + 1
            do ii = 0, nx + 1
              if(testCase == 'wavePacket') then
                var%chi(ii, jj, kk) = rho(ii, jj, kk) * (tracerprime(ii, jj, &
                    &kk) + alphaTracer * zTFC(ii, jj, kk))
              else
                var%chi(ii, jj, kk) = rho(ii, jj, kk) * alphaTracer &
                    &* zTFC(ii, jj, kk)
              end if
              initialtracer(ii, jj, kk) = alphaTracer * zTFC(ii, jj, kk)
              initialtracerrho(ii, jj, kk) = alphaTracer * zTFC(ii, jj, &
                  &kk) * rho(ii, jj, kk)
            end do
          end do
        end do
      else
        do kk = 0, nz + 1
          if(testCase == 'wavePacket') then
            var%chi(:, :, kk) = rho(:, :, kk) * (tracerprime(:, :, kk) &
                &+ alphaTracer * (z(kk) - z(1)))
          else
            var%chi(:, :, kk) = rho(:, :, kk) * alphaTracer * (z(kk) - z(1))
          end if
          initialtracer(:, :, kk) = alphaTracer * (z(kk) - z(1))
          initialtracerrho(:, :, kk) = alphaTracer * (z(kk) - z(1)) * rho(:, &
              &:, kk)
        end do
      end if

    case default
      stop "tracer.f90: unkown tracer distribution setup"
    end select

  end subroutine setup_tracer

end module tracer_module
