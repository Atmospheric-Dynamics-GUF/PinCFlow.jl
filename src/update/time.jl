struct Time{A<:Integer,B<:Vector{<:AbstractFloat}}
    nstages::A
    alphark::B
    betark::B
    stepfrac::B
end

function Time()

  # Set Runge-Kutta parameters.
  nstages = 3
  alphark = [0.0, -5.0 / 9.0, -153.0 / 128.0]
  betark = [1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0]
  stepfrac = [1.0 / 3.0, 5.0 / 12.0, 1.0 / 4.0]

  # Return a Time instance.
  return Time(nstages, alphark, betark, stepfrac)
end
