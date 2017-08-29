program timeoutput   

  implicit none
  
  
  integer    :: days, hours, mins, secs
  real       :: timeVar, cpuTime
  character (len = 40)  :: cpuTimeChar
  
  cpuTime = 129382.0
  
  ! calc cpu-time in days/hours/minutes/seconds
  timeVar = cpuTime   
  days = floor(timeVar / 86400.0)
  timeVar = timeVar - 86400.0 * days
  hours = floor(timeVar / 3600.0)
  timeVar = timeVar - 3600.0 * hours
  mins = floor(timeVar / 60.0)
  timeVar = timeVar - 60.0 * mins
  secs = int(timeVar)


  ! write tecplot header
  write(unit=cpuTimeChar, fmt="(3(i2,a),i2)") &
       & days, " ", hours, ":",mins,":",secs

  write(*,fmt="(a25,a25)")  "CPU time = ", cpuTimeChar



end program timeoutput

