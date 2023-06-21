program test
integer ::  v(3) 
integer a

a = 3
v =(/1, 2, 3/)

if (any(a == v)) then

print*, v

end if

end program test
