module sizeof_module
   ! module is provided as a replacement for sizeof() and aquires the sizes of int, single and double types
   implicit none
   public

   ! test variables to get sizes in terms of blocks
   double precision :: doubletest = 1.0
   real :: realtest = 1.0
   real*4 :: real4test = 1.0
   integer :: inttest = 1

   integer :: sizeofdouble, sizeofreal, sizeofreal4, sizeofint

contains
   ! subroutine to get the size of a double in terms of block size
   subroutine getsize
      inquire(iolength = sizeofdouble) doubletest
      inquire(iolength = sizeofreal) realtest
      inquire(iolength = sizeofreal4) real4test
      inquire(iolength = sizeofint) inttest
   end subroutine getsize
end module sizeof_module
