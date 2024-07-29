module opt_field_mod

  implicit none
  !integer, parameter :: numBscVar = 8 ! number of "basic" PincFlow fields

  contains

  subroutine set_opt_field

    use type_module, ONLY:include_tracer, include_ice2, master, nVar, nBscVar, &
        iVarIce, nVarIce, inN, inQ, inQv, iVarT, iVarO, include_testoutput, &
        model, iVarP

    integer :: numOfld, indM1, numOfldVar, allocstat !auxillary variables

    nBscVar = nVar !

    if(model == "compressible") then
      nVar = nVar + 1
      iVarP = nVar
    end if

    if(include_tracer) then
      nVar = nVar + 1
      iVarT = nVar
    end if

    if(include_ice2) then
      nVarIce = 3
      nVar = nVar + nVarIce
      inN = nVar - 2
      inQ = nVar - 1
      inQv = nVar

      allocate(iVarIce(nVarIce), stat = allocstat)
      if(allocstat /= 0) stop "set_dim_opt_field.f90: Could not allocate &
          iVarIce."

      iVarIce = (/inN, inQ, inQv/)
    end if

    if(include_testoutput) then
      nVar = nVar + 3
      iVarO = nVar - 2
    end if

  end subroutine set_opt_field

  subroutine read_nml_opt_field

    use type_module, ONLY:master, include_tracer, include_ice2, iceList2, &
        tracerList, nVar, nBscVar, varOut, varIn, iVarIce, nVarIce, inN, inQ, &
        inQv, iVarT, varOut, include_testoutput, iVarO, model, iVarP, varOut, &
        varIn
    integer i, j, del

    if(model == "compressible") then
      varOut(iVarP) = 0
      varIn(iVarP) = 0
    end if

    if(include_ice2) then
      ! read ice physics parametrization
      read(unit = 10, nml = iceList2)

      do j = 1, nVarIce
        i = iVarIce(j)
        varOut(i) = 1
        varIn(i) = 1
      end do
    end if

    if(include_tracer) then
      read(unit = 10, nml = tracerList)
      varOut(iVarT) = 1
      varIn(iVarT) = 1
    end if

    if(include_testoutput) then
      varOut(iVarO:iVarO + 2) = 1
      varIn(iVarO:iVarO + 2) = 1
    end if

  end subroutine read_nml_opt_field

  subroutine write_index_opt_field

    use type_module, ONLY:master, include_tracer, include_ice2, master, nVar, &
        nBscVar, iVarIce, nVarIce, inN, inQ, inQv, iVarT, varOut, &
        include_testoutput, iVarO, varOut, varIn
    integer i, j, del

    if(master) then
      open(555, file = 'index_opt_field.py')

      !check if all basic fields are output
      j = 0
      do i = 1, nBscVar
        if(varOut(i) == 1) then
          j = j + 1
        end if
      end do

      ! python index starts at 0 ( + 1 shift in del)
      del = nBscVar - j + 1

      if(include_tracer) then
        write(555, "(a,i2.1)") 'iTr = ', iVarT - del
      end if

      if(include_ice2) then
        write(555, "(a,i2.1)") 'iIce1 = ', inN - del
        write(555, "(a,i2.1)") 'iIce2 = ', inQ - del
        write(555, "(a,i2.1)") 'iIce3 = ', inQv - del
      end if

      if(include_testoutput) then
        write(555, "(a,i2.1)") 'iO1 = ', iVarO - del
        write(555, "(a,i2.1)") 'iO2 = ', iVarO - del + 1
        write(555, "(a,i2.1)") 'iO3 = ', iVarO - del + 2
      end if

      close(555)
    end if

  end subroutine write_index_opt_field

end module opt_field_mod

module opt_field_old_mod

  implicit none
  integer, parameter :: MaxNumVarOptFld = 5
  integer, parameter :: numBscVar = 8 ! number of "basic" PincFlow fields

  type :: opt_field

    integer :: numVar
    logical :: output
    integer :: indVarPF(MaxNumVarOptFld)

    !contains
    !  procedure :: init_opt_field
  end type opt_field

  type(opt_field) :: ofTracer, ofIce, ofSomething

  contains

  subroutine init_opt_field(ofld, nvar, out, numOfld, indOfldVarM1)

    class(opt_field), intent(inout) :: ofld
    integer :: nvar !number of variables for opt_field
    !tracer : 1; ice2 : 3 (n,q,q_v)
    logical :: out ! True if output required

    integer :: numOfld, indOfldVarM1 !auxillary counters
    integer :: i

    if(nvar .gt. MaxNumVarOptFld) then
      print *, 'opt_field%numVar > MaxNumVarOptFld !'
      print *, 'stop'
      stop
    end if

    ofld%numVar = nvar
    ofld%output = out

    ofld%indVarPF(:) = 0 !init
    do i = 1, nvar
      ofld%indVarPF(i) = indOfldVarM1 + i
    end do

    numOfld = numOfld + 1
    indOfldVarM1 = indOfldVarM1 + nvar

    print *, 'Init opt_field'

  end subroutine init_opt_field

  subroutine print_opt_field(ofld)
    class(opt_field), intent(in) :: ofld
    integer :: i

    print *, 'produce output ', ofld%output
    print *, 'number variables ', ofld%numVar

    do i = 1, ofld%numVar
      print *, 'index of variable ', i, ' is', ofld%indVarPF(i)
    end do
  end subroutine print_opt_field

  subroutine set_opt_field_old

    use type_module, ONLY:include_tracer, include_ice2, master, nVar
    !logical, parameter :: include_tracer=.true., include_ice2=.false.
    integer :: numOfld, indOfldVarM1 !auxillary counters
    !class(opt_field), intent(inout) :: ofTracer, ofIce, ofSomething

    open(555, file = 'index_opt_field.py')

    numOfld = 0 ! first no optional fields
    !numBscVar = nVar
    indOfldVarM1 = numBscVar

    if(include_tracer) then
      call init_opt_field(ofTracer, 1, .true., numOfld, indOfldVarM1)
      nVar = nVar + ofTracer%numVar
      if(master) then
        !call print_opt_field(ofTracer)
        write(555, "(a,i2.1)") 'iTr = ', ofTracer%indVarPF(1)
      end if
    end if

    if(include_ice2) then
      call init_opt_field(ofIce, 3, .true., numOfld, indOfldVarM1)
      nVar = nVar + ofIce%numVar
      if(master) then
        !call print_opt_field(ofIce)
        write(555, "(a,i2.1)") 'iIce1 = ', ofIce%indVarPF(1)
        write(555, "(a,i2.1)") 'iIce2 = ', ofIce%indVarPF(2)
        write(555, "(a,i2.1)") 'iIce3 = ', ofIce%indVarPF(3)
      end if
    end if

    close(555)

  end subroutine set_opt_field_old

  subroutine write_index_opt_field_old

    use type_module, ONLY:include_tracer, include_ice2, master, nVar, varOut
    !logical, parameter :: include_tracer=.true., include_ice2=.false.
    integer :: numOfld, indOfldVarM1 !auxillary counters
    !class(opt_field), intent(inout) :: ofTracer, ofIce, ofSomething
    integer i, j, numBsc, del

    if(master) then
      open(555, file = 'index1_opt_field.py')

      !check if all basic fields are output
      j = 0
      do i = 1, numBscVar
        if(varOut(i) == 1) then
          j = j + 1
        end if
      end do

      del = numBscVar - j - 1

      if(include_tracer) then
        write(555, "(a,i2.1)") 'iTr = ', ofTracer%indVarPF(1) - del
      end if

      if(include_ice2) then
        write(555, "(a,i2.1)") 'iIce1 = ', ofIce%indVarPF(1) - del
        write(555, "(a,i2.1)") 'iIce2 = ', ofIce%indVarPF(2) - del
        write(555, "(a,i2.1)") 'iIce3 = ', ofIce%indVarPF(3) - del
      end if

      close(555)
    end if

  end subroutine write_index_opt_field_old

end module opt_field_old_mod
