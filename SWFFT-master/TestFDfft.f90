!                  Copyright (C) 2017, UChicago Argonne, LLC
!                             All Rights Reserved
!
!            Hardware/Hybrid Cosmology Code (HACC), Version 1.0
!
!  Salman Habib, Adrian Pope, Hal Finkel, Nicholas Frontiere, Katrin Heitmann,
!       Vitali Morozov, Jeffrey Emberson, Thomas Uram, Esteban Rangel
!                         (Argonne National Laboratory)
!
!   David Daniel, Patricia Fasel, Chung-Hsing Hsu, Zarija Lukic, James Ahrens
!                       (Los Alamos National Laboratory)
!
!                                George Zagaris
!                                  (Kitware)
!
!                             OPEN SOURCE LICENSE
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!   1. Redistributions of source code must retain the above copyright notice,
!      this list of conditions and the following disclaimer. Software changes,
!      modifications, or derivative works, should be noted with comments and
!      the author and organization’s name.
!
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!
!   3. Neither the names of UChicago Argonne, LLC or the Department of Energy 
!      nor the names of its contributors may be used to endorse or promote 
!      products derived from this software without specific prior written 
!      permission.
!
!   4. The software and the end-user documentation included with the
!      redistribution, if any, must include the following acknowledgment:
!
!     "This product includes software produced by UChicago Argonne, LLC under
!      Contract No. DE-AC02-06CH11357 with the Department of Energy."
!
! *****************************************************************************
!                                DISCLAIMER
! THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND. NEITHER THE
! UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR 
! UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE
! ACCURARY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS,
! PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
! PRIVATELY OWNED RIGHTS.
!
! *****************************************************************************

program main

  use, intrinsic :: iso_c_binding
  use FDistribution
  use FDfft
  use mpi
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  include 'fftw3.f03'

  integer :: rank, nproc, ierr
  integer :: repetitions
  integer :: ng(3)
  integer :: omt
  real(8) :: t1, t2 

  !
  ! Initialize MPI communicators
  !

  call mpi_initialize

  !
  ! Assign user arguments
  ! 

  call read_user_arguments

  !
  ! Initialize fftw3 openmp threads if neccessary
  !

#ifdef _OPENMP
  ierr = fftw_init_threads()
  if (ierr == 0) then 
    write(*,*) "fftw_init_threads() failed!"
    call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
  endif
  omt = omp_get_max_threads()
  call fftw_plan_with_nthreads(omt)
  if (rank == 0) write(*,*) "Threads per process: ", omt
#endif

  !
  ! Testing subroutine 
  !

  t1 = MPI_Wtime()
  call test
  t2 = MPI_Wtime()
  if (rank == 0) write(*,*) " TEST TIME: ", t2-t1

  !
  ! Finalize MPI communicators and exit 
  !

  call MPI_Finalize(ierr)

contains

!! -------------------------------------------------------------------------- !!

subroutine mpi_initialize
    !
    ! Initializes MPI communication 
    !

    implicit none

    !
    ! Setup communications
    !

    call MPI_Init(ierr)
    if (ierr .ne. MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    if (ierr .ne. MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (ierr .ne. MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)

    return

end subroutine mpi_initialize

!! -------------------------------------------------------------------------- !!

subroutine read_user_arguments
  !
  ! Read in user supplied arguments
  !

  implicit none

  character(len=200) :: argument

  if (command_argument_count() < 2 .and. rank == 0) then
    call get_command_argument(0, argument)
    write(*,*) "USAGE: ", trim(argument), " <n_repetitions> <ngx> [ngy ngz]"
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  endif

  call get_command_argument(1, argument)
  read(argument, "(i4)") repetitions 

  call get_command_argument(2, argument)
  read(argument, "(i4)") ng(1)
  ng(2) = ng(1) ; ng(3) = ng(1)
  if (command_argument_count() > 2) then
    call get_command_argument(3, argument)
    read(argument, "(i4)") ng(2)
    call get_command_argument(4, argument)
    read(argument, "(i4)") ng(3)
  endif

  return

end subroutine read_user_arguments  

!! -------------------------------------------------------------------------- !!

subroutine test
  !
  ! Test the fft code.
  !

  implicit none

  type(dist_type) :: d
  type(dfft_type) :: dfft

  complex(C_double_complex), pointer :: a_arr(:), b_arr(:)
  type(C_ptr) :: pa, pb

  integer :: nlocal
  integer :: cartComm, cartRank
  integer :: i
  real(8) :: tstart, tstop
  real(8) :: zero, one, numg

  !! Instantiate distribution and dfft object
  call newDistribution(d, MPI_COMM_WORLD, ng) 
  call newDfft(dfft, d)

  !! Determine size of 1D local data on this rank
  nlocal = localSize(dfft)

  !! Use fftw allocation to obtain proper memory alignment
  pa = fftw_alloc_complex(int(nlocal,C_size_t))
  pb = fftw_alloc_complex(int(nlocal,C_size_t))

  !! Associate C pointer with the Fortran array
  call c_f_pointer(pa, a_arr, [nlocal])
  call c_f_pointer(pb, b_arr, [nlocal])

  !! Have Dfft make FFTW plans based on the allocated arrays
  call makePlans(dfft, pa, pb, pa, pb)

  !! Grab a copy of the 3D Cartesian communicator that Distribution is using
  cartComm = cart3D(d) 
  call MPI_Comm_rank(cartComm, cartRank, ierr)

  !! Write out hex representation of the three relevent numbers to the problem
  if (cartRank == 0) then
    zero = 0.0
    one  = 1.0
    numg = 1.0*ng(1)*ng(2)*ng(3)
    write(*,*)
    write(*,*) "Hex representations of double precision floats"
    write(*,fmt="(F18.3,A,Z18)") zero, " = ", zero
    write(*,fmt="(F18.3,A,Z18)") one, " = ", one 
    write(*,fmt="(F18.3,A,Z18)") numg, " = ", numg 
    write(*,*)
  endif

  !! Forward and backward transform a delta function many times
  do i = 1, repetitions

    if (cartRank == 0) write(*,*) "TESTING ", i-1 

    call assign_delta_function(dfft, a_arr, nlocal)

    tstart = MPI_Wtime()
    call forward(dfft, pa)
    tstop = MPI_Wtime()

    call check_kspace(dfft, a_arr)

    tstart = MPI_Wtime()
    call backward(dfft, pa)
    tstop = MPI_Wtime()

    call check_rspace(dfft, a_arr)

  enddo

  !! Free data arrays
  call fftw_free(pa)
  call fftw_free(pb)

  !! Free distributiona and dfft objects
  call delDfft(dfft)
  call delDistribution(d)

  return 

end subroutine test

!! -------------------------------------------------------------------------- !!

subroutine assign_delta_function(dfft, a, n)

  implicit none

  type(dfft_type), intent(in) :: dfft  
  complex(C_double_complex), intent(inout), dimension(:) :: a
  integer, intent(in) :: n  

  integer, dimension(3) :: self
  integer, dimension(3) :: local_ng

  integer :: local_indx
  integer :: i, j, k
  integer :: global_i, global_j, global_k 

  !! Determine location of rank in r-space
  call selfRspace(dfft, self)

  !! Determine local grid dimensions in r-space
  call nlocalRspace(dfft, local_ng)

  !! Fill in the delta function.
  !! NOTE: We are filling in one-dimensional memory using the row-major C convention.
  local_indx = 1
  do i = 1, local_ng(1)
    global_i = local_ng(1)*self(1) + i
    do j = 1, local_ng(2)
      global_j = local_ng(2)*self(2) + j
      do k = 1, local_ng(3)
        global_k = local_ng(3)*self(3) + k
        if (global_i == 1 .and. global_j == 1 .and. global_k == 1) then
          a(local_indx) = cmplx(1.,0.)
        else
          a(local_indx) = cmplx(0.,0.)
        endif
        local_indx = local_indx + 1
      enddo
    enddo
  enddo

end subroutine assign_delta_function

!! -------------------------------------------------------------------------- !!

subroutine check_kspace(dfft, a)

  implicit none

  type(dfft_type), intent(in) :: dfft
  complex(C_double_complex), intent(in), dimension(:) :: a

  real(8) :: LocalRealMin, LocalRealMax, LocalImagMin, LocalImagMax
  real(8) :: GlobalRealMin, GlobalRealMax, GlobalImagMin, GlobalImagMax
  real(8) :: re, im

  integer :: parComm, parRank
  integer :: i, nlocal

  nlocal = localSize(dfft)

  LocalRealMin = real(a(2))  ; LocalRealMax = real(a(2))
  LocalImagMin = aimag(a(2)) ; LocalImagMax = aimag(a(2))

  do i = 1, nlocal
    re = real(a(i))
    im = aimag(a(i))
    if (re < LocalRealMin) LocalRealMin = re
    if (re > LocalRealMax) LocalRealMax = re
    if (im < LocalImagMin) LocalImagMin = im
    if (im > LocalImagMax) LocalImagMax = im
  enddo

  parComm = parentComm(dfft)
  call MPI_Comm_rank(parComm, parRank, ierr)  

  call MPI_Allreduce(LocalRealMin, GlobalRealMin, 1, MPI_DOUBLE, MPI_MIN, parComm, ierr)
  call MPI_Allreduce(LocalRealMax, GlobalRealMax, 1, MPI_DOUBLE, MPI_MAX, parComm, ierr)
  call MPI_Allreduce(LocalImagMin, GlobalImagMin, 1, MPI_DOUBLE, MPI_MIN, parComm, ierr)
  call MPI_Allreduce(LocalImagMax, GlobalImagMax, 1, MPI_DOUBLE, MPI_MAX, parComm, ierr)

  if (parRank == 0) then
    write(*,*) 
    write(*,*) "k-space:"
    write(*,fmt="(A,F18.3,F18.3,A,Z18,Z18,A)") "real in [", GlobalRealMin, GlobalRealMax, "] = [", GlobalRealMin, GlobalRealMax, "]"
    write(*,fmt="(A,F18.3,F18.3,A,Z18,Z18,A)") "imag in [", GlobalImagMin, GlobalImagMax, "] = [", GlobalImagMin, GlobalImagMax, "]"
    write(*,*)
  endif

end subroutine check_kspace

!! -------------------------------------------------------------------------- !!

subroutine check_rspace(dfft, a)

  implicit none

  type(dfft_type), intent(in) :: dfft
  complex(C_double_complex), intent(in), dimension(:) :: a

  integer, dimension(3) :: self
  integer, dimension(3) :: local_ng

  real(8) :: LocalRealMin, LocalRealMax, LocalImagMin, LocalImagMax
  real(8) :: GlobalRealMin, GlobalRealMax, GlobalImagMin, GlobalImagMax
  real(8) :: re, im

  integer :: parComm, parRank
  integer :: global_i, global_j, global_k
  integer :: local_indx, i, j, k

  call selfRspace(dfft, self)
  call nlocalRspace(dfft, local_ng)

  LocalRealMin = real(a(2))  ; LocalRealMax = real(a(2))
  LocalImagMin = aimag(a(2)) ; LocalImagMax = aimag(a(2))

  parComm = parentComm(dfft)
  call MPI_Comm_rank(parComm, parRank, ierr)

  if (parRank == 0) then
    write(*,*)
    write(*,*) "r-space:"
  endif

  local_indx = 1
  do i = 1, local_ng(1)
    global_i = local_ng(1)*self(1) + i
    do j = 1, local_ng(2)
      global_j = local_ng(2)*self(2) + j
      do k = 1, local_ng(3)
        global_k = local_ng(3)*self(3) + k
        if (global_i == 1 .and. global_j == 1 .and. global_k == 1) then
            write(*,fmt="(A,F18.3,F18.3,A,Z18,Z18,A)") "a[0,0,0] = ", real(a(local_indx)), aimag(a(local_indx)), &
                                                                     "= (", real(a(local_indx)), aimag(a(local_indx)), ")"
        else
          re = real(a(local_indx))
          im = aimag(a(local_indx))
          if (re < LocalRealMin) LocalRealMin = re
          if (re > LocalRealMax) LocalRealMax = re
          if (im < LocalImagMin) LocalImagMin = im
          if (im > LocalImagMax) LocalImagMax = im
        endif
        local_indx = local_indx + 1
      enddo
    enddo
  enddo

  call MPI_Allreduce(LocalRealMin, GlobalRealMin, 1, MPI_DOUBLE, MPI_MIN, parComm, ierr)
  call MPI_Allreduce(LocalRealMax, GlobalRealMax, 1, MPI_DOUBLE, MPI_MAX, parComm, ierr)
  call MPI_Allreduce(LocalImagMin, GlobalImagMin, 1, MPI_DOUBLE, MPI_MIN, parComm, ierr)
  call MPI_Allreduce(LocalImagMax, GlobalImagMax, 1, MPI_DOUBLE, MPI_MAX, parComm, ierr)

  if (parRank == 0) then
    write(*,fmt="(A,F18.3,F18.3,A,Z18,Z18,A)") "real in [", GlobalRealMin, GlobalRealMax, "] = [", GlobalRealMin, GlobalRealMax, "]"
    write(*,fmt="(A,F18.3,F18.3,A,Z18,Z18,A)") "imag in [", GlobalImagMin, GlobalImagMax, "] = [", GlobalImagMin, GlobalImagMax, "]"
    write(*,*)
  endif

end subroutine check_rspace

!! -------------------------------------------------------------------------- !!

end program main

