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

module FDfft 

  use, intrinsic :: iso_c_binding
  use FDistribution
  implicit none

  include 'fftw3.f03'

  private
  type dfft_type
    private
    type(C_ptr) :: object = C_NULL_ptr
  end type dfft_type

  interface

    function C_Dfft__new(dist) result(this) bind(C, name="Dfft__new")
      import
      type(C_ptr) :: this
      type(C_ptr), value :: dist
    end function C_Dfft__new 

    subroutine C_Dfft__makePlans(this, fo, fs, bi, bs, flag) bind(C, name="Dfft__makePlans")
      import
      type(C_ptr), value :: this
      type(C_ptr), value :: fo, fs, bi, bs
      integer(kind=C_size_t), value :: flag
    end subroutine C_Dfft__makePlans

    subroutine C_Dfft__forward(this, ain) bind(C, name="Dfft__forward")
      import
      type(C_ptr), value :: this
      type(C_ptr), value :: ain
    end subroutine C_Dfft__forward

    subroutine C_Dfft__backward(this, aout) bind(C, name="Dfft__backward")
      import
      type(C_ptr), value :: this
      type(C_ptr), value :: aout
    end subroutine C_Dfft__backward

    function C_Dfft__global_size(this) result(n) bind(C, name="Dfft__global_size")
      import
      type(C_ptr), value :: this
      integer(kind=C_size_t) :: n
    end function C_Dfft__global_size

    function C_Dfft__local_size(this) result(n) bind(C, name="Dfft__local_size")
      import
      type(C_ptr), value :: this
      integer(kind=C_size_t) :: n
    end function C_Dfft__local_size

    subroutine C_Dfft__global_ng(this, n) bind(C, name="Dfft__global_ng")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__global_ng

    subroutine C_Dfft__self_rspace(this, n) bind(C, name="Dfft__self_rspace")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__self_rspace

    subroutine C_Dfft__nproc_rspace(this, n) bind(C, name="Dfft__nproc_rspace")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__nproc_rspace

    subroutine C_Dfft__local_ng_rspace(this, n) bind(C, name="Dfft__local_ng_rspace")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__local_ng_rspace

    subroutine C_Dfft__self_kspace(this, n) bind(C, name="Dfft__self_kspace")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__self_kspace

    subroutine C_Dfft__nproc_kspace(this, n) bind(C, name="Dfft__nproc_kspace")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__nproc_kspace

    subroutine C_Dfft__local_ng_kspace(this, n) bind(C, name="Dfft__local_ng_kspace")
      import
      type(C_ptr), value :: this
      integer(C_int) :: n(3)
    end subroutine C_Dfft__local_ng_kspace

    function C_Dfft__parent_comm(this) result(comm) bind(C, name="Dfft__parent_comm")
      import
      type(C_ptr), value :: this
      integer :: comm
    end function C_Dfft__parent_comm

    function C_Dfft__cartcomm_rspace(this) result(comm) bind(C, name="Dfft__cartcomm_rspace")
      import
      type(C_ptr), value :: this
      integer :: comm
    end function C_Dfft__cartcomm_rspace

    function C_Dfft__cartcomm_kspace(this) result(comm) bind(C, name="Dfft__cartcomm_kspace")
      import
      type(C_ptr), value :: this
      integer :: comm
    end function C_Dfft__cartcomm_kspace

    subroutine C_Dfft__delete(this) bind(C, name="Dfft__delete")
      import
      type(C_ptr), value :: this
    end subroutine C_Dfft__delete

  end interface

  interface newDfft
    module procedure Dfft__new
  end interface newDfft

  interface makePlans
    module procedure Dfft__makePlans
  end interface makePlans

  interface forward
    module procedure Dfft__forward
  end interface forward

  interface backward
    module procedure Dfft__backward
  end interface backward

  interface globalSize
    module procedure Dfft__globalSize
  end interface globalSize

  interface localSize
    module procedure Dfft__localSize
  end interface localSize

  interface nglobal 
    module procedure Dfft__global_ng
  end interface nglobal

  interface selfRspace
    module procedure Dfft__self_rspace
  end interface selfRspace

  interface nprocRspace
    module procedure Dfft__nproc_rspace
  end interface nprocRspace

  interface nlocalRspace
    module procedure Dfft__local_ng_rspace
  end interface nlocalRspace

  interface selfKspace
    module procedure Dfft__self_kspace
  end interface selfKspace

  interface nprocKspace
    module procedure Dfft__nproc_kspace
  end interface nprocKspace

  interface nlocalKspace
    module procedure Dfft__local_ng_kspace
  end interface nlocalKspace

  interface parentComm
    module procedure Dfft__parent_comm
  end interface parentComm

  interface cartCommRspace 
    module procedure Dfft__cartcomm_rspace
  end interface cartCommRspace 

  interface cartCommKspace
    module procedure Dfft__cartcomm_kspace
  end interface cartCommKspace

  interface delDfft
    module procedure Dfft__delete
  end interface delDfft

  public :: newDfft, makePlans, forward, backward, globalSize, localSize, nglobal, selfRspace, nprocRspace, nlocalRspace, &
                selfKspace, nprocKspace, nlocalKspace, parentComm, cartCommRspace, cartCommKspace, delDfft, dfft_type

contains

  subroutine Dfft__new(this, dist) 
    type(dfft_type), intent(out) :: this
    type(dist_type), intent(in) :: dist
    this%object = C_Dfft__new(dist%object)
  end subroutine Dfft__new

  subroutine Dfft__makePlans(this, fo, fs, bi, bs, opt_flag)
    type(dfft_type), intent(in) :: this
    type(C_ptr), intent(in) :: fo, fs, bi, bs
    integer, intent(in), optional :: opt_flag
    integer :: flag
    if (present(opt_flag)) then
      flag = opt_flag
    else
      flag = FFTW_MEASURE
    endif
    call C_Dfft__makePlans(this%object, fo, fs, bi, bs, int(flag,C_size_t))
  end subroutine Dfft__makePlans

  subroutine Dfft__forward(this, ain)
    type(dfft_type), intent(in) :: this
    type(C_ptr), intent(inout) :: ain
    call C_Dfft__forward(this%object, ain)
  end subroutine

  subroutine Dfft__backward(this, aout)
    type(dfft_type), intent(in) :: this
    type(C_ptr), intent(inout) :: aout
    call C_Dfft__backward(this%object, aout)
  end subroutine

  integer function Dfft__globalSize(this)
    type(dfft_type), intent(in) :: this
    integer(C_size_t) :: n
    n = C_Dfft__global_size(this%object)
    Dfft__globalSize = int(n)
  end function Dfft__globalSize

  integer function Dfft__localSize(this)
    type(dfft_type), intent(in) :: this
    integer(C_size_t) :: n
    n = C_Dfft__local_size(this%object)
    Dfft__localSize = int(n) 
  end function Dfft__localSize

  subroutine Dfft__global_ng(this, self)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: self(3)
    integer(C_int) :: n(3)
    call C_Dfft__global_ng(this%object, n)
    self = int(n)
  end subroutine Dfft__global_ng

  subroutine Dfft__self_rspace(this, self)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: self(3)
    integer(C_int) :: n(3)    
    call C_Dfft__self_rspace(this%object, n)
    self = int(n)
  end subroutine Dfft__self_rspace

  subroutine Dfft__nproc_rspace(this, ng)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: ng(3)
    integer(C_int) :: n(3)
    call C_Dfft__nproc_rspace(this%object, n)
    ng = int(n)
  end subroutine Dfft__nproc_rspace

  subroutine Dfft__local_ng_rspace(this, ng)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: ng(3)
    integer(C_int) :: n(3)
    call C_Dfft__local_ng_rspace(this%object, n)
    ng = int(n)
  end subroutine Dfft__local_ng_rspace

  subroutine Dfft__self_kspace(this, self)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: self(3)
    integer(C_int) :: n(3)
    call C_Dfft__self_kspace(this%object, n)
    self = int(n)
  end subroutine Dfft__self_kspace

  subroutine Dfft__nproc_kspace(this, ng)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: ng(3)
    integer(C_int) :: n(3)
    call C_Dfft__nproc_kspace(this%object, n)
    ng = int(n)
  end subroutine Dfft__nproc_kspace

  subroutine Dfft__local_ng_kspace(this, ng)
    type(dfft_type), intent(in) :: this
    integer, intent(out) :: ng(3)
    integer(C_int) :: n(3)
    call C_Dfft__local_ng_kspace(this%object, n)
    ng = int(n)
  end subroutine Dfft__local_ng_kspace

  integer function Dfft__parent_comm(this)
    type(dfft_type), intent(in) :: this
    Dfft__parent_comm = C_Dfft__parent_comm(this%object)
  end function Dfft__parent_comm

  integer function Dfft__cartcomm_rspace(this)
    type(dfft_type), intent(in) :: this
    Dfft__cartcomm_rspace = C_Dfft__cartcomm_rspace(this%object)
  end function Dfft__cartcomm_rspace

  integer function Dfft__cartcomm_kspace(this)
    type(dfft_type), intent(in) :: this
    Dfft__cartcomm_kspace = C_Dfft__cartcomm_kspace(this%object)
  end function Dfft__cartcomm_kspace

  subroutine Dfft__delete(this)
    type(dfft_type), intent(inout) :: this
    call C_Dfft__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine Dfft__delete

end module FDfft 

