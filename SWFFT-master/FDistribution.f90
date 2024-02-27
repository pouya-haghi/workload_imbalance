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

module FDistribution

  use, intrinsic :: iso_c_binding, only: C_int, C_bool, C_ptr, C_NULL_ptr
  use mpi
  implicit none

  private
  type dist_type
    type(C_ptr) :: object = C_NULL_ptr
  end type dist_type

  interface

    function C_Distribution__new(comm, n, debug) result(this) bind(C, name="Distribution__new")
      import
      type(C_ptr) :: this
      integer :: comm
      integer(C_int) :: n(*)
      logical(C_bool), value :: debug
    end function C_Distribution__new 

    function C_Distribution__Cart_3D(this) result(comm) bind(C, name="Distribution__Cart_3D")
      import
      type(C_ptr), value :: this
      integer :: comm
    end function C_Distribution__Cart_3D

    subroutine C_Distribution__delete(this) bind(C, name="Distribution__delete")
      import
      type(C_ptr), value :: this
    end subroutine C_Distribution__delete

  end interface

  interface newDistribution
    module procedure Distribution__new
  end interface newDistribution

  interface cart3D
    module procedure Distribution__Cart_3D
  end interface cart3D

  interface delDistribution
    module procedure Distribution__delete 
  end interface delDistribution

  public :: newDistribution, cart3D, delDistribution, dist_type

contains

  subroutine Distribution__new(this, comm, n, opt_debug)
    type(dist_type), intent(out) :: this
    integer, intent(in) :: comm
    integer, intent(in) :: n(*)
    logical, intent(in), optional :: opt_debug
    logical :: debug
    if (present(opt_debug)) then
      debug = opt_debug
    else
      debug = .false.
    endif
    this%object = C_Distribution__new(comm, n, logical(debug,C_bool))
  end subroutine Distribution__new

  integer function Distribution__Cart_3D(this)
    type(dist_type), intent(in) :: this
    Distribution__Cart_3D = C_Distribution__Cart_3D(this%object)
  end function Distribution__Cart_3D 

  subroutine Distribution__delete(this)
    type(dist_type), intent(inout) :: this
    call C_Distribution__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine Distribution__delete

end module FDistribution

