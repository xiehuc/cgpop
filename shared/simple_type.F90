!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This data structure describes the how blocks are distributed accross tasks.
!<
module simple_type

  use kinds_mod, only: i4
! !PUBLIC TYPES:

   !>
   !! Distribution data type (describes how blocks map to processes)
   !<
   type, public :: distrb  ! 

      !>
      !! number of processors in this dist
      !<
      integer (i4) :: nprocs

      !>
      !! communicator to use in this dist
      !<
      integer (i4) :: communicator

      !>
      !! processor location for this block
      !<
      integer (i4), dimension(:), pointer :: proc

      !>
      !! block position in local array on proc
      !<
      integer (i4), dimension(:), pointer :: local_block         
   end type

end module simple_type 

