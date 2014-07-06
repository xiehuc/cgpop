!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!>
!! Routine to provide a graceful means of exiting from POP when
!! encountering an error.
!<
 module exit_mod

! !USES:

   use kinds_mod, only: i4
   use communicate, only: my_task, master_task, &
	 abort_message_environment, exit_message_environment
   use constants, only: blank_fmt, delim_fmt

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: exit_POP

! !DEFINED PARAMETERS:

   integer (i4), parameter, public :: &
      sigExit  =  0,    &! signal for normal exit
      sigAbort = -1      ! signal for aborting (exit due to error)

!***********************************************************************

 contains

!***********************************************************************

!>
!! This routine prints a message, exits any message environment
!! and cleans up before stopping
!!
!! @param exit_mode        method for exiting (normal exit or abort)
!! @param exit_message     message to print before stopping
!<
 subroutine exit_POP(exit_mode, exit_message)

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
     exit_mode    

   character (*), intent(in) :: &
     exit_message

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: ierr  ! error flag

!-----------------------------------------------------------------------
!
!  print message - must use unit 6 in place of stdout to
!  prevent circular dependence with io module
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write (6,delim_fmt)
      write (6,blank_fmt)

      select case(exit_mode)
      case(sigExit)
         write (6,'(a14)') 'POP exiting...'
      case(sigAbort)
         write (6,'(a15)') 'POP aborting...'
      case default
         write (6,'(a37)') 'POP exiting with unknown exit mode...'
      end select

      write (6,*) exit_message
      write (6,blank_fmt)
      write (6,delim_fmt)
   endif

!-----------------------------------------------------------------------
!
!  exit or abort the message-passing environment if required
!
!-----------------------------------------------------------------------

   select case(exit_mode)
   case(sigExit)
      call exit_message_environment(ierr)
   case(sigAbort)
      call abort_message_environment(ierr)
   case default
   end select

!-----------------------------------------------------------------------
!
!  now we can stop
!
!-----------------------------------------------------------------------

   stop

!-----------------------------------------------------------------------

 end subroutine exit_POP

!***********************************************************************

 end module exit_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
