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
!! I/O unit manager for tracking, assigning and reserving I/O unit numbers.
!!
!! There are three reserved I/O units set as parameters in this
!! module.  The default units for standard input (stdin), standard
!! output (stdout) and standard error (stderr).  These are currently
!! set as units 5,6,6, respectively as that is the most commonly
!! used among vendors. However, the user may change these if those
!! default units are conflicting with other models or if the
!! vendor is using different values.
!!
!! The maximum number of I/O units per node is currently set by
!! the parameter POP\_IOMaxUnits.
!<
 module IOUnitsMod

! !USES:

   use kinds_mod, only: i4, log_kind

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: IOUnitsGet,                &
             IOUnitsRelease,            &
             IOUnitsFlush

! !PUBLIC DATA MEMBERS:

   integer (i4), parameter, public :: &
      stdin  =  5,  &! reserved unit for standard input
      stdout =  6,  &! reserved unit for standard output
      stderr =  6,  &! reserved unit for standard error
      nmlin  = 10

   ! common formats for writing to stdout, stderr

   character (9), parameter, public ::   &
      delimFormat    = "(72('-'))",  &
      delimFormatNew = "(72('='))"

   character (5), parameter, public :: &
      blankFormat = "(' ')" 

   character (7), parameter, public :: &
      nml_filename = 'pop_in'  ! namelist input file name


!-----------------------------------------------------------------------
!
!  private io unit manager variables
!
!-----------------------------------------------------------------------

   integer (i4), parameter :: &
      IOUnitsMinUnits = 11,   & ! do not use unit numbers below this
      IOUnitsMaxUnits = 99      ! maximum number of open units

   logical (log_kind) :: &
      IOUnitsInitialized = .false.

   logical (log_kind), dimension(IOUnitsMaxUnits) :: &
      IOUnitsInUse       ! flag=.true. if unit currently open

!***********************************************************************

contains

!***********************************************************************

!>
!! This routine returns the next available i/o unit and marks it as
!! in use to prevent any later use.
!! Note that {\em all} processors must call this routine even if only
!! the master task is doing the i/o.  This is necessary insure that
!! the units remain synchronized for other parallel I/O functions.
!!
!! @param iunit  next free i/o unit
!<
 subroutine IOUnitsGet(iunit)

   integer (i4), intent(out) :: &
      iunit

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: n  ! dummy loop index

   logical (log_kind) :: alreadyInUse

!-----------------------------------------------------------------------
!
!  check to see if units initialized and initialize if necessary
!
!-----------------------------------------------------------------------

   if (.not. IOUnitsInitialized) then
      IOUnitsInUse = .false.
      IOUnitsInUse(stdin) = .true.
      IOUnitsInUse(stdout) = .true.
      IOUnitsInUse(stderr) = .true.

      IOUnitsInitialized = .true.
   endif

!-----------------------------------------------------------------------
!
!  find next free unit
!
!-----------------------------------------------------------------------

   srch_units: do n=IOUnitsMinUnits, IOUnitsMaxUnits
      if (.not. IOUnitsInUse(n)) then   ! I found one, I found one

         !*** make sure not in use by library or calling routines
         INQUIRE (unit=n,OPENED=alreadyInUse)

         if (.not. alreadyInUse) then
            iunit = n        ! return the free unit number
            IOUnitsInUse(iunit) = .true.  ! mark iunit as being in use
            exit srch_units
         else
            !*** if inquire shows this unit in use, mark it as
            !***    in use to prevent further queries
            IOUnitsInUse(n) = .true.
         endif
      endif
   end do srch_units

   if (iunit > IOUnitsMaxUnits) stop 'IOUnitsGet: No free units'

!-----------------------------------------------------------------------
!EOC

 end subroutine IOUnitsGet

!***********************************************************************

!>
!! This routine releases an i/o unit (marks it as available).
!! Note that {\em all} processors must call this routine even if only
!! the master task is doing the i/o.  This is necessary insure that
!! the units remain synchronized for other parallel I/O functions.
!!
!! @param  iunit   i/o unit to be released
!<
 subroutine IOUnitsRelease(iunit)

! !INPUT PARAMETER:

   integer (i4), intent(in) :: &
      iunit                    ! i/o unit to be released

!-----------------------------------------------------------------------
!
!  check for proper unit number
!
!-----------------------------------------------------------------------

   if (iunit < 1 .or. iunit > IOUnitsMaxUnits) then
      stop 'IOUnitsRelease: bad unit'
   endif

!-----------------------------------------------------------------------
!
!  mark the unit as not in use
!
!-----------------------------------------------------------------------

   IOUnitsInUse(iunit) = .false.  !  that was easy...

!-----------------------------------------------------------------------

 end subroutine IOUnitsRelease

!***********************************************************************

!>
!! This routine enables a user to flush the output from an IO unit
!! (typically stdout) to force output when the system is buffering
!! such output.  Because this system function is system dependent,
!! we only support this wrapper and users are welcome to insert the
!! code relevant to their local machine.  In the case where the CCSM
!! libraries are available, the shared routine for sys flush can be
!! used (and is provided here under a preprocessor option).
!!
!! @param iunit  i/o unit to be flushed
!<
 subroutine IOUnitsFlush(iunit)

! !INPUT PARAMETER:

   integer (i4), intent(in) :: &
      iunit 

!-----------------------------------------------------------------------
!
!  insert your system code here
!
!-----------------------------------------------------------------------

 end subroutine IOUnitsFlush

!***********************************************************************

 end module IOUnitsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
