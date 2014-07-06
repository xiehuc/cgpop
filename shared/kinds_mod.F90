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
!! Default numerical data types for all common data
!! types like integer, character, logical, real4 and real8.
!<
 module kinds_mod

! !USES:
!  uses no other modules

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer, parameter, public ::          &
      char_len  = 100                    ,&
      log_kind  = kind(.true.)           ,&
      i4        = selected_int_kind(6)   ,&
      i8        = selected_int_kind(13)  ,&
      r4        = selected_real_kind(6)  ,&
      r8        = selected_real_kind(13)

!***********************************************************************

 end module kinds_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
