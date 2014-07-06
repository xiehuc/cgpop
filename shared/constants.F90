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
!! Physical and numerical constants used throughout the Parallel Ocean
!! Program.
!<
 module constants

   use kinds_mod, only: r8, i4, r4, char_len

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_constants

! !DEFINED PARAMETERS:

   ! numbers
   
   real (r8), parameter, public :: &
      c0     =    0.0_r8   ,&
      c1     =    1.0_r8   ,&
      c2     =    2.0_r8   ,&
      c3     =    3.0_r8   ,&
      c4     =    4.0_r8   ,&
      c5     =    5.0_r8   ,&
      c8     =    8.0_r8   ,&
      c10    =   10.0_r8   ,&
      c16    =   16.0_r8   ,&
      c1000  = 1000.0_r8   ,&
      c10000 =10000.0_r8   ,&
      c1p5   =    1.5_r8   ,&
      p33    = c1/c3       ,&
      p5     = 0.500_r8    ,&
      p25    = 0.250_r8    ,&
      p125   = 0.125_r8    ,&
      p001   = 0.001_r8    ,&
      eps    = 1.0e-10_r8  ,&
      eps2   = 1.0e-20_r8  ,&
      bignum = 1.0e+30_r8
   
   integer, public, parameter :: &
        ALG_MPI2S_1D              = 1,  &
        ALG_MPI2S_2D              = 2,  &
        ALG_CAF_MULTI_PULL_1D     = 3,  &
        ALG_CAF_SINGLE_PULL_1D    = 4,  &
        ALG_CAF_SINGLE_PUSH_1D    = 5,  &
        ALG_CAF_SINGLE_PUSH_2D    = 6,  &
        ALG_CAF_SINGLE_PULL_2D    = 7,  &
        ALG_MPI1S_MULTI_PULL_1D   = 8,  &
        ALG_MPI1S_SINGLE_PULL_1D  = 9,  &
        ALG_MPI1S_SINGLE_PUSH_1D  = 10, &
        ALG_MPI1S_2D              = 11, &
        ALG_HYBRID_MULTI_PULL_1D  = 12, &
        ALG_HYBRID_SINGLE_PULL_1D = 13, &
        ALG_HYBRID_SINGLE_PUSH_1D = 14
   
   integer(i4), public, parameter  :: solv_max_iters = 124
   integer(i4), public, parameter  :: ntrials = 226
   
   integer, public :: boundary_exchange_algorithm = ALG_MPI2S_1D
   
   real (r4), parameter, public :: undefined = 1.0e35
   integer(i4), parameter, public :: undefined_nf_int = 100000

   real (r8), public :: &
      pi, pih, pi2            ! pi, pi/2 and 2pi

   !*** location of fields for staggered grids

   integer (i4), parameter, public ::   &
      field_loc_unknown  =  0, &
      field_loc_noupdate = -1, &
      field_loc_center   =  1, &
      field_loc_NEcorner =  2, &
      field_loc_Nface    =  3, &
      field_loc_Eface    =  4

   !*** field type attribute - necessary for handling
   !*** changes of direction across tripole boundary

   integer (i4), parameter, public ::   &
      field_type_unknown  =  0, &
      field_type_noupdate = -1, &
      field_type_scalar   =  1, &
      field_type_vector   =  2, &
      field_type_angle    =  3

   !  common formats for formatted output

   character (1), parameter, public :: &
      char_delim = ','

   character (9), parameter, public :: &
      delim_fmt = "(72('-'))",         &
     ndelim_fmt = "(72('='))"


   character (5), parameter, public :: &
      blank_fmt = "(' ')"

!  !PUBLIC DATA MEMBERS:

   character (char_len), public ::  &
      char_blank          ! empty character string

!***********************************************************************

 contains

!***********************************************************************

!>
!! This subroutine initializes constants that are best defined
!! at run time (e.g. pi).
!<
 subroutine init_constants

!-----------------------------------------------------------------------

   integer (i4) :: n

!-----------------------------------------------------------------------
!
!  more numbers and character constants
!
!-----------------------------------------------------------------------

   pi  = c4*atan(c1)
   pi2 = c2*pi
   pih = p5*pi


   do n=1,char_len
     char_blank(n:n) = ' '
   end do

!EOC
!-----------------------------------------------------------------------

 end subroutine init_constants

!***********************************************************************

 end module constants

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
