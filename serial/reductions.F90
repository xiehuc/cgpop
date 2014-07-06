!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module reductions

    ! !MODULE: reductions
    ! !DESCRIPTION:
    !  This module contains all the routines for performing global
    !  reductions like global sums, minvals, maxvals, etc.
    use kinds_mod, only : i4, r4, r8, log_kind 
    use simple_type, only: distrb
    use communicate, only :
    use constants, only : field_loc_NEcorner, field_loc_Nface, bignum, c0
    use blocks, only : get_block_parameter, nblocks_tot, nblocks_x, nblocks_y
    use domain_size, only : nx_global

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: global_sum,      &
        global_maxval,   &
        global_minval,   &
        init_global_reductions

    !-----------------------------------------------------------------------
    !
    !  generic interfaces for module procedures
    !
    !-----------------------------------------------------------------------
    interface global_sum
        module procedure  global_sum_vector_dbl,       &
            global_sum_vector_int,       &
            global_sum_scalar_dbl,       &
            global_sum_scalar_real,      &
            global_sum_scalar_int
    end interface

    interface global_maxval
        module procedure global_maxval_scalar_dbl
        module procedure global_maxval_scalar_real
        module procedure global_maxval_scalar_int
    end interface

    interface global_minval
        module procedure global_minval_scalar_dbl
        module procedure global_minval_scalar_real
        module procedure global_minval_scalar_int
    end interface

    !-----------------------------------------------------------------------
    !
    !  module variables
    !
    !-----------------------------------------------------------------------

    contains

    subroutine init_global_reductions

        ! !DESCRIPTION:
        !  Initializes necessary buffers for global reductions.

    end subroutine init_global_reductions


    function global_sum_scalar_dbl(local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the sum of scalar value across processors
        !
        !---------------------------------------------------------------
        real (r8), intent(in) :: &
            local_scalar                ! local scalar to be compared

        real (r8) :: &
            global_sum_scalar_dbl   ! sum of scalars across processors

        !---------------------------------------------------------------
        !
        !  no operation needed for serial execution
        !
        !---------------------------------------------------------------
        global_sum_scalar_dbl = local_scalar
    end function global_sum_scalar_dbl


    function global_sum_vector_dbl(local_vector, dist)
        !---------------------------------------------------------------
        !
        !  this function returns the sum of vector value across processors
        !
        !---------------------------------------------------------------
        type (distrb), intent(in) :: &
            dist                 ! distribution from which this is called

        real (r8), intent(in) :: &
            local_vector(:)   ! local scalar to be compared

        real (r8), dimension(size(local_vector,dim=1)) :: &
            global_sum_vector_dbl   ! sum of vectors across processors

        !---------------------------------------------------------------
        !
        !  no operation needed for serial execution
        !
        !---------------------------------------------------------------
        global_sum_vector_dbl = local_vector
    end function global_sum_vector_dbl


    function global_sum_vector_int(local_vector, dist)
        !---------------------------------------------------------------
        !
        !  this function returns the sum of vector value across processors
        !
        !---------------------------------------------------------------
        type (distrb), intent(in) :: &
            dist                 ! distribution from which this is called

        integer (i4), intent(in) :: &
            local_vector(:)                ! local scalar to be compared

        integer (i4), dimension(size(local_vector,dim=1)) :: &
            global_sum_vector_int   ! sum of vectors across processors

        !---------------------------------------------------------------
        !
        !  no operation needed for serial execution
        !
        !---------------------------------------------------------------
        global_sum_vector_int = local_vector
    end function global_sum_vector_int


    function global_sum_scalar_real(local_scalar, dist)
        !---------------------------------------------------------------
        !
        !  this function returns the sum of scalar value across processors
        !
        !---------------------------------------------------------------
        real (r4), intent(in) :: &
        local_scalar                ! local scalar to be compared

        type (distrb), intent(in) :: &
        dist                 ! distribution from which this is called

        real (r4) :: &
        global_sum_scalar_real   ! sum of scalars across processors

        !---------------------------------------------------------------
        !
        !  no operation needed for serial execution
        !
        !---------------------------------------------------------------
        global_sum_scalar_real = local_scalar
    end function global_sum_scalar_real


    function global_sum_scalar_int(local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the sum of scalar value across processors
        !
        !---------------------------------------------------------------
        integer (i4), intent(in) :: &
            local_scalar                ! local scalar to be compared

        integer (i4) ::      &
            global_sum_scalar_int   ! sum of scalars across processors

        !---------------------------------------------------------------
        !
        !  no operation needed for serial execution
        !
        !---------------------------------------------------------------
        global_sum_scalar_int = local_scalar
    end function global_sum_scalar_int


    function global_maxval_scalar_dbl (local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the maximum scalar value across processors
        !
        !---------------------------------------------------------------

        real (r8), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r8) :: &
        global_maxval_scalar_dbl   ! resulting global max

        !---------------------------------------------------------------
        !
        !  no operations required for serial execution - return input value
        !
        !---------------------------------------------------------------

        global_maxval_scalar_dbl = local_scalar
    end function global_maxval_scalar_dbl


    function global_maxval_scalar_real (local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the maximum scalar value across processors
        !
        !---------------------------------------------------------------

        real (r4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r4) :: &
        global_maxval_scalar_real   ! resulting global max

        !---------------------------------------------------------------
        !
        !  no operations required for serial execution - return input value
        !
        !---------------------------------------------------------------

        global_maxval_scalar_real = local_scalar
    end function global_maxval_scalar_real


    function global_maxval_scalar_int (local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the maximum scalar value across processors
        !
        !---------------------------------------------------------------

        integer (i4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        integer (i4) :: &
        global_maxval_scalar_int   ! resulting global max

        !---------------------------------------------------------------
        !
        !  no operations required for serial execution - return input value
        !
        !---------------------------------------------------------------

        global_maxval_scalar_int = local_scalar
    end function global_maxval_scalar_int


    function global_minval_scalar_dbl (local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !---------------------------------------------------------------

        real (r8), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r8) :: &
        global_minval_scalar_dbl   ! resulting global min

        !---------------------------------------------------------------
        !
        !  no operations required for serial execution - return input value
        !
        !---------------------------------------------------------------
        global_minval_scalar_dbl = local_scalar
    end function global_minval_scalar_dbl


    function global_minval_scalar_real (local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !-------------------------------------------------------------------
        real (r4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r4) :: &
        global_minval_scalar_real   ! resulting global min

        !-------------------------------------------------------------------
        !
        !  no operations required for serial execution - return input value
        !
        !-------------------------------------------------------------------
        global_minval_scalar_real = local_scalar
    end function global_minval_scalar_real


    function global_minval_scalar_int (local_scalar)
        !---------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !---------------------------------------------------------------
        integer (i4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        integer (i4) :: &
        global_minval_scalar_int   ! resulting global min

        !---------------------------------------------------------------
        !
        !  no operations required for serial execution - return input value
        !
        !---------------------------------------------------------------
        global_minval_scalar_int = local_scalar
    end function global_minval_scalar_int

end module reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
