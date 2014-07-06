!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! Routuines to verify that the miniapp produced the correct answer.
!<
module check
    use kinds_mod, only: i4, r8
    use simple_blocks, only: nx_block,ny_block,get_block_parameter
    use reductions, only: global_sum
    use communicate, only: my_task, master_task
    use domain_size, only: max_blocks_tropic
    use simple_domain, only: blocks_tropic, nblocks_tropic

    implicit none

    public :: CheckAnswers

    contains

    !>
    !! This subroutine compares the expected solution read in from 
    !! the input file to that which is calculated by the minapp.
    !<
    subroutine CheckAnswers(algorithm, PRESSF, PRESS)
        character(*) :: algorithm
        integer(i4) :: nscan
        real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: PRESSF,PRESS

        real(r8) :: sum_diff,gdiff

        integer(i4) :: iblock,ib,ie,jb,je

        !----------------------------------
        ! check the accuracy of the solution
        !----------------------------------
        sum_diff = 0.0_r8
        do iblock=1,nblocks_tropic
            call get_block_parameter( &
                blocks_tropic(iblock),ib=ib,ie=ie,jb=jb,je=je)
            sum_diff = sum_diff + &
                SUM(PRESSF(ib:ie,jb:je,iblock)-PRESS(ib:ie,jb:je,iblock))
        enddo
        gdiff = global_sum(sum_diff)
        if(my_task == master_task) then 
          write(*,202) 'CheckAnswers: ',algorithm, ' SUM(calc-read): ',gdiff
!            print *,'CheckAnswers: ',algorithm, &
!                ' SUM(calc-read): ',gdiff
        endif
 202 format(a14,a,a16,f20.14)

    end subroutine CheckAnswers
end module check
