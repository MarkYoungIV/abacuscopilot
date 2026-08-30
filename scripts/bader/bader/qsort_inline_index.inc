! FILE:qsort_inline_index.inc
! PURPOSE:
! Common internal procedures for sorting by index, for
! use with "qsort_inline.inc".

! set up initial index:
  subroutine init()
    integer :: i
    do i=1,array_size
      index(i)=i
    end do
  end subroutine init

! swap indices a,b
  subroutine swap(a,b)
    integer, intent(in) :: a,b
    integer :: hold
    hold=index(a)
    index(a)=index(b)
    index(b)=hold
  end subroutine swap

! circular shift-right by one:
  subroutine rshift(left,right)
    implicit none
    integer, intent(in) :: left, right
    integer :: hold, i
    hold=index(right)
    ! This sytnax is valid, but has poor optimization in GFortran:
    ! index(left+1:right)=index(left:right-1)
    do i=right,left+1,-1
      index(i)=index(i-1)
    end do
    index(left)=hold
  end subroutine rshift
