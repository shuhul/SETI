program makesigmadot
    implicit none
    integer,parameter :: dp = kind(-1.0d0)
    real(dp) :: r,s
    integer :: i,j,k,n,nr


    open(2,file='sigmadot.inp')
    open(1,file='rgrid.inp')
    read(1,*) n
    write(2,*) n
    do i=1,n
      read(1,*) r
      call sdot(r,s)
      write(2,'(ES12.3)') s
    enddo
    close(1)
    close(2)

end program makesigmadot


subroutine sdot(rad,rate)
    implicit none
    integer,parameter :: dp = kind(-1.0d0)
    real(dp) :: rad, rate
    integer :: func

!    func = 1 ! constant rate of sdot
!    func = 2 ! even depletion of sigma
!    func = 3 ! inside out dispersal
!    func = 4 ! outside in dispersal

    func=4

    ! numbers correspond to ~ 3 Myr disk lifetime
    select case (func)
    case(1)
         rate = 3.0e-15
    case(2)
         rate = 5.0e-16*(rad/50.0/1.5d13)**(-1)
    case(3)
         rate = 4.0e-16*(rad/50.0/1.5e13)**(-2)
    case(4)
         rate = 5.0e-14*(rad/50.0/1.5e13)**(2)
    end select 

    return
end subroutine sdot



    



