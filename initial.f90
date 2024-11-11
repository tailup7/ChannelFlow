    program initial
    implicit none
    integer,parameter :: NX=64, NY=64, NZ=64
    character(3)::ct
    character(4)::file1
    integer :: i, j, k, istart
    double precision :: tstep
    double precision :: u(1:NZ, 1:NX, 1:NY), v(1:NZ, 1:NX, 0:NY), w(1:NZ, 1:NX, 1:NY) ,p(1:NZ, 1:NX, 1:NY)

    ct='000'
    file1='dns2'
    istart = 0
    tstep = 0.0000

    do j = 1, NY
        do i = 1, NX
          do k = 1, NZ
            U(k, i, j)=1
            W(k, i, j)=0.01*(65-I) 
            P(k, i, j)=100-I
          end do
        end do
    end do

    do j = 1, NY-1
        do  i = 1, NX
            do k = 1, NZ
              V(k, i, j)=0              
            end do
        end do 
    end do

    OPEN(10,file=FILE1//CT//'u.d',status='replace',form='formatted')            
      WRITE(10,1000) istart,tstep
      WRITE(10,2000) U
      CLOSE(10)
    OPEN(10,file=FILE1//CT//'v.d',status='replace',form='formatted')
      WRITE(10,1000) istart,tstep
      WRITE(10,2000) V
      CLOSE(10)
    OPEN(10,file=FILE1//CT//'w.d',status='replace',form='formatted')
      WRITE(10,1000) istart, tstep
      WRITE(10,2000) W
      CLOSE(10)
    OPEN(10,file=FILE1//CT//'p.d',status='replace',form='formatted')
      WRITE(10,1000) istart, tstep
      WRITE(10,3000) P
      CLOSE(10)

1000 FORMAT(I10,F12.6)
2000 FORMAT(8F10.6)
3000 FORMAT(8F10.5)

      write(*,'(a12)')'make initial'
    return
    end