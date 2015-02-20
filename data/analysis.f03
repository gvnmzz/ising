!###############################################################################
!   ASSIGNMENT 4 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS
!###############################################################################

program data_analysis
    
    integer, parameter :: L = 40, nr=1000000, ntau = 50000 
    real, parameter    :: T0 = 1.0, T1 = 4.0, dT = 0.1
    integer, parameter :: ns = nint((T1-T0)/dT)
    real*8    :: mts(nr),ets(nr)
    integer :: i,t,tau
    integer :: autocorrelationtime
    character(len=10) :: namefile,taufile
    
    external autocorrelationtime
    write(taufile,"(A4,I2,A4)") "taul",L,".res"
    open(12,file=taufile,action="write")
    
    do t = 0,ns
        print*,t
        !Open the file
        write(namefile,"(A1,I2,A1,I2,A4)") "l",L,"t",nint(T0*10)+t,".res"
        open(11,file=namefile, action="read")
        !Read the file
        do i=1,nr
            read(11,*) mts(i),ets(i)
        enddo

        tau = autocorrelationtime(mts(1:ntau),ntau)
        
        write(12,*) float(nint(T0*10)+t)/10, tau
        !Close the file
        close(11)
    enddo
    close(12)
end

!###############################################################################

function autocorrelationtime(obs,ntau) result (intau)
    
    integer :: ntau,intau
    real*8    :: obs(ntau)
    real*8    :: chi0, chi, a, b, c, tau
    integer :: i,j
    
    chi0 = sum(obs**2)/ntau - (sum(obs)/ntau)**2
    tau = 0.
    print*,chi0
    do i = 1, ntau
        a = 0
        b = 0
        c = 0
        do j = 1, ntau - i
            a = a + obs(j)*obs(j+i)
            b = b + obs(j)
            c = c + obs(j+i)
        enddo
        a = a / (ntau - i)
        b = b / (ntau - i)
        c = c / (ntau - i)
        
        chi = (a - b * c)/chi0
        tau = tau + chi
        
        if (chi<1.e-3) then
            intau = int(tau) + 1
            return
        endif
    enddo
        
    intau = int(tau) + 1
    return
    
end function
