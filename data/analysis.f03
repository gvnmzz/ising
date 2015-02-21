!###############################################################################
!   ASSIGNMENT 4 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS
!###############################################################################

program data_analysis

    use mtmod
    
    integer, parameter :: L = 20, nr=1000000, ntau = 50000, nboot = 1000 
    real, parameter    :: T0 = 1.0, T1 = 4.0, dT = 0.1
    integer, parameter :: ns = nint((T1-T0)/dT)
    real*8    :: mts(nr),ets(nr),prefact
    integer :: i,t,tau
    real*8    :: mmagn,mener,vmagn,vener,mchi,vchi,mspc,vspc
    integer :: autocorrelationtime
    character(len=10) :: namefile,datafile
    
    external autocorrelationtime
    
    !RNG Initialisation
    call sgrnd(time())
    
    write(datafile,"(A4,I2,A4)") "datl",L,".res"
    open(12,file=datafile,action="write")
    
    do t = 0,ns
        !print*,t
        !Open the file
        write(namefile,"(A1,I2,A1,I2,A4)") "l",L,"t",nint(T0*10)+t,".res"
        open(11,file=namefile, action="read")
        !Read the file
        do i=1,nr
            read(11,"(2F14.7)") mts(i),ets(i)
        enddo

        tau = autocorrelationtime(mts(1:ntau),ntau)
        call statistics(mts,nr,tau,mmagn,vmagn)
        call statistics(ets,nr,tau,mener,vener)
        prefact = 1.d0 / (dfloat(nint(T0*10)+t)/10) / (L*L)
        call bootstrap(mts,nr,tau,mchi,vchi,nboot,prefact)
        prefact = prefact / (dfloat(nint(T0*10)+t)/10)
        call bootstrap(ets,nr,tau,mspc,vspc,nboot,prefact)
        
        write(12,"(F14.7,I7,10(2x,E14.7))") float(nint(T0*10)+t)/10, tau, &
                    & mmagn,vmagn,mener,vener,mchi,vchi,mspc,vspc
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

!###############################################################################

subroutine statistics(obs,nr,tau,mean,var)
    
    integer :: nr,tau,i
    real*8  :: obs(nr),s,q,var,mean
    
    s = 0.d0
    q = 0.d0
    
    do i = 1,nr
        s = s + obs(i)
    enddo
    mean = s/nr
    do i = 1,nr
        q = (obs(i)-mean)**2
    enddo
    var = dfloat(1+2*tau)/(nr-1)*q
    
    return
end subroutine

!###############################################################################

subroutine bootstrap(obs,nr,tau,mean,var,nboot,prefact)
    
    use mtmod
       
    integer :: nr,tau,nboot,i,j
    real*8  :: obs(nr),samp(nr),var,mean,s,q
    real*8  :: conj(nboot),prefact
    
    !Create Sample
    do i = 1,nboot
        do j = 1,nr
            samp(j) = obs(int(grnd()*nr)+1)
        enddo         
        call statistics(samp,nr,tau,mean,conj(i))
    enddo
    conj = conj * prefact
    
    s = 0
    q = 0
    
    do i = 1,nboot
        s = s + conj(i)
    enddo
    mean = s/nboot
    do i = 1,nboot
        q = (conj(i)-mean)**2
    enddo
    var = q/nboot
    
    return
end subroutine

!###############################################################################
            
    
    
    
    
    
