!###############################################################################
!   ASSIGNMENT 4 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS
!###############################################################################

program data_analysis

    use mtmod
    
    integer, parameter :: L = 40, nr=1000000, ntau = 50000 
    real*8, parameter    :: T0 = 1.0, T1 = 4.0, dT = 0.1
    integer, parameter :: ns = nint((T1-T0)/dT)
    real*8    :: mts(nr),ets(nr),prefact
    integer :: i,t
    real*8    :: mmagn,mener,vmagn,vener,mchi,vchi,mspc,vspc,mamagn,vamagn
    character(len=10) :: namefile,datafile
    
    external autocorrelationtime
    
    !RNG Initialisation
    call sgrnd(time())
    
    write(datafile,"(A4,I2,A4)") "datl",L,".res"
    open(12,file=datafile,action="write")
    
    do t = 0,ns
        print*,t
        !Open the file
        write(namefile,"(A1,I2,A1,I2,A4)") "l",L,"t",nint(T0*10)+t,".res"
        open(11,file=namefile, action="read")
        !Read the file
        do i=1,nr
            read(11,"(2F14.7)") mts(i),ets(i)
        enddo
        
        call statistics(mts,nr,0,mamagn,vamagn)
        
        mts = dabs(mts)
        
        call statistics(mts,nr,0,mmagn,vmagn)
        call statistics(ets,nr,0,mener,vener)
        
        mchi = mmagn
        vchi = vmagn*nr
        prefact = dfloat(L*L) / (dfloat(nint(T0*10)+t)/10)
        call jacknife(mts,nr,mchi,vchi,prefact)

        mspc = mener
        vspc = vener*nr
        prefact = prefact / (dfloat(nint(T0*10)+t)/10)
        call jacknife(ets,nr,mspc,vspc,prefact)
        
        write(12,"(F14.7,10(2x,E14.7))") dfloat(nint(T0*10)+t)/10, &
                    & mmagn,vmagn,mener,vener,mchi,vchi,mspc,vspc,mamagn,vamagn
        !Close the file
        close(11)
    enddo
    close(12)
end


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
        q = q + (obs(i)-mean)**2
    enddo
    var = dfloat(1+2*tau)/(nr-1)*q/nr
    
    return
end subroutine

!###############################################################################
            
subroutine jacknife(obs,nr,mchi,vchi,prefact)
    
    integer :: nr,i
    real*8  :: mchi,vchi,obs(nr),conj(nr),prefact
    
    do i = 1,nr
        conj(i) = dfloat(nr-1)/(nr-2)*vchi-1.d0/(nr-2)*(obs(i)-mchi)**2
    enddo
    conj = conj * prefact
    
    call statistics(conj,nr,0,mchi,vchi)
    vchi = vchi*(nr-1)*nr
    
    return
end subroutine
    
    
    
    
