!###############################################################################
!   ASSIGNMENT 4 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS                   
!###############################################################################

program ising_model

    !Mersenne Twister RNG
    use mtmod
    
    integer, parameter :: L = 20, ns = L*L
    integer, parameter :: nr=1000000, nequ = 10000, ntau = 50000
    real, parameter    :: T0 = 1.0, T1 = 4.0, dT = 0.1, m0 = 0.5
    
    integer :: S(0:L+1,0:L+1) !This is our spin matrix
    real*8  :: M,E            !Magnetisation and energy
    real*8  :: mts(ntau)      !Time series of the magnetisation
    
    
    
    integer :: j,t,i
    integer :: tic,tau
    
    character(len=20) :: taufile,namefile
    
    !Initialise RNG with time
    call sgrnd(time())
    
    tic = time()              !Performance evaluation
    
    !Open file to store the autocorrelation time curve
    !write(taufile,"(A6,I2,A7)") "data/l",L,"nau.res"
    !open(12,file=taufile,action="write")
    
    
    !Start the cycle over the temperature
    do t = 15,15!nint((T1-T0)/dT)
        !It is easy to start from a configuration of all spins up
        S = 1
        M = ns          !Total Magnetisation
        E = -2*ns       !Total Energy
        
        !Surpass equilibration time
        do j = 1,nequ
            call advance_metropolis(S,L,T0+t*dT,M,E)
        enddo
        
        !Find Autocorrelation Time
        do j = 1,ntau
            call advance_metropolis(S,L,T0+t*dT,M,E)
            mts(j) = dabs(M)/ns
        enddo
        
        call autocorrelationtime(mts,ntau,tau)
        !write(12,"(F7.3,I5)") float(nint(T0*10)+t)/10,tau
        
        !Then start to sample
        write(namefile,"(A6,I2,A1,I2,A4)") "data/l",L,"t",nint(T0*10)+t,".res"
        open(11,file=namefile, action="write")
        do j=1,nr
            do i=1,2*tau
                call advance_metropolis(S,L,T0+t*dT,M,E)
            enddo
            write(11,"(2F14.7)") M/ns,E/ns
        enddo
        close(11)
        
    enddo
    
    !close(12)
    
end program

!###############################################################################

subroutine advance_metropolis(S,L,T,M,E)
    
    !Mersenne Twister RNG
    use mtmod

    integer :: L,x,y
    real*8  :: M,E,dE
    integer :: ns
    real    :: T,beta
    integer :: S(0:L+1,0:L+1)
    logical :: flipflag
    
    ns = L*L
    beta = 1.0/T
    
    !Start Metropolis sweep
    do i = 1,ns
        !Select a spin
        x = int(grnd()*L)+1
        y = int(grnd()*L)+1
        !Calculate dE. Use + because consider flipped spin in x,y
        dE = dfloat(2*S(x,y)*(S(x+1,y)+S(x-1,y)+S(x,y+1)+S(x,y-1)))
        if (dE<=0) then
            S(x,y) = -S(x,y)
            flipflag = .TRUE.
        else
            if (grnd()<exp(-beta*dE)) then
                S(x,y) = -S(x,y)
                flipflag = .TRUE.
            else
                flipflag = .FALSE.
            endif
        endif
        !If it was on the border, need to update the ghost cell
        if (flipflag .eqv. .TRUE.) then
            if (x==1) S(L+1,y) = -S(L+1,y)
            if (x==L) S(0,y) = -S(0,y)
            if (y==1) S(x,L+1) = -S(x,L+1)
            if (y==L) S(x,0) = -S(x,0)
            M = M + dfloat(2*S(x,y))
            E = E + dE
        endif
    enddo
    !End of a sweep of Metropolis
    
end subroutine

!###############################################################################

subroutine autocorrelationtime(obs,ntau,intau)
    
    integer :: ntau,intau
    real*8  :: obs(ntau)
    real*8  :: chi0, chi, a, b, c, tau
    integer :: i,j,cap
    
    cap = 1000 !We do not want to go past this
    
    chi0 = sum(obs**2)/ntau - (sum(obs)/ntau)**2
    tau = 0.d0
    do i = 1, cap
        a = 0.d0
        b = 0.d0
        c = 0.d0
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
        
        if (chi<1.d-3) then
            intau = int(tau) + 1
            return
        endif
    enddo
        
    intau = int(tau) + 1
    return
    
end subroutine
