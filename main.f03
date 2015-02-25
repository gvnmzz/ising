!###############################################################################
!   ASSIGNMENT 4 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS                   
!###############################################################################

program ising_model

    !Mersenne Twister RNG
    use mtmod
    
    integer, parameter :: L = 20, ns = L*L
    integer, parameter :: nr=1000000, nequ = 10000, ntau = 50000
    real*8, parameter    :: T0 = 1.0, T1 = 4.0, dT = 0.1, m0 = 0.5
    
    integer :: S(0:L+1,0:L+1) !This is our spin matrix
    real*8  :: M,En,dEn,beta  !Magnetisation and energy
    real*8  :: mts(ntau)      !Time series of the magnetisation
    real*8  :: lfact
    
    
    integer :: j,t,i,ii,x,y
    integer :: tic,tau
    
    character(len=20) :: taufile,namefile
    
    !Initialise RNG with time
    call sgrnd(time())
    
    tic = time()              !Performance evaluation
    
    !Open file to store the autocorrelation time curve
    write(taufile,"(A6,I2,A7)") "data/l",L,"tau.res"
    open(12,file=taufile,action="write")
    
    
    !Start the cycle over the temperature
    do t = 0,nint((T1-T0)/dT)
    
        !It is easy to start from a configuration of all spins up
        S = 1
        M = ns          !Total Magnetisation
        En = -2*ns       !Total Energy
        beta = 1.d0/(T0+t*dT)
        lfact = L/(1.d0+epsilon(0.d0))
        
        !Surpass equilibration time
        do j = 1,nequ
            !This is an entire sweep of the Metropolis Algorithm
            do i = 1,ns
            !Select a spin
            x = int(grnd()*lfact)+1
            y = int(grnd()*lfact)+1
            !Calculate dE. Use + because consider flipped spin in x,y
            dEn = dfloat(2*S(x,y)*(S(x+1,y)+S(x-1,y)+S(x,y+1)+S(x,y-1)))
            if (grnd()<dexp(-beta*dEn)) then
                S(x,y) = -S(x,y)
                !If it was on the border, need to update the ghost cell
                if (x==1) S(L+1,y) = -S(L+1,y)
                if (x==L) S(0,y) = -S(0,y)
                if (y==1) S(x,L+1) = -S(x,L+1)
                if (y==L) S(x,0) = -S(x,0)
                M = M + dfloat(2*S(x,y))
                En = En + dEn
            endif
            enddo
        enddo
        
        !Find Autocorrelation Time
        do j = 1,ntau
            !This is an entire sweep of the Metropolis Algorithm
            do i = 1,ns
            !Select a spin
            x = int(grnd()*lfact)+1
            y = int(grnd()*lfact)+1
            !Calculate dE. Use + because consider flipped spin in x,y
            dEn = dfloat(2*S(x,y)*(S(x+1,y)+S(x-1,y)+S(x,y+1)+S(x,y-1)))
            if (grnd()<dexp(-beta*dEn)) then
                S(x,y) = -S(x,y)
                !If it was on the border, need to update the ghost cell
                if (x==1) S(L+1,y) = -S(L+1,y)
                if (x==L) S(0,y) = -S(0,y)
                if (y==1) S(x,L+1) = -S(x,L+1)
                if (y==L) S(x,0) = -S(x,0)
                M = M + dfloat(2*S(x,y))
                En = En + dEn
            endif
            enddo
            mts(j) = dabs(M)/ns
        enddo
        
        call autocorrelationtime(mts,ntau,tau)
        write(12,"(F7.3,I5)") float(nint(T0*10)+t)/10,tau
        cycle

        !Then start to sample
        write(namefile,"(A6,I2,A1,I2,A4)") "data/l",L,"t",nint(T0*10)+t,".res"
        open(11,file=namefile, action="write")
        do j=1,nr
          do ii=1,2*tau
            !This is an entire sweep of the Metropolis Algorithm
            do i = 1,ns
            !Select a spin
            x = int(grnd()*lfact)+1
            y = int(grnd()*lfact)+1
            !Calculate dE. Use + because consider flipped spin in x,y
            dEn = dfloat(2*S(x,y)*(S(x+1,y)+S(x-1,y)+S(x,y+1)+S(x,y-1)))
            if (grnd()<dexp(-beta*dEn)) then
                S(x,y) = -S(x,y)
                !If it was on the border, need to update the ghost cell
                if (x==1) S(L+1,y) = -S(L+1,y)
                if (x==L) S(0,y) = -S(0,y)
                if (y==1) S(x,L+1) = -S(x,L+1)
                if (y==L) S(x,0) = -S(x,0)
                M = M + dfloat(2*S(x,y))
                En = En + dEn
            endif
            enddo
          enddo
          write(11,"(4F14.7)") M/ns, En/ns
        enddo
        close(11)
        
    enddo
    
    close(12)
    
end program

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
