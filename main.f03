!###############################################################################
!   ASSIGNMENT 4 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS                   
!###############################################################################

program ising_model

!The code is composed of several subroutines that do different kind of stuff
!The subroutine ADVANCE performs an entire sweep of the system using Metropolis   
    
    !Mersenne Twister RNG
    use mtmod
    
    integer, parameter :: L = 40, ns = L*L, nr=1000000
    real, parameter    :: T0 = 1.0, T1 = 4.0, dT = 0.1, m0 = 0.5
    integer :: S(0:L+1,0:L+1) 
    integer :: M,E
    integer :: i,j,t
    integer :: tic
    character(len=20) :: namefile

    call sgrnd(1)!(time())
    
    tic = time()
    
    do t = 0,nint((T1-T0)/dT)
    
        !First generate a random initial configuration
        do j = 1,L
            do i = 1,L
                if (grnd()<m0) then
                    S(i,j) = +1
                else
                    S(i,j) = -1
                endif
            enddo
        enddo
        !This part will make the rest of the code easier to write
        S(0,:) = S(L,:)
        S(L+1,:) = S(1,:)
        S(:,0) = S(:,L)
        S(:,L+1) = S(:,1)
    
        !Initialise variables
        M = sum(S(1:L,1:L))
        E = 0
        do j = 1,L
            do i = 1,L
                E = E - S(j,i)*(S(j,i+1)+S(j-1,i))
            enddo    
        enddo
        !End of Initialisation
    
        !Surpass equilibration time
        do j = 1,10000
            call advance_metropolis(S,L,T0+t*dT,M,E)
        enddo
        
        !Then start to sample
        write(namefile,"(A6,I2,A1,I2,A4)") "data/l",L,"t",nint(T0*10)+t,".res"
        open(11,file=namefile, status="unknown")
        do j=1,nr
            call advance_metropolis(S,L,T0+t*dT,M,E)
            write(11,*) j,float(M)/ns,float(E)/ns
        enddo
     enddo
     
     print *, "End",t,time()-tic

end

!###############################################################################

subroutine advance_metropolis(S,L,T,M,E)

    use mtmod

    integer :: L,M,E,dE,x,y
    integer :: ns
    real    :: T,beta
    integer :: S(0:L+1,0:L+1)
    logical :: flipflag

    ns = L*L
    beta = 1./T

    !Start Metropolis
    do i = 1,ns
        !Select a spin
        x = int(grnd()*L)+1
        y = int(grnd()*L)+1
        !Calculate dE. Use + because consider flipped spmin in x,y
        dE = 2*S(x,y)*(S(x+1,y)+S(x-1,y)+S(x,y+1)+S(x,y-1))
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
            M = M + 2*S(x,y)
            E = E + dE
        endif
    enddo
    !End of a sweep of Metropolis
    
end subroutine

!###############################################################################
            
    
