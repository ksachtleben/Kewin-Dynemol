! Program for computing Surface Hopping forces and transition probabilities from eHuckel Hamiltonian

module Surface_Hopping

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only  : driver , verbose , n_part , QMMM
    use Structure_Builder       , only  : Unit_Cell 
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: SH_Force , PES  , states_eh

    private

!    type( MM_atomic ) , allocatable :: d_vec(:,:,:)
    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]
    real*8  , parameter :: delta      = 1.d-8
    real*8  , parameter :: ua_2_kg = 1.660540199d-27
    real*8  , parameter :: ev_2_J  = 1.60218d-19
    logical , parameter :: T_ = .true. , F_ = .false.


    !module variables ...
    integer                     :: mm , newPES(2) , fermistate
    integer     , allocatable   :: BasisPointer(:) , DOS(:), PES(:)
    real*8      , allocatable   :: Kernel(:,:)  , d_vec(:,:,:,:)
    real*8      , allocatable   :: grad_S(:,:) , F_vec(:) , F_mtx(:,:,:)  
    logical     , allocatable   :: mask(:,:)

    real*8      , allocatable   :: Omega(:,:) , pastQR(:,:) , newQR(:,:) , g_switch(:,:)
    logical     , save          :: flip , done = F_


contains
!
!
!
!==================================================================
 subroutine SH_Force( system , basis , MO_bra , MO_ket ,  QM , dt ) 
!==================================================================
 use MD_read_m              , only  : atom
 use Semi_empirical_parms   , only  : atomx
 use parameters_m           , only  : electron_state , hole_state
 implicit none
 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 type(R_eigen)              , intent(in)    :: QM
 real*8                     , intent(in)    :: dt
 complex*16                 , intent(in)    :: MO_bra(:,:) , MO_ket(:,:)

! local variables ...
integer                  :: i , j , xyz , jL , L , indx , nn , iS
integer                  :: k , ik , DOS_atom_k , BasisPointer_k 
integer                  :: j1 , j2 , dima , dimb

! local arrays ...
integer , allocatable :: pairs(:)
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) , R1(:,:) , R2(:,:) , R3(:,:), VN(:,:) 
real*8                :: Force(3) , tmp_coord(3) , delta_b(3) 
real*8  , allocatable :: X_(:,:) , A1(:,:) , A2(:,:) , A3(:,:) , A4(:,:) , A5(:,:) , A6(:,:) , A7(:,:)

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 integer :: seed(5)

 seed = [10051965,27092004,2092002,22021967,-76571]

!================================================================================================
! some preprocessing ...
!================================================================================================
if( .NOT. allocated( PES ) ) then
allocate( PES( 2 ) )
PES(1) = electron_state
PES(2) = hole_state

call random_seed( put=seed(1:5) )
endif

!call init_random_seed()

mm = size(basis)
nn = n_part

fermistate = QM% Fermi_state

allocate( F_mtx     (system%atoms,system%atoms,3) , source=D_zero )
allocate( F_vec     (system%atoms)                , source=D_zero )
allocate( d_vec     (system%atoms,mm,mm,3)         , source=D_zero )
allocate( VN        (system%atoms,3)              , source= D_zero )
allocate( g_switch  (mm,mm)                        , source= D_zero )
If( .NOT. allocated(Kernel) ) then
    allocate( Kernel  (mm,mm) , source = D_zero )
end if

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis )
CALL preprocess    ( system )

CALL Huckel_stuff( basis , X_ )

! set all forces to zero beforehand ...
! using %Ehrenfest(:) to store the SH force ...
forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero

verbose     = .false.

!force on atom site ...
do k = 1 , system% atoms
If( system%QMMM(k) == "MM" .OR. system%flex(k) == F_ ) cycle

        DOS_atom_k     =  atomx( system% AtNo(k) )% DOS
        BasisPointer_k =  system% BasisPointer(k) 

        allocate( grad_S  (mm,DOS_Atom_K) )
        
        grad_S      = D_zero
        dima        = size( grad_S( : , 1 ) )
        dimb        = size( grad_S( 1 , : ) )

        ! temporary arrays ...
        allocate( A1(dima,dimb) , A2(dima,dima) , A3(dima,dimb) , A4(dima,dima) , A5(dima,dimb) , A6(dima,dimb) , A7(dima,dima) ) 
        allocate( R1(dimb,dima) , R2(dima,dima) , R3(dima,dimb) )


        ! save coordinate ...
        tmp_coord = system% coord(k,:)

        do xyz = 1 , 3

                VN(k,xyz) = atom(k)% vel(xyz) * mts_2_angs / sec_2_pico

                delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )
 
                system% coord (k,:) = tmp_coord + delta_b
                CALL Overlap_Matrix( system , basis , S_fwd , purpose = "Pulay" , site = K )

                system% coord (k,:) = tmp_coord - delta_b
                CALL Overlap_Matrix( system , basis , S_bck , purpose = "Pulay" , site = K )

                forall( j=1:DOS_Atom_K ) grad_S(:,j) = ( S_fwd( : , BasisPointer_K+j ) - S_bck( : , BasisPointer_K+j ) ) / (TWO*delta) 

                forall( j=1:DOS_Atom_K ) A1(:,j) = grad_S(:,j) * X_( : , BasisPointer_K + j )   ! N_{\alpha n }
                forall( j=1:mm )         A7(j,:) = QM%L(j,:)*QM%erg(j)                          ! Q'_{\varphi \alpha }
                forall( j=1:DOS_Atom_K ) A3(:,j) = QM%L(:, BasisPointer_K + j )                 ! Q_{\varphi n }
                forall( j=1:mm )         A2(j,:) = QM%L(j,:)                    ! \Phi_{\alpha}    
                forall( j=1:mm )         A4(j,:) = QM%L(j,:)*QM%erg(j)    ! Q'_{\phi \alpha}
                forall( j=1:DOS_Atom_K ) A5(:,j) = QM%L(: , BasisPointer_K+j )  ! Q_{\phi n}          

                    !==========================================================
                    call gemm( A2, A1 , R3 )   ! Qphi_{alpha}* N_{\alpha n} result vector size n
                    call gemm( R3 , A3 , R2 , transb="T" )
                    d_vec(k,:,:,xyz) = d_vec(k,:,:,xyz) - R2

                    !==========================================================
                    call gemm( QM%L , A1 , R3 )             ! Q_{\varphi \alpha} * N_{\alpha n} result matrix size \varphi n
                    call gemm( A5 , R3 , R2 , transb="T" )
                    d_vec(k,:,:,xyz) = d_vec(k,:,:,xyz) - R2

                    !==========================================================
                    call gemm( A4 , grad_S , R3 ) ! (<\alpha|\nabla_{xyz}n>)^{T} * Q'_{\phi \alpha} result matrix size phi n
                    call gemm( R3 , A3 , R2 , transb="T" )
                    d_vec(k,:,:,xyz) = d_vec(k,:,:,xyz) + R2

                    !==========================================================
                    call gemm( A7 , grad_S , R3 )           ! Q'_{\varphi \alpha } * <\alpha|\nabla_{xyz}n>
                    call gemm( A5 , R3 , R2 , transb="T" )
                    d_vec(k,:,:,xyz) = d_vec(k,:,:,xyz) + R2

                Force(xyz) = d_vec(k,PES(1),PES(1),xyz) - d_vec(k,PES(2),PES(2),xyz)

                forall( i=1:mm , is=1:mm , i .ne. is ) d_vec(k,i,is,xyz) = d_vec(k,i,is,xyz) / ( QM%erg(is) - QM%erg(i) ) 
                forall( i=1:mm , j=1:mm , i == j ) d_vec(k,i,j,xyz) = d_zero

                forall( i = 1 : mm ) g_switch(i,:) = g_switch(i,:) + VN(k,xyz) * d_vec(k,i,:,xyz)  
        end do 

        atom(k)% Ehrenfest = Force * eVAngs_2_Newton 

        ! recover original system ...
        system% coord (K,:) = tmp_coord

       deallocate(A1,A2,A3,A4,A5,A6,A7,R1,R2,R3 , grad_S )

end do 

call Hopping( system , basis , MO_bra , MO_ket , QM , dt )

deallocate( mask , X_ , F_vec , F_mtx , kernel , d_vec , VN , g_switch )

include 'formats.h'

end subroutine SH_Force
!
!
!
!=====================================
 subroutine Huckel_stuff( basis , Xi ) 
!=====================================
use Hamiltonians , only : X_ij
implicit none
type(STO_basis)  , intent(in) :: basis(:)
real*8           , allocatable , intent(out) :: Xi(:,:)

!local variables ...
integer :: i , j

allocate ( Xi(mm,mm) )

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

do j = 1 , mm
do i = j , mm

         Xi(i,j) = X_ij( i , j , basis )

         Xi(j,i) = Xi(i,j) 

end do
end do    

end subroutine Huckel_stuff
!
!
!
!
!============================
 subroutine Preprocess( sys ) 
!============================
use Semi_empirical_parms , only: atom
implicit none
type(structure) , intent(in) :: sys

!local variables ...
real*8                :: R_LK
integer               :: K , L
logical               :: flag1 , flag2 , flag3

If( .NOT. allocated(BasisPointer) ) allocate( BasisPointer(sys%atoms) , DOS(sys%atoms) )

Allocate( mask(sys%atoms,sys%atoms) , source = .false. )

do K = 1   , sys% atoms
   do L = K+1 , sys% atoms
   
       R_LK = sqrt(sum( (sys%coord(K,:)-sys%coord(L,:))**2 ) )
   
       flag1 = R_LK < cutoff_Angs  
        
       flag2 = sys% flex(K) .AND. sys% flex(L)
   
       flag3 = (sys% QMMM(L) == "QM") .AND. (sys% QMMM(K) == "QM")
   
       mask(L,K) = flag1 .AND. flag2 .AND. flag3

   end do
   BasisPointer(K) = sys% BasisPointer(K) 
   DOS(K)          = atom( sys% AtNo(K) )% DOS
end do    

end subroutine Preprocess
!
!
!================================================================
 subroutine Hopping( system , basis , MO_bra , MO_ket , QM , dt )
!================================================================
 use MD_read_m              , only  : atom
 use parameters_m           , only  : electron_state , hole_state
 implicit none

 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 type(R_eigen)              , intent(in)    :: QM
 complex*16                 , intent(in)    :: MO_bra(:,:)
 complex*16                 , intent(in)    :: MO_ket(:,:)
 real*8                     , intent(in)    :: dt

 real*8     , allocatable   :: rho_eh(:,:) , par_hop(:,:) , g_switche(:,:) , g_switchh(:,:) 
 real*8     , allocatable   :: A1(:,:) , A2(:,:) , A3(:,:) , A4(:), A5(:) , aux(:,:)
 integer    , allocatable   :: states_eh(:) 
 integer                    :: j , i , k ,l , ehs , h_en_c , lumo
 real*8                     :: Ergk , hop_to , hop_fr , rn , a_r , b_r , c_r , en_c , dv 
 real*8                     :: direction(3) , t(2) 
 integer                    :: eh_switch(2)
 logical                    :: p_t 


mm = size(basis)
allocate( rho_eh (mm,mm) )
allocate( aux (mm,mm) )
allocate( g_switche (mm,mm) )
allocate( g_switchh (mm,mm) )
allocate( par_hop(0:mm,2) , source=D_zero )

!temporary
allocate( A1  (system%atoms,3) , source= D_zero )
allocate( A2  (system%atoms,3) , source= D_zero )
allocate( A3  (system%atoms,3) , source= D_zero )
allocate( A4  (system%atoms)   , source= D_zero )
allocate( A5  (system%atoms)   , source= D_zero )

do k = 1 , mm
rho_eh(:,k) = real( MO_ket(:,1) * MO_bra(k,1) ) - real( MO_ket(:,2) * MO_bra(k,2) ) 
rho_eh(:,k) = rho_eh(:,k) / ( real( MO_ket(:,1) * MO_bra(:,1) ) - real( MO_ket(:,2) * MO_bra(:,2) ) )
enddo

g_switche = two * dt * rho_eh * g_switch
g_switche = max( d_zero , g_switche )
g_switchh = two * dt * rho_eh * g_switch
g_switchh = max( d_zero , -1.0d0*g_switchh )

call random_number( rn )

par_hop = 0.d0

forall( i = 1:mm ) par_hop(i,1) = sum( g_switche( PES(1) , 1:i ) ) 
forall( i = 1:mm ) par_hop(i,2) = sum( g_switchh( PES(2) , 1:i ) ) 

eh_switch = .false.
newPES    = PES

do j = 1 , 2
        states_eh = func_eh_state(j)

        check: do i = 1 , size( states_eh )
                k = states_eh( i )

                hop_fr  = par_hop( k-1 , j )
                hop_to  = par_hop( k , j )

                if( hop_fr < rn .and. rn < hop_to ) then
                        newPES(j)       = k
                        eh_switch(j)    = .true.
                        print*, j , k
                end if

                if( eh_switch(j) == .true.) exit check

         enddo check

enddo

switch: if( any( eh_switch == .true. ) ) then
        do k = 1 , system% atoms
        
                A2(k,:) = d_vec(k,PES(1),newPES(1),:) - d_vec(k,PES(2),newPES(2),:)                         ! d 
                A4(k)   = sqrt( sum( A2(k,:)**2 ) ) 
                       
                if( A4(k) .ne. 0.d0 ) direction(:) = A2(k,:) / A4(k)
                if( A4(k)  ==  0.d0 ) direction(:) = 0.d0
                        
                A1(k,:) = direction(:) * ( atom(k)%mass * imol )                 ! md
                A2(k,:) = direction(:)
                A3(k,:) = atom(k)% vel(:)                                         ! v  
                
                A4(k)   = dot_product( A1(k,:) , A2(k,:) )
                A5(k)   = dot_product( A1(k,:) , A3(k,:) )
        
        enddo

        a_r     = sum( A4 ) / 2.0d0
        c_r     = ( ( QM%erg(newPES(1)) - QM%erg(newPES(2)) ) - ( QM%erg(PES(1)) - QM%erg(PES(2)) ) )* ev_2_J
        b_r     = sum( A5 ) 

        en_c    = b_r*b_r - 4.0d0 * a_r * c_r
        ehs     = sum( abs( eh_switch ) )
        h_en_c  = heaviside( en_c )

        erg_c: select case( h_en_c )
        case( 0 ) ! dont have energy enough to hop
                 
                hop_eh: select case( ehs ) ! verify which hop occurs
                case( 1 ) ! only h or e
                        do k = 1, system% atoms
                        If( system%QMMM(k) == "QM" .AND. system%flex(k) == T_ ) then                                
                                atom(k)% vel(:) = atom(k)% vel(:) - 2.0d0 * ( dot_product(atom(k)% vel(:) , A2(k,:) ) ) * A2(k,:)
                                endif
                                enddo

                        exit switch

                case( 2 ) ! h and e

                        exit switch

                end select hop_eh
                        
        case( 1 ) ! energy conserv

                dv   = roots_square( a_r , b_r , c_r  )
                do k = 1 , system%atoms
                        If( system%QMMM(k) == "QM" .AND. system%flex(k) == T_ ) then
                                atom(k)% vel(:) = atom(k)% vel(:) - dv * A2(k,:) 
                                endif
                                end do
                PES = newPES
                exit switch

        end select erg_c
                
endif switch

electron_state = PES(1)
hole_state = PES(2)

deallocate( rho_eh , par_hop , A4 , A1 , A2 , A3 , A5 , g_switche , g_switchh )

end subroutine

!===========================================
 function func_eh_state(k) result(states_eh)
!===========================================
 use parameters_m           , only  : electron_state , hole_state
 integer , intent(in) :: k 

 integer               :: i , Neh ,j
 integer , allocatable :: temp(:) , states_eh(:)

 allocate( temp( mm ) , source= 0 )

 select case( k )
        case( 2 )
                do i = 1 , fermistate
                        temp(i) = i
                enddo
                temp( fermistate+1 ) = electron_hole
        case( 1 )
                j = 2 
                temp(1) = hole_state
                do i = fermistate+1 , mm
                        temp(j) = i
                        j = j + 1
                enddo
        end select

 Neh  = count( temp/=0 ) 
 allocate( states_eh(Neh) , source = temp )
 deallocate( temp )

end function
!
!
!
!===========================================
 function roots_square( a , b , c )
!===========================================
 real*8  ,       intent(in)      :: a , b , c


 real*8  ,       dimension(2)    :: min_value_abs
 real*8                          :: roots_square

 min_value_abs( 1 ) = -b/( two*a ) - sqrt( b*b - four*a*c )/( two*a )
 min_value_abs( 2 ) = -b/( two*a ) + sqrt( b*b - four*a*c )/( two*a )

 roots_square = minval( abs( min_value_abs ) )

end function
!
!===========================================
  function heaviside( x )
!===========================================
  real*8  ,       intent(in)      :: x

  heaviside = sign( 0.5d0 , x ) + 0.5d0

end function
!  
!
!==============================
 subroutine init_random_seed ()
!==============================
implicit none

!local variables ...
integer :: seed(5)

seed = [10051965,27092004,2092002,22021967,-76571]

call random_seed(put=seed(1:5))

end subroutine init_random_seed
!
!
end module Surface_Hopping
