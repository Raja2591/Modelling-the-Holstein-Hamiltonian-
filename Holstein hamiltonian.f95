!******************************************************************************
!This module contains all the variables that the main program will need
!******************************************************************************


module common_variables 
implicit none
 integer conf_max, conf_process
        real*8 dwidth, ycorrlen, xcorrlen, time_start, time_finish,    &
               delta_da
        real*8, allocatable :: offset(:,:)
          !lattice numbers
    integer :: lattice_x = 1
    integer :: lattice_y = 1
    integer :: lattice_z = 1
    
    !mpi
    integer :: dim2 = 10
  
    !vibration
    real*8 :: hw = 1400.d0
    real*8 :: s  = 1.d0
    real*8 :: s_cat = 0.d0
    real*8 :: s_ani = 0.d0
    integer:: vibmax = 0
    
    !coupling
    real*8 :: jo_x = 0.d0
    real*8 :: jo_y = 0.d0
    real*8 :: jo_z = 0.d0
    real*8 :: ct_u = 0.d0
    real*8 :: ct_v = 0.d0
    real*8 :: te_x = 0.d0
    real*8 :: te_y = 0.d0
    real*8 :: te_z = 0.d0
    real*8 :: th_x = 0.d0
    real*8 :: th_y = 0.d0
    real*8 :: th_z = 0.d0

    !doping
    real*8 :: dintra = 0.d0
    real*8 :: dinter = 0.d0
    real*8 :: r = 0.d0
 
    !coupling
    logical :: extended_cpl = .false.
    real*8, allocatable :: jo_ex(:,:,:)
    character*256 extended_cpl_file

    !task title
    character*256 task_title

    !hamiltonian counters
    integer :: kount = 0
    integer :: kount_1p = 0
    integer :: kount_2p = 0
    integer :: kount_3p = 0
    integer :: kount_ct = 0
    integer :: kount_ct2 = 0
    integer :: kount_lattice = 0

    !indexes
   integer, allocatable :: nx_lattice(:,:,:)
    integer, allocatable :: nx_1p(:,:)
    integer, allocatable :: nx_2p(:,:,:,:)
    real *8, allocatable :: disorder(:,:,:)
    real *8, allocatable :: Force(:,:,:) 
    real, allocatable :: pdisorder(:,:,:,:,:)  
    real *8, allocatable :: te (:,:,:)
    real *8, allocatable :: D(:,:,:) 
    
    !the hamiltonian
    real*8, allocatable :: h(:,:)
    real*8, allocatable :: eval(:)
    
    !franck condon factors
    real*8, allocatable :: fc_gf(:,:)
    real*8, allocatable :: fc_gc(:,:)
    real*8, allocatable :: fc_ga(:,:)
    real*8, allocatable :: fc_af(:,:)
    real*8, allocatable :: fc_cf(:,:)

    !constants
    real*8, parameter :: pi = 4.d0*datan(1.d0)
    real*8, parameter :: ev = 8065.d0
    real*8, parameter :: hc = 1.23984193d3 * ev !(plancks constant times the speed of light in nm*hw )
    real*8, parameter :: kb = 0.6956925d0        !units of cm-1 k
    real*8, parameter :: angstrom = 0.0000000001 !units of m
    real*8, parameter :: beta = 2.35
    !empty
    integer, parameter :: empty = -1

    !multiparticle states
    logical :: one_state    =.true.
    logical :: two_state    =.false.
    logical :: LL    =.true.
    logical :: SS    =.false.
    logical :: LS = .false.
    logical :: offdiagonal = .true.
    logical :: diagonal = .false.
    logical :: dopant = .false.
    logical :: three_state  =.false.
    logical :: ct_state     =.false.
    logical :: ct2_state    =.false.
    
    !periodic
    logical periodic
    
    !oscillator strength
    real*8, allocatable :: osc_x(:), osc_y(:), osc_z(:)
    real*8, allocatable :: plosc_x(:,:), plosc_y(:,:), plosc_z(:,:)
    real*8, allocatable :: ux(:), uy(:), uz(:)
    real*8 :: abs_lw = 100.d0
    real*8 :: mon_tr_e = 0.d0
    integer, parameter  :: spec_step = 2600        !for 10 cm-1 resolution

    !number of threads
    integer :: numthreads = 1
    
    !dielectric
    real*8 :: dielectric=1.d0
    
    !ct states trunctation
    logical :: ct_truncate = .false.
    logical :: cty_truncate = .false.
    integer :: ct_range = 1000
    integer :: cty_range = 1000
    logical :: two_truncate = .false.
    integer :: two_range = 1000
    
    logical :: lorentzian = .false.
    logical :: linewidthprogress = .false.
    real*8 :: delta_sigma = 0.d0
    real*8, allocatable :: vibsum_evec(:), vibsum_state(:)
    logical :: abs_freq_dep = .false.
    
    
    integer :: kount2
    integer :: iu=1
    character(1) :: rrange='A'
    end module


! Module ends
!****************************************************************************************************
!__________________________________________________________________________________________________________________________________________________!

       
!****************************************************************************************************
                                !Main Program starts here!
!****************************************************************************************************
program cms
use common_variables
implicit none
include 'mpif.h'
integer ista,iend,ierr,iproc,nproc
integer config, complete

        integer :: bin1 = 0
        integer :: mbin1=0
        integer :: bin2 = 0
        integer :: mbin2=0
        integer :: bin3 = 0
        integer :: mbin3=0
        integer :: bin4 = 0
        integer :: mbin4=0
        integer :: bin5 = 0
        integer :: mbin5=0
        integer :: bin6 = 0
        integer :: mbin6=0
        integer :: bin7 = 0
        integer :: mbin7=0
        integer :: bin8 = 0
        integer :: mbin8=0
        integer :: bin9 = 0
        integer :: mbin9=0
        
    real*8  start, finish
    real*8 :: ab_x(spec_step) = 0.D0
    real*8 :: mab_x(spec_step) = 0.D0
    real*8 :: ab_y(spec_step) = 0.D0
    real*8 :: mab_y(spec_step) = 0.D0
    real*8 :: sumoscx = 0.d0, sumoscy = 0.d0
    real*8 :: msumoscx = 0.d0, msumoscy = 0.d0
    real*8, allocatable :: cohfxn(:,:,:), cohtmp(:,:,:)
    real*8, allocatable :: mcohfxn(:,:,:)    
    character*128 sstring
    
  
!****************************************************************************************************
  !We need to first set the parameters for the hamiltonian
    call set_para()
    
!****************************************************************************************************
  !We need to index the lattice and the basis state, that is the one particle and two particle states.
  
    call index_lattice()
    if ( one_state )   call index_1p()
    if ( two_state )   call index_2p()    
!****************************************************************************************************      


 !**************************output the total no. of one particle and two particle basis sets for our check*************!
    print*, '***************************'
    print*, 'kount   :', kount
    print*, 'kount 1p:', kount_1p
    print*, 'kount 2p:', kount_2p
    print*, '***************************'
!****************************************************************************************************      
  
    
!****************************************************************************************************
   !We need to set the Frank Condon factors for the vibrational overlap factors
    call set_fctable()
!****************************************************************************************************        
    
    !output the parameters (local)
    call para_out_cms()
    
    !initialize the hamiltonian and allocate (exciton_common_local)
      call allocate_hev()

  !*******************************************************************************************************************************!
    !We build the hamiltonian here; the one particle hamiltonian, two particle hamiltonian and the one particle -  two particle hamiltonian
       
        if ( one_state ) call hamiltonian_1p(config)
        if ( two_state ) call hamiltonian_2p(config)
        if ( one_state .and. two_state ) call hamiltonian_1p2p(config)
 !*******************************************************************************************************************************!        
        

 !*******************************************************************************************************************************!
                        !We diagonalize the hamiltonian and find the eigen values and the eigen vectors
  
       call cpu_time(start)
       call dsyevr_diagonalize( h, kount, eval, kount2, rrange, iu ) 
       call cpu_time(finish)
        print*, 'Diagonalization Time::', finish-start
 !*******************************************************************************************************************************!        

end program


!****************************************************************************************************
                               !main program ends
!****************************************************************************************************


!****************************************************************************************************
                  !Subroutines used for the program
!****************************************************************************************************


!****************************************************************************************************
       !Subroutine set_para()  sets the parameters that are needed for the hamiltonian. 
                                   
!****************************************************************************************************

subroutine set_para()
     use common_variables
    implicit none

    logical         exists
    character*100   buffer, label, fname
    integer         fno, ios, line, pos, errstat, num_threads
    parameter       ( fno = 90 )

    !Set the default parameters
    print*, 'Setting the default parameters.'
    lattice_x       = 4
    lattice_y       = 4
    vibmax          = 4
    hw              = 1400.d0
    s               = 1.d0
    jo_x            = -0.20d0        !These are really te or th, but using exciton subroutines
    jo_y            = -0.15d0       !I have to use the jo variables
    delta_da        = 0.d0        !donor - acceptor energy difference. Set to 0 for homopolymer
    task_title      = 'test'
    abs_lw          = 350.d0
    lorentzian      = .false.
    rrange          = 'A'           !By default, find all of the eigenstates
    iu              = 1           !disorder conf parameters
    conf_max        = 5000 
    dwidth          = 0.d0
    xcorrlen        = 0.d0
    ycorrlen        = 0.d0
    one_state       = .true.
    two_state       = .true.
    dintra = 0.4
    dinter = 0.4
     r = 1
    LL = .false.
    LS = .false.
    SS = .true.
     dopant = .false.
    diagonal = .true.
   offdiagonal = .false.
    two_truncate    = .false.
    two_range       = 3

      

    !Read the name of the input file
    call get_command_argument(1, fname, status=errstat )
    if ( errstat .ne. 0 ) goto 1010
    inquire( file = trim(fname), exist = exists)
    if ( .not. exists ) then
        print*, 'Input file not found...aborting'
        stop
    end if
    
    !Open the file to read in
    open( unit = fno, file = fname, status = 'old', action = 'read' )
    !-----------------------------------!
    !   This reads the control file
    !   The ios changes if end of record
    !   or end of file is reached
    !-----------------------------------!
    ios = 0
    line = 0
    print*, 'Reading the input file...'
    do while ( ios == 0 )
        read( fno, '(a)', iostat=ios ) buffer
        if( ios == 0 ) then
            line = line + 1

            !Find first instance of whitespace
            pos = scan( buffer, ' ' )
            label = buffer( 1:pos )
            buffer = buffer( pos + 1:)

            select case ( label )
            case('lattice_x')
                read( buffer, *, iostat=ios ) lattice_x
                print*, '    Setting lattice_x to: ', lattice_x
            case('lattice_y')
                read( buffer, *, iostat=ios ) lattice_y
                print*, '    Setting lattice_y to: ', lattice_y
            case('vibmax')
                read( buffer, *, iostat=ios ) vibmax
                print*, '    Setting vibmax to: ', vibmax
            case('hw')
                read( buffer, *, iostat=ios ) hw
                print*, '    Setting vibrational energy to (cm-1): ', hw
            case('s')
                read( buffer, *, iostat=ios ) s
                print*, '    Setting s to : ', s
            case('te_x')
                read( buffer, *, iostat=ios) jo_x
                print*, '    Setting te_x to (eV): ', jo_x
            case('th_x')
                read( buffer, *, iostat=ios) jo_x
                print*, '    Setting th_x to (eV): ', jo_x    
            case('delta_da')
                read( buffer, *, iostat=ios) delta_da
                print*, '    Setting delta_da to (cm-1): ', delta_da           
            case('te_y')
                read( buffer, *, iostat=ios) jo_y
                print*, '    Setting te_y to (eV): ', jo_y
            case('th_y')
                read( buffer, *, iostat=ios) jo_y
                print*, '    Setting th_y to (eV): ', jo_y
            case('task_title')
                read( buffer, *, iostat=ios) task_title
                print*, '    Setting task_title to: ', trim(task_title)
            case('abs_lw')
                read( buffer, *, iostat=ios) abs_lw
                print*, '    Setting the linewidth to (cm-1): ', abs_lw
            case('two_range')
                read( buffer, *, iostat=ios) two_range
                print*, '    Setting the vibrational radius to: ', two_range
            case('two_truncate')
                read( buffer, *, iostat=ios) two_truncate
                if ( two_truncate ) print*, '    Vibrational radius is on'
                if ( .not. two_truncate ) print*, '    Vibrational radius is off'
            case('lorentzian')
                read( buffer, *, iostat=ios) lorentzian
                if ( lorentzian ) print*, '    Setting the lineshape to:', &
                                          ' lorentzian'
                if ( .not.lorentzian ) print*, '    Setting the lineshape to:', &
                                          ' gaussian'
            case('one_state')
                read( buffer, *, iostat=ios) one_state
                if ( one_state ) print*, '    One particle states are turned on.'
                if ( .not.one_state ) print*, '    One particle states are turned off.'
            case('two_state')
                read( buffer, *, iostat=ios) two_state
                if ( two_state ) print*, '    Two particle states are turned on.'
                if ( .not.two_state ) print*, '    Two particle states are turned off.'
            case('rrange')
                read( buffer, *, iostat=ios) rrange
                if ( rrange .eq. 'A' ) print*, '    dsyevr will find all  eigenvectors'
                if ( rrange .eq. 'I' ) print*, '    dsyevr will find a range of eigenvectors'
            case('iu')
                read( buffer, *, iostat=ios) iu
                print*, '    desyevr will look for the lowest', iu, ' eigenvectors'           
            case('conf_max')
                read( buffer, *, iostat=ios) conf_max
                print*, '    Setting conf_max to: ', conf_max
            case('dwidth')
                read( buffer, *, iostat=ios) dwidth
                print*, '    Setting disorder width to (cm-1): ', dwidth
            case('xcorrlen')
                read( buffer, *, iostat=ios) xcorrlen
                print*, '    Setting xcorrlen length to: ', xcorrlen            
            case('ycorrlen')
                read( buffer, *, iostat=ios) ycorrlen
                print*, '    Setting ycorrlen length to: ', ycorrlen
            case default
                print*, '    invalid label at line, ', line
            end select
        end if

    end do

    close(fno)

1010 continue
    print*, 'Calculating derived parameters.'
    !Normalize parameters to hamiltonian units of vib quanta
    jo_x = jo_x * eV / hw
    jo_y = jo_y * eV / hw
    abs_lw = abs_lw / hw
    dwidth = dwidth / hw
    delta_da = delta_da / hw
    if ( vibmax == 0 ) s = 0.d0
    
!    call omp_set_num_threads(1)

end subroutine



!***************************************************************************************!
!    Write the parameters to a file. This is optional. Its upto you.
!***************************************************************************************!


subroutine para_out_cms()
    use common_variables
        implicit none
    
    character*64 fname
    integer fno
    parameter( fno = 16 )
    character*8 ddd
    character*10 ttt
    
    !create the file
    fname = trim(task_title)//'_para.csv'
    open( unit = fno, file = fname, action='write' )
    call date_and_time( ddd, ttt )
    
    write( fno, * ) 'parameter file for, cms_local_disorder.f95'
    write( fno, * ) 'task title, ', trim(task_title)
    write( fno, * ) 'date and time, ', ddd, ' ', ttt
    write( fno, * ) 'x dimension, ', lattice_x
    write( fno, * ) 'y dimension, ', lattice_y
    write( fno, * ) 'vibmax, ', vibmax
    write( fno, * ) 'vib energy (cm-1), ', hw 
    write( fno, * ) 'huang-rhys, ', s
    write( fno, * ) 'th x (eV), ', jo_x*hw/eV
    write( fno, * ) 'th y (eV), ', jo_y*hw/eV    
    write( fno, * ) 'abs lw (cm-1), ', abs_lw * hw
    if ( lorentzian ) write( fno, * ) 'lineshape, lorentzian'
    if ( .not. lorentzian ) write( fno, * ) 'lineshape, gaussian'    
    write( fno, * ) 'conformations, ', conf_max
    write( fno, * ) 'disorder width (cm-1), ', dwidth*hw
    write( fno, * ) 'x correlation length, ', xcorrlen
    write( fno, * ) 'y correlation length, ', ycorrlen
    write( fno, * ) 'one state, ', one_state
    write( fno, * ) 'two state, ', two_state
    write( fno, * ) 'two truncate, ', two_truncate
    write( fno, * ) 'two range, ', two_range   
    
    close( fno )
    
end subroutine



 
 
!**********************************************!
!                        Index the lattice	                                     !
!**********************************************!
subroutine index_lattice()
    use common_variables
    implicit none 

    integer lx, ly, lz
    
    if ( .not. allocated( nx_lattice ) ) then 
        allocate( nx_lattice( lattice_x, lattice_y, lattice_z ) )
    end if
    nx_lattice = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
        kount_lattice = kount_lattice + 1
        nx_lattice( lx, ly, lz ) = kount_lattice
    end do
    end do
    end do
        
end subroutine

!**********************************************!
!                        Index the 1p states                                   !
!**********************************************!
subroutine index_1p()
    use common_variables
    implicit none 
    
    integer lx, ly, lz, lxyz, vib

    !allocate the index
    if ( .not. allocated( nx_1p ) ) then 
            allocate( nx_1p ( kount_lattice, 0:vibmax ) )
    end if
    nx_1p = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
    do vib = 0, vibmax
        kount = kount + 1
        kount_1p = kount_1p + 1
        lxyz = nx_lattice( lx, ly, lz )
        nx_1p( lxyz, vib ) = kount
        write (9,*) 'nx_1p(', lxyz, vib,' )', kount
!                print*, lxyz, vib, kount
    end do
    end do
    end do
    end do
        
end subroutine
!**********************************************!
!                        Index the 2 p states                                   !
!**********************************************!

subroutine index_2p()
    use common_variables
    implicit none 
    
    integer lx, ly, lz, lxyz, vib
    integer lxv, lyv, lzv, lxyzv, vibv

    !allocate the index
    if ( .not. allocated( nx_2p ) ) then 
            allocate( nx_2p ( kount_lattice, 0:vibmax, kount_lattice, 1:vibmax ) )
    end if
    nx_2p = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
    do vib = 0, vibmax
        do lxv = 1, lattice_x
        do lyv = 1, lattice_y
        do lzv = 1, lattice_z
!                do vib = 0, vibmax
        do vibv = 1, vibmax
            lxyz = nx_lattice( lx, ly, lz )
            lxyzv = nx_lattice( lxv, lyv, lzv ) 
            if ( lxyz == lxyzv ) cycle                !not a 2 p state
            if ( vibv + vib > vibmax ) cycle	!truncate vibs at vibmax
            if ( two_truncate .and. two_range*1.d0 < dsqrt( 1.d0*((lx-lxv)**2+(ly-lyv)**2 + (lzv-lzv)**2)) ) cycle	!truncate at two_range
            kount = kount + 1
            kount_2p = kount_2p + 1
            nx_2p( lxyz, vib, lxyzv, vibv ) = kount
            write (9,*) 'nx_2p(', lxyz, vib, lxyzv, vibv,' )', kount
            !write (7,*) 'nx_2p(',lxyz,vib,lxyzv,vibv,')', kount
!                        print*, lxyz, vib,lxyzv, vibv, kount
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
        
end subroutine


!**********************************************!
!     create the one particle Hamiltonian                                 !
!**********************************************!
subroutine hamiltonian_1p(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    real*8 get_j

    if ( .not. allocated( vibsum_state ) ) allocate( vibsum_state(kount) )
    vibsum_state = 0.d0

    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
      
        hx1 = nx_1p( lxyz1, vib1 )
       
        if ( hx1 == empty ) cycle
            !count number of vibrations for this state
            vibsum_state( hx1 ) = vib1*1.d0
        !diagonal
        h( hx1, hx1 ) = vib1*1.d0
        write(7,*) 'h(',hx1,hx1,')', h(hx1,hx1)
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
        do vib2 = 0, vibmax
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            hx2 = nx_1p( lxyz2, vib2 )
            if ( hx2 == empty ) cycle
            if ( hx2 == hx1 ) cycle

            !off diagonal
            h( hx1, hx2 ) = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                            fc_gf( 0, vib1 ) * fc_gf( 0, vib2 )
            !it is hermitian                
            h( hx2, hx1 ) = h( hx1, hx2 )
                
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do

end subroutine

!***********************************************************!
!    create the one particle and two particle hamiltonian                                 
!***********************************************************!



subroutine hamiltonian_1p2p(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    integer lx2v, ly2v, lz2v, lxyz2v, vib2v
    real*8 get_j

    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
        hx1 = nx_1p( lxyz1, vib1 )
        if ( hx1 == empty ) cycle
                        
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
        do vib2 = 0, vibmax
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            
            !these are the only ones that are non-zero
            lx2v = lx1 
            ly2v = ly1
            lz2v = lz1
            lxyz2v = nx_lattice( lx2v, ly2v, lz2v )
            do vib2v = 1, vibmax

                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                if ( hx2 == empty ) cycle

                h( hx1, hx2 ) = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                                fc_gf( vib2v, vib1 ) * fc_gf( 0, vib2 )
                !it is hermitian                
                h( hx2, hx1 ) = h( hx1, hx2 )

            end do
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
end subroutine

!**********************************************!
!     create the two particle  hamiltonian     !
!**********************************************!

subroutine hamiltonian_2p(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx1v, ly1v, lz1v, lxyz1v, vib1v
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    integer lxyz2v, vib2v
    real*8 get_j
    
    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
    
        do lx1v = 1, lattice_x
        do ly1v = 1, lattice_y
        do lz1v = 1, lattice_z
        do vib1v = 1, vibmax
            lxyz1v = nx_lattice( lx1v, ly1v, lz1v )
    
            
            hx1 = nx_2p( lxyz1, vib1, lxyz1v, vib1v )
            
           ! write (7,*) 'nx_2p(',lxyz1,vib1,lxyz1v,vib1v,')', nx_2p( lxyz1, vib1, lxyz1v, vib1v )
                !count number of vibrations for this state
                vibsum_state( hx1 ) = (vib1 +vib1v)*1.d0
            if ( hx1 == empty ) cycle

            !diagonal
            h( hx1, hx1 ) = 1.d0*( vib1 + vib1v )                         
     write (7,*) 'h(', hx1, hx1, ')' , h(hx1,hx1)
            do lx2 = 1, lattice_x
            do ly2 = 1, lattice_y
            do lz2 = 1, lattice_z
                lxyz2 = nx_lattice( lx2, ly2, lz2 )
            do vib2 = 0, vibmax
            
                !linker type
                lxyz2v = lxyz1v
                vib2v = vib1v
                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                if ( hx2 == empty .or. hx2 == hx1 ) then
                    continue
                else
                    h( hx1, hx2 ) = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                    fc_gf( 0, vib1 ) * fc_gf( 0, vib2 )
                    h( hx2, hx1 ) = h( hx1, hx2 )
                end if
        
                !exchange type
                if ( lxyz1v == lxyz2 ) then
                    lxyz2v = lxyz1
                    do vib2v = 1, vibmax
                        hx2 = nx_2p(  lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty .or. hx2 == hx1 ) cycle
                        
                        h( hx1, hx2 ) = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                        fc_gf( vib2v, vib1 ) * fc_gf( vib1v, vib2 )
                        h( hx2, hx1 ) = h( hx1, hx2 )
                    end do
                end if
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
                
end subroutine

!**********************************************!
!  get the coupling function                   !
!**********************************************!
pure real*8 function get_j( lx1, lx2, ly1, ly2, lz1, lz2,config )
    use common_variables
    implicit none

    integer, intent(in) :: lx1, lx2, ly1, ly2, lz1, lz2,config
    integer dx, dy, dz
    
    dx = abs( lx1 - lx2 )
    dy = abs( ly1 - ly2 )
    dz = abs( lz1 - lz2 )
    if (periodic) then
        dx = min( dx, lattice_x - dx )
        dy = min( dy, lattice_y - dy )
        dz = min( dz, lattice_z - dz )
    end if
    
    !initialize
    get_j = 0.d0

    if ( .not. extended_cpl ) then
        if ( dy == 0 .and. dz == 0 ) then
                if (dx == 1 ) get_j = jo_x
                return
        end if
        if ( dx == 0 .and. dz == 0 ) then
                 if ( dy == 1 ) get_j = jo_y 
            return
          ! if ( dy == 1 .and. offdiagonal ) get_j = jo_y + ((jo_y*exp(-beta*pdisorder(config,lx1,ly1,lx2,ly2)))-jo_y)
           !      return
        end if
        if ( dx == 0 .and. dy == 0 ) then
                if ( dz == 1 ) get_j = jo_z
                return
        end if
    else if ( extended_cpl ) then
        get_j = jo_ex( dx, dy, dz )
    end if
    
    return

end function

!********************************************************!
!  This subroutine sets the Frank Condon Factors         !
!********************************************************!


subroutine set_fctable()
        use common_variables
        implicit none

	integer vib_g, vib_exc, vib_ion, vib_n
	! vib_g is for ground state, vib_exc for excited, vib_ion is for charged mol, and vib_n is for neutral mol.
	real*8 fc, s_diff

        !allocate matrices
        if ( .not. allocated( fc_gf ) ) then 
                allocate( fc_gf ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_gc ) ) then 
                allocate( fc_gc ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_ga ) ) then 
                allocate( fc_ga ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_cf ) ) then 
                allocate( fc_cf ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_af ) ) then 
                allocate( fc_af ( 0:vibmax, 0:vibmax ) )
        end if

	!----- generate fc table ----!
	do vib_g =0,vibmax
	do vib_exc =0,vibmax
		call fcfac( vib_g,vib_exc,s,fc)
		fc_gf(vib_g, vib_exc) = fc
	enddo
	enddo

        !ground to cation
        do vib_n  =0,vibmax
	do vib_ion=0,vibmax
		call fcfac(vib_n, vib_ion,s_cat,fc)
		fc_gc(vib_n, vib_ion) = fc
	enddo
        enddo
        
        !ground to anion
        do vib_n  =0,vibmax
	do vib_ion=0,vibmax
		call fcfac(vib_n, vib_ion,s_ani,fc)
		fc_ga(vib_n, vib_ion) = fc
	enddo
        enddo

	!---- generate fc table for ct ----!
	! the procedure is designed in a way that the first potential well is inside
	! the second. if this order is reversed, it may give wrong signe.
	! for example, if s_frenkel > se, then the order should be se then s_frenkel.
	! however, if s_frenkel < se, then the order becomes s_frenkel and se.
	! 
	! 1/27/2009
	! up to this point, s(1) = 1 and nuclear displacement factor of
	! ct states (lamda) are set such that the sum of the lamda is s(1), which is 1
	! since lamda^2 is s, root sct_h + root sct_e = root s(1)=1.
	
	!cation to frenkel	
	do vib_ion=0,vibmax
	do vib_n  =0,vibmax
		s_diff = dsqrt(s) - dsqrt(s_cat)
		s_diff = s_diff**2
		if( s <= s_cat ) then
			call fcfac(vib_n,vib_ion,dabs( s_diff ),fc)
		else
			call fcfac(vib_ion,vib_n,dabs( s_diff ),fc)
		endif
		fc_cf(vib_ion, vib_n) = fc
	enddo
	enddo
	
	!anion to frenkel
	do vib_ion=0,vibmax
	do vib_n  =0,vibmax
		s_diff = dsqrt(s) - dsqrt(s_ani)
		s_diff = s_diff**2
		if( s <= s_ani ) then
			call fcfac(vib_n,vib_ion,dabs( s_diff ),fc)
		else
			call fcfac(vib_ion,vib_n,dabs( s_diff ),fc)
		endif
		fc_af(vib_ion, vib_n) = fc
	enddo
	enddo

end subroutine

!********************************************************************
                                    ! frank-condon Factors !
!********************************************************************
subroutine fcfac(n,m,s,fc)
	implicit none
	integer n,m,k
	real*8 s,fc,f_m,f_n,f_k,f_nmk,f_mk,facin
	real*8 fact

	fc = 0.d0

	do k = 0,m
		if(n-m+k < 0) go to 100	! if n-m+k is negative, factorial is not calculatable.

		f_mk  = fact(m-k)
		f_nmk = fact(n-m+k)
		f_k   = fact(k)
		facin = 1.d0/(1.d0*f_k*f_mk*f_nmk)

		fc = fc + facin*s**(k*0.5d0)*s**(1.0d0*(n-m+k)*0.5d0)*(-1)**(n-m+k)
100		continue
	enddo

	f_n = fact(n)
	f_m = fact(m)
!	print*,'f_n,f_m=',f_n,f_m
	fc = fc*dsqrt(1.d0*f_m*f_n)*dexp(-s/2.d0)

	return
end subroutine

!********************************************************************
!	calculating factorial
!********************************************************************
real*8 function fact( n )
	implicit none
	integer	n
	integer i

	fact = 1.0d0
	if( n .ge. 0 ) then
		do i = 2, n
			fact = fact * i
		enddo
	endif
end function


 

!********************************************************************************!
!        This subroutine diagonalizes the hamiltonian to give us the eigen values and eigen vectors
!********************************************************************************!

subroutine dsyevr_diagonalize( a, n, w, m, rrange, iu)
    implicit none
    
    character*1, intent(in):: rrange
    integer, intent(in)   :: n, iu
    integer, intent(out)  :: m
    real*8, intent(inout) :: a(n,n)
    real*8, intent(out)   :: w(n)
    character*1     jobz, uplo
    parameter       ( jobz = 'v', uplo = 'u' )
    
    integer :: lda
    real*8  :: vl = -3.d0    !bounds, but are not referenced
    real*8  :: vu = 20.d0    !bounds, but are not referenced
    integer :: il = 1
    real*8  :: abstol
    real*8, allocatable  :: z(:,:)
    integer :: ldz
    integer isuppz(2*max(1,n))
    real*8, allocatable :: work(:)
    real*8  :: workdim(1)
    integer :: lwork
    integer, allocatable :: iwork(:)
    integer :: iworkdim(1)
    integer :: liwork
    integer :: info
    real*8  time_start, time_finish
    
    real*8, external :: dlamch

    abstol = dlamch( 'safe minimum')    !allow for high relative accuracy of eigenvalues
    lda = n
    ldz = n
    allocate( z( n, n ) )
    !==============================================!
    !    query for work dimensions
    !==============================================!
    lwork =  -1
    liwork = -1
    call dsyevr( jobz, rrange, uplo, n, a, &
                 lda, vl, vu, il, iu,      &
                 abstol, m, w, z, ldz,     &
                 isuppz, workdim, lwork,      &
                 iworkdim, liwork, info )

    lwork = workdim(1)
    liwork = iworkdim(1)
    
    allocate( work( lwork ) )
    allocate( iwork( liwork ) )
    
    !==============================================!
    !    ask lapack to fine eigenspectrum
    !==============================================!
!    print*, 'the size is :: ', n
!    print*, '//                call to dsyevr               \\'
    call cpu_time( time_start )
    call dsyevr( jobz, rrange, uplo, n, a, &
                 lda, vl, vu, il, iu,      &
                 abstol, m, w, z, ldz,     &
                 isuppz, work, lwork,      &
                 iwork, liwork, info )
    call cpu_time( time_finish )
!    print*, '\\                     done                    //'
!    print'(a, f14.4, a)', ' diagonalization time is :: ', time_finish - time_start, ' seconds'
    !reassign kount for abs for number of found eigenvalues
!    print*, 'found:: ', m, ' eigenstates.'
    a(:,1:m) = z(:,1:m)

    !deallocate allocated so it can come back and allocate next time
    deallocate( z, work, iwork )

end subroutine   


!*************************************************!
!        alocate h and eval
!*************************************************!
subroutine allocate_hev()
    use common_variables
    implicit none

    if ( .not. allocated( h ) ) then
        allocate( h( kount, kount ) )
    end if
    if ( .not. allocated( eval ) ) then
        allocate( eval( kount ) )
    end if
    
    h=0.d0
    eval = 0.d0
        
end subroutine

