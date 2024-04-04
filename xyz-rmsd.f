!---------------------------------------------------------------------!
!  xyz-align:                                                         !
!  program so align one set of coordinates onto another               !
!  prerequisite - you coordinates have to be in .xyz format and they  !
!                 MUST have the atoms listed in the same order        !
!                 This code assumes that atom_1 in the first system   !
!                 is the same as atom_1 in the second system          !
!  usage - if you have two .xyz files (say, 1.xyz and 2.xyz), run the !
!          code with:                                                 !
!          cat 1.xyz 2.xyz | xyz-align > out.ani                      !
!  output - standard output will contain three xyz coordinates:       !
!           1. the first system (centred)                             !
!           2. the second system (centred)                            !
!           3. the second system (aligned with the first)             !
!           along with some general system info in the comment lines  !
!---------------------------------------------------------------------!
      module mytype
      type prod
        double precision :: xx, xy, xz
        double precision :: yx, yy, yz
        double precision :: zx, zy, zz
        double precision :: xx2, yy2, zz2
        double precision :: xy2, yz2, xz2
        double precision :: yx2, zy2, zx2
      endtype prod
      endmodule mytype
      module functions
      use mytype
      character(2),     dimension(118) :: al=(/" H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", &
      & "Na", "Mg", "Al", "Si", " P", " S", "Cl", "Ar", " K", "Ca", "Sc", "Ti", " V", "Cr", "Mn", &
      & "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", &
      & "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Ca", "In", "Sn", "Sb", "Te", " I", "Xe", "Cs", &
      & "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
      & "Lu", "Hf", "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", &
      & "Rn", "Fr", "Ra", "Ac", "Th", "Pa", " U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
      & "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", &
      & "Lv", "Ts", "Og"/)
      double precision, dimension(118) :: am=(/1.0080, 4.0026, 6.9400, 9.0122, 10.8100, 12.0110, 14.0070, &
      & 015.9990, 018.9984, 020.1797, 022.9898, 024.3050, 026.9815, 028.0850, 030.9738, 032.0600, 035.4500, & 
      & 039.9500, 039.0983, 040.0780, 044.9559, 047.8670, 050.9415, 051.9961, 054.9380, 055.8450, 058.9332, & 
      & 058.6934, 063.5460, 065.3800, 069.7230, 072.6300, 074.9216, 078.9710, 079.9040, 083.7980, 085.4678, &
      & 087.6200, 088.9058, 091.2240, 092.9064, 095.9500, 097.0000, 101.0700, 102.9055, 106.4200, 107.8682, &
      & 112.4140, 114.8180, 118.7100, 121.7600, 127.6000, 126.9045, 131.2930, 132.9055, 137.3270, 138.9055, &
      & 140.1160, 140.9077, 144.2420, 145.0000, 150.3600, 151.9640, 157.2500, 158.9254, 162.5000, 164.9303, &
      & 167.2590, 168.9342, 173.0450, 174.9668, 178.4860, 180.9479, 183.8400, 186.2070, 190.2300, 192.2170, &
      & 195.0840, 196.9666, 200.5920, 204.3800, 207.2000, 208.9804, 209.0000, 210.0000, 222.0000, 223.0000, &
      & 226.0000, 227.0000, 232.0377, 231.0359, 238.0289, 237.0000, 244.0000, 243.0000, 247.0000, 247.0000, &
      & 251.0000, 252.0000, 257.0000, 258.0000, 259.0000, 262.0000, 267.0000, 270.0000, 269.0000, 270.0000, &
      & 270.0000, 278.0000, 281.0000, 281.0000, 285.0000, 286.0000, 289.0000, 289.0000, 293.0000, 293.0000, &
      & 294.0000/)
      double precision                 :: deg=57.29577951
      type(prod)                       :: S
      double precision                 :: AG, BG, CG, lambda_final
      type coor
        integer                                     :: n
        double precision, allocatable, dimension(:) :: x,y,z,w
        character(2),     allocatable, dimension(:) :: l
      contains
        procedure :: alloc      => coor_allocate
        procedure :: centre     => coor_set_at_com
        procedure :: G          => coor_getG
        procedure :: coor_getS
        generic   :: operator(.getS.)  => coor_getS
        procedure :: coor_get_eig
        generic   :: operator(.getE.)  => coor_get_eig
      endtype coor
      contains      
!---------------------------------------------------------------------!
!  allocate the coor data type arrays                                 !       
!---------------------------------------------------------------------!
      subroutine coor_allocate(A)
      implicit none
      class(coor) :: A
      allocate(A%x(A%n),A%y(A%n),A%z(A%n),A%w(A%n),A%l(A%n))
      endsubroutine coor_allocate
!---------------------------------------------------------------------!
!  calculate the inner product matrix of a group of atoms             !       
!---------------------------------------------------------------------!
      function coor_getS(A,B) result(S_out)
      implicit none
      class(coor), intent(in) :: A,B
      type(prod)              :: S_out
      double precision, dimension(A%n) :: w_dum
      w_dum     = A%w * A%x
      S_out%xx  = sum(w_dum*B%x)
      S_out%xy  = sum(w_dum*B%y)
      S_out%xz  = sum(w_dum*B%z)
      w_dum     = A%w * A%y
      S_out%yx  = sum(w_dum*B%x)
      S_out%yy  = sum(w_dum*B%y)
      S_out%yz  = sum(w_dum*B%z)
      w_dum     = A%w * A%z
      S_out%zx  = sum(w_dum*B%x)
      S_out%zy  = sum(w_dum*B%y)
      S_out%zz  = sum(w_dum*B%z)
      S_out%xx2 = S_out%xx**2
      S_out%yy2 = S_out%yy**2
      S_out%zz2 = S_out%zz**2
      S_out%xy2 = S_out%xy * S_out%xy; S_out%yx2 = S_out%yx * S_out%yx;
      S_out%yz2 = S_out%yz * S_out%yz; S_out%zy2 = S_out%zy * S_out%zy;
      S_out%xz2 = S_out%xz * S_out%xz; S_out%zx2 = S_out%zx * S_out%zx;
      endfunction coor_getS
!---------------------------------------------------------------------!
!  calculate the center of mass of a group of atoms                   !       
!---------------------------------------------------------------------!
      subroutine coor_set_at_com(A)
      implicit none
      class(coor) :: A
      integer     :: i
      double precision :: mass, comx, comy, comz
      mass = sum(A%w)
      comx = sum(A%w * A%x) / mass
      comy = sum(A%w * A%y) / mass
      comz = sum(A%w * A%z) / mass
      A%x  = A%x - comx
      A%y  = A%y - comy
      A%z  = A%z - comz
      endsubroutine coor_set_at_com
!---------------------------------------------------------------------!
!  get the atomic mass for an atom label                              !       
!---------------------------------------------------------------------!
      function get_mass(lab) result(mass)
      implicit none
      character        :: lab
      double precision :: mass
      integer          :: i
      mass=0
      do i=1,10
       if(trim(adjustl(al(i))).eq.trim(adjustl(lab))) mass=am(i)
      enddo
      if(mass.eq.0) stop 'oops: one of your atoms is broken'
      endfunction get_mass
!---------------------------------------------------------------------!
!  get the inner product of a structure                               !       
!---------------------------------------------------------------------!
      function coor_getG(A) result(G)
      implicit none
      class(coor) :: A
      double precision :: G
      G = sum(A%w*(A%x**2+A%y**2+A%z**2))
      endfunction coor_getG
!---------------------------------------------------------------------!
!  get eigenvalue via Newton-Raphson QCP algorithm                    !       
!---------------------------------------------------------------------!
      function coor_get_eig(A,B) result(lambda)
      implicit none
      class(coor), intent(in) :: A, B
! Characteristic Equation Stuff
      double precision        :: A0, A1, A2
      double precision        :: xzpzx, xypyx, yzpzy, xxpyy
      double precision        :: xzmzx, xymyx, yzmzy, xxmyy
      double precision        :: C0, C1, C2
! Newton-Raphson stuff
      double precision        :: l2, la, lb, ld
      double precision        :: lambda, lambda_old, zero=1d-11
! set up some coefficients for the characteristic equation
      A0 = 2.0 * (S%yz * S%zy - S%yy * S%zz)
      A1 = S%xy2 + S%xz2 - S%yx2 - S%zx2
      A2 = S%yy2 + S%zz2 - S%xx2 + S%yz2 + S%zy2
      xzpzx = S%xz + S%zx; xzmzx = S%xz - S%zx
      xypyx = S%xy + S%yx; xymyx = S%xy - S%yx
      yzpzy = S%yz + S%zy; yzmzy = S%yz - S%zy
      xxpyy = S%xx + S%yy; xxmyy = S%xx - S%yy
      C0 = A1**2 + (A2 + A0) * (A2 - A0) &
      &  + (-(xzpzx * yzmzy) + (xymyx * (xxmyy - S%zz))) * (-(xzmzx * yzpzy) + (xymyx * (xxmyy + S%zz))) &
      &  + (-(xzpzx * yzpzy) - (xypyx * (xxpyy - S%zz))) * (-(xzmzx * yzmzy) - (xypyx * (xxpyy + S%zz))) &
      &  + (+(xypyx * yzpzy) + (xzpzx * (xxmyy + S%zz))) * (-(xymyx * yzmzy) + (xzpzx * (xxpyy + S%zz))) &
      &  + (+(xypyx * yzmzy) + (xzmzx * (xxmyy - S%zz))) * (-(xymyx * yzpzy) + (xzmzx * (xxpyy - S%zz)))
      C1 = 8.d0 * (S%xx * (S%yz * S%zy - S%yy * S%zz) + S%zx * (S%yy * S%xz - S%yz * S%xy) + S%yx * (S%zz * S%xy - S%zy * S%xz))
      C2 =-2.d0 * (S%xx2 + S%yy2 + S%zz2 + S%xy2 + S%yx2 + S%xz2 + S%zx2 + S%yz2 + S%zy2)
! perform iterative Newton-Raphson QCP algorithm:
      lambda_old = AG + BG
      lambda     = lambda_old * 0.5d0
      do while(abs(lambda-lambda_old).gt.zero*lambda)  
        lambda_old = lambda
        l2         = lambda**2
        lb         = (l2 + C2) * lambda;
        la         = lb + C1;
        ld         = ((la * lambda + C0)/(2.d0 * l2 * lambda + lb + la))
        lambda     = lambda_old - ld
      enddo
      endfunction coor_get_eig
!---------------------------------------------------------------------!
!  get the determinant of the rotation matrix                         !       
!---------------------------------------------------------------------!      
      function get_det(A) result(det)
      implicit none
      double precision, dimension(3,3) :: A
      double precision                 :: det
      det = A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2)) &
         & - A(1,2) * (A(2,1) * A(3,3) - A(2,3) * A(3,1)) &
          & + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1)) 
      endfunction get_det
!---------------------------------------------------------------------!
!  get rotation matrix from the unit quaternion                       !       
!---------------------------------------------------------------------!
      function coor_get_rot() result(R)
      implicit none
! K matrix stuff
      double precision                 :: k01, k02, k03, k12, k13, k23
      double precision                 :: k11, k22, k33
      double precision                 :: k2233_3223, k2133_3123, k2132_3122
      double precision                 :: k2032_3022, k2033_3023, k2031_3021
! Quaternion stuff
      double precision                 :: nrm
      double precision, dimension(0:3) :: q
! Rotation Matrix stuff
      double precision                 :: xy, wz, zx, wy, yz, wx, diag
      double precision, dimension(3,3) :: R
! build the elements of the K matrix that we need
      k11 = + S%xx - S%yy - S%zz - lambda_final
      k22 = - S%xx + S%yy - S%zz - lambda_final
      k33 = - S%xx - S%yy + S%zz - lambda_final
      k01 =   S%yz - S%zy; k23 =   S%yz + S%zy
      k02 =   S%zx - S%xz; k13 =   S%zx + S%xz
      k03 =   S%xy - S%yx; k12 =   S%xy + S%yx
! find the adjugate and, thus, the quaternion
      k2233_3223 = k22 * k33 - k23 * k23; k2133_3123 = k12 * k33 - k13 * k23
      k2132_3122 = k12 * k23 - k13 * k22; k2032_3022 = k02 * k23 - k03 * k22
      k2033_3023 = k02 * k33 - k03 * k23; k2031_3021 = k02 * k13 - k03 * k12
      q(0) =   k11 * k2233_3223 - k12 * k2133_3123 + k13 * k2132_3122
      q(1) = - k01 * k2233_3223 + k12 * k2033_3023 - k13 * k2032_3022
      q(2) =   k01 * k2133_3123 - k11 * k2033_3023 + k13 * k2031_3021
      q(3) = - k01 * k2132_3122 + k11 * k2032_3022 - k12 * k2031_3021
      nrm  = sum(q**2)
      q    = q / sqrt(nrm)
! build the rotation matrix from the quaternion
      xy     = q(1) * q(2); wz     = q(0) * q(3);
      zx     = q(3) * q(1); wy     = q(0) * q(2);
      yz     = q(2) * q(3); wx     = q(0) * q(1);
      diag   = 2 * q(0)**2 -1
      R(:,1) = (/2 * q(1)**2 + diag, 2 * (xy - wz),      2 * (zx + wy)/)
      R(:,2) = (/2 * (xy + wz),      2 * q(2)**2 + diag, 2 * (yz - wx)/)
      R(:,3) = (/2 * (zx - wy),      2 * (yz + wx),      2 * q(3)**2 + diag/)
      endfunction coor_get_rot
      endmodule functions
!---------------------------------------------------------------------!
!  main program                                                       !       
!---------------------------------------------------------------------!
      program xyz_rmsd
      use functions
      implicit none
      type(coor)                       :: A,B,C
      double precision                 :: rmsd, theta
      double precision, dimension(3,3) :: R,R1,R2
      integer                          :: i,j
      character(100)                   :: dummy
      logical                          :: do_align
      call moses(do_align)
      read(*,*) A%n
      call A%alloc
      read(*,*) dummy
      do i=1,A%n; read(*,*) A%l(i),A%x(i),A%y(i),A%z(i); enddo
      do i=1,A%n; A%w(i)=get_mass(A%l(i)); enddo
      call A%centre()
      AG = A%G()
      read(*,*) B%n
      if(A%n.ne.B%n) stop 'ERROR in n'
      call B%alloc
      read(*,*) dummy
      do i=1,B%n; read(*,*) B%l(i),B%x(i),B%y(i),B%z(i); enddo
      B%w = A%w
      call B%centre()
      BG = B%G()
      C%n = A%n; call C%alloc
      if(do_align) then
       write(*,*) A%n
       write(*,100) AG
       do i=1,A%n
         write(*,200) A%l(i), A%x(i), A%y(i), A%z(i)
       enddo
       write(*,*) A%n
       write(*,100) BG
       do i=1,A%n
         write(*,200) B%l(i), B%x(i), B%y(i), B%z(i)
       enddo
      endif
      S            = A.getS.B
      lambda_final = A.getE.B
      R            = coor_get_rot()
      rmsd         = sqrt( abs(AG + BG - 2 * lambda_final) / A%n )
      theta        = acos((R(1,1)+R(2,2)+R(3,3)-1)/2.d0) * deg
      if(do_align) then
       C%x          = R(1,1) * B%x + R(1,2) * B%y + R(1,3) * B%z
       C%y          = R(2,1) * B%x + R(2,2) * B%y + R(2,3) * B%z
       C%z          = R(3,1) * B%x + R(3,2) * B%y + R(3,3) * B%z
       write(*,*) A%n
       write(*,300) rmsd, theta, get_det(R)
       do i=1,A%n
         write(*,200) A%l(i), C%x(i), C%y(i), C%z(i)
       enddo
      else
       write(*,300) rmsd, theta, get_det(R)
      endif
100   format(t1,"G = ",t5,f15.4)
200   format(t1,a,t5,f8.4,t15,f8.4,t25,f8.4)
300   format(t1,"RMSD = ",t8,f8.4,t18,"; Angle = ",t28,f8.4,t38,"; DetR = ",t47,f8.4)
      endprogram xyz_rmsd
!---------------------------------------------------------------------!
      subroutine moses(do_align)
      implicit none
      integer                                :: narg, i
      character(50),allocatable,dimension(:) :: arg
      logical                                :: do_align
      do_align = .false.
      narg = iargc()
      allocate(arg(narg))
      do i=1,narg
        call getarg(i,arg(i))
        if(trim(arg(i)).eq."-h".or.trim(arg(i)).eq."--help") then
         write(*,100)
         write(*,200)
         write(*,300) "Find the RMSD and angle between two geometries"
         write(*,200)
         write(*,300) "Usage - Consise Mode:"
         write(*,300) "       if you have two .xyz files (say, 1.xyz and 2.xyz), run the"
         write(*,300) "       code with:"
         write(*,200)
         write(*,300) "       cat 1.xyz 2.xyz | xyz-rmsd"
         write(*,200)
         write(*,300) "Usage - Alignment Mode:"
         write(*,300) "       if you have two .xyz files (say, 1.xyz and 2.xyz), run the"
         write(*,300) "       code with:"
         write(*,200)
         write(*,300) "       cat 1.xyz 2.xyz | xyz-align -a"
         write(*,200)
         write(*,300) "Requirements: Your coordinates have to be in .xyz format and they"
         write(*,300) "              MUST have the atoms listed in the same order"
         write(*,300) "              This code assumes that atom_1 in the first system"
         write(*,300) "              is the same as atom_1 in the second system"
         write(*,200)
         write(*,300) "Output - Consise Mode:" 
         write(*,300) "         RMSD, Angle and the determinant of the rotation matrix"
         write(*,200)
         write(*,300) "Output - Alignment Mode:"
         write(*,300) "         Standard output will contain three sets of xyz coordinates:"
         write(*,300) "         1. the first system (centred)"
         write(*,300) "         2. the second system (centred)"
         write(*,300) "         3. the second system (aligned with the first)"
         write(*,300) "         along with some general system info in the comment lines"
         write(*,200)
         write(*,100)
         write(*,200)
         write(*,300) "Author:"
         write(*,300) "Thomas Pope"
         write(*,300) "thomas.pope2@newcastle.ac.uk"
         write(*,300) "Chemistry, School of Natural and Environmental Sciences,"
         write(*,300) "Newcastle University, NE1 7RU, Newcastle Upon Tyne, UK"
         write(*,200)
         write(*,100)
         stop
        elseif(trim(arg(i)).eq."-a") then
         do_align = .true.
        endif
       enddo
100    format("+----------------------------------------------------------------------+")
200    format(t1,"|",t72,"|")
300    format(t1,"|",1x,a,t72,"|")
      endsubroutine moses

      
