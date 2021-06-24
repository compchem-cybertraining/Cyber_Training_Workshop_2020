

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
       222822400 of real*8 words ( 1700.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-5,
   tol(2)=1.e-5,
   tol(1)=1.e-8,
   NSTATE=0,
   npath=1,3,9,10,13,17,19,21,-11,12, 2,
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  1,7,20,1,8,20,1,9,20
   NAVST(1) = 3,
   WAVST(1,1)=1 ,
   WAVST(1,2)=1 ,
   WAVST(1,3)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /panfs/panfs.cbls.ccr.buffalo.edu/scratch/grp-cyberwksp21
 /el

 Integral file header information:
 Hermit Integral Program : SIFS version  srv-p22-13.cbls.c 15:10:26.633 22-Jun-21

 Core type energy values:
 energy( 1)=  3.858870424817E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   38.588704248


   ******  Basis set information:  ******

 Number of irreps:                  1
 Total number of basis functions:  36

 irrep no.              1
 irrep label           A  
 no. of bas.fcions.    36


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-08. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-05. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-05. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.0000E-06. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  1 DRT.

 DRT  first state   no.of aver.states   weights
  1   ground state          3             0.333 0.333 0.333

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        4

 orbital coefficients are optimized for the ground state (nstate=0).

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       7(  7)    20
       1       8(  8)    20
       1       9(  9)    20

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a  
 Total number of electrons:   16
 Spin multiplicity:            1
 Number of active orbitals:    3
 Number of active electrons:   4
 Total number of CSFs:         6
 

 faar:   0 active-active rotations allowed out of:   3 possible.


 Number of active-double rotations:        18
 Number of active-active rotations:         0
 Number of double-virtual rotations:      162
 Number of active-virtual rotations:       81
 lenbfsdef=                131071  lenbfs=                   729
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #          21
 class  2 (pq|ri):         #         108
 class  3 (pq|ia):         #         972
 class  4 (pi|qa):         #        1458
 class  5 (pq|ra):         #         486
 class  6 (pq|ij)/(pi|qj): #         342
 class  7 (pq|ab):         #        2268
 class  8 (pa|qb):         #        4374
 class  9 p(bp,ai)         #       13122
 class 10p(ai,jp):        #        2916
 class 11p(ai,bj):        #       15309

 Size of orbital-Hessian matrix B:                    37395
 Size of the orbital-state Hessian matrix C:           4698
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:          42093


 Source of the initial MO coeficients:

 Input MO coefficient file: /panfs/panfs.cbls.ccr.buffalo.edu/scratch/grp-cyberwksp21/el
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58622 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     1
 ciiter=   4 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3970483261     -132.9857525743        0.0000000000        0.0000010000
    2       -94.0395145268     -132.6282187750        0.0000000000        0.0000010000
    3       -94.0086169091     -132.5973211573        0.0000000000        0.0000010000
    4       -93.7330751539     -132.3217794020        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.825368875234238E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.98885683 pnorm= 0.0000E+00 rznorm= 1.7310E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.839774  -23.209412   -3.117436   -2.365828   -2.063345   -1.951360

 qvv(*) eigenvalues. symmetry block  1
    -0.046699    0.029521    0.082826    0.260919    0.398846    0.993398    0.999124    1.021952    1.181037    1.393786
     1.394987    1.584273    1.797026    1.811615    1.841617    2.469267    2.678596    2.815955    2.910274    3.591675
     3.769080    4.093437    4.393313    4.752156    4.857009    5.481647    5.564782

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.118 walltime=     0.595

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -94.1483932540 demc= 9.4148E+01 wnorm= 1.4603E-01 knorm= 1.4887E-01 apxde= 5.9151E-03    *not conv.*     

               starting mcscf iteration...   2
 !timer:                                 cpu_time=     0.128 walltime=     0.637

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58622 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3982437836     -132.9869480317        0.0000000000        0.0000100000
    2       -94.0523459828     -132.6410502310        0.0000000000        0.0000100000
    3       -94.0141078220     -132.6028120701        0.0000000000        0.0000100000
    4       -93.7596461045     -132.3483503526        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.178916636609998E-003
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99836019 pnorm= 0.0000E+00 rznorm= 6.5429E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.799707  -23.167099   -3.093061   -2.349542   -2.044295   -1.932803

 qvv(*) eigenvalues. symmetry block  1
    -0.037905    0.038888    0.085783    0.266883    0.405149    0.999084    1.010023    1.032703    1.187397    1.406100
     1.410939    1.598344    1.808564    1.831042    1.854801    2.485495    2.688993    2.834971    2.930068    3.609634
     3.786795    4.111322    4.413251    4.770150    4.874389    5.503625    5.583601

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=    -94.1548991961 demc= 6.5059E-03 wnorm= 1.7431E-02 knorm= 5.7244E-02 apxde= 2.4315E-04    *not conv.*     

               starting mcscf iteration...   3
 !timer:                                 cpu_time=     0.184 walltime=     0.928

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58619 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3983140646     -132.9870183128        0.0000000000        0.0000010000
    2       -94.0534919706     -132.6421962188        0.0000000000        0.0000010000
    3       -94.0139316850     -132.6026359332        0.0000000000        0.0000010000
    4       -93.7671404724     -132.3558447206        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.666757514520606E-003
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99976661 pnorm= 0.0000E+00 rznorm= 2.9981E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.796336  -23.165527   -3.091190   -2.348597   -2.038938   -1.932038

 qvv(*) eigenvalues. symmetry block  1
    -0.036639    0.039779    0.085205    0.266848    0.404980    0.998103    1.011226    1.032849    1.186522    1.407394
     1.413180    1.598501    1.808906    1.832683    1.855919    2.486818    2.689661    2.836142    2.931120    3.610537
     3.787528    4.111977    4.414234    4.771109    4.875438    5.505218    5.585043

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=    -94.1552459067 demc= 3.4671E-04 wnorm= 1.3334E-02 knorm= 2.1604E-02 apxde= 4.6411E-05    *not conv.*     

               starting mcscf iteration...   4
 !timer:                                 cpu_time=     0.240 walltime=     1.165

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58617 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3984326560     -132.9871369041        0.0000000000        0.0000010000
    2       -94.0535923194     -132.6422965676        0.0000000000        0.0000010000
    3       -94.0139140105     -132.6026182587        0.0000000000        0.0000010000
    4       -93.7697644840     -132.3584687321        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.814357123816976E-004
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99995776 pnorm= 0.0000E+00 rznorm= 1.2939E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.793104  -23.167754   -3.090448   -2.348550   -2.036118   -1.932145

 qvv(*) eigenvalues. symmetry block  1
    -0.036189    0.040016    0.084783    0.266663    0.404698    0.997797    1.011407    1.032617    1.185813    1.408365
     1.414020    1.597961    1.808653    1.833515    1.856299    2.487432    2.689843    2.836383    2.931222    3.610619
     3.787422    4.111750    4.414242    4.771272    4.875575    5.505602    5.585415

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=    -94.1553129953 demc= 6.7089E-05 wnorm= 6.2515E-03 knorm= 9.1917E-03 apxde= 8.9823E-06    *not conv.*     

               starting mcscf iteration...   5
 !timer:                                 cpu_time=     0.300 walltime=     1.473

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58619 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3984894709     -132.9871937191        0.0000000000        0.0000010000
    2       -94.0535843606     -132.6422886088        0.0000000000        0.0000010000
    3       -94.0139035264     -132.6026077746        0.0000000000        0.0000010000
    4       -93.7707533161     -132.3594575643        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.393641658690922E-004
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999274 pnorm= 0.0000E+00 rznorm= 5.3208E-07 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.791579  -23.168902   -3.090132   -2.348550   -2.034798   -1.932212

 qvv(*) eigenvalues. symmetry block  1
    -0.036003    0.040110    0.084599    0.266576    0.404567    0.997646    1.011464    1.032500    1.185496    1.408815
     1.414370    1.597695    1.808512    1.833881    1.856472    2.487697    2.689916    2.836464    2.931245    3.610637
     3.787350    4.111623    4.414221    4.771328    4.875615    5.505750    5.585561

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=    -94.1553257860 demc= 1.2791E-05 wnorm= 2.7149E-03 knorm= 3.8094E-03 apxde= 1.5969E-06    *not conv.*     

               starting mcscf iteration...   6
 !timer:                                 cpu_time=     0.357 walltime=     1.737

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58622 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985144890     -132.9872187371        0.0000000000        0.0000010000
    2       -94.0535717825     -132.6422760306        0.0000000000        0.0000010000
    3       -94.0138978619     -132.6026021100        0.0000000000        0.0000010000
    4       -93.7711414647     -132.3598457129        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.420718596791048E-004
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999878 pnorm= 0.0000E+00 rznorm= 2.1741E-07 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790903  -23.169411   -3.089995   -2.348550   -2.034225   -1.932240

 qvv(*) eigenvalues. symmetry block  1
    -0.035926    0.040148    0.084522    0.266539    0.404512    0.997581    1.011487    1.032451    1.185362    1.409011
     1.414517    1.597580    1.808450    1.834041    1.856551    2.487812    2.689948    2.836499    2.931254    3.610645
     3.787319    4.111568    4.414211    4.771353    4.875633    5.505813    5.585624

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    6 emc=    -94.1553280444 demc= 2.2584E-06 wnorm= 1.1366E-03 knorm= 1.5625E-03 apxde= 2.7284E-07    *not conv.*     

               starting mcscf iteration...   7
 !timer:                                 cpu_time=     0.416 walltime=     1.961

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58619 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985250221     -132.9872292703        0.0000000000        0.0000010000
    2       -94.0535650125     -132.6422692607        0.0000000000        0.0000010000
    3       -94.0138952530     -132.6025995011        0.0000000000        0.0000010000
    4       -93.7712970275     -132.3600012757        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  5.854585049740314E-005
 Total number of micro iterations:    5

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999980 pnorm= 0.0000E+00 rznorm= 8.0514E-07 rpnorm= 0.0000E+00 noldr=  5 nnewr=  5 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790617  -23.169626   -3.089938   -2.348550   -2.033985   -1.932252

 qvv(*) eigenvalues. symmetry block  1
    -0.035894    0.040165    0.084490    0.266524    0.404489    0.997554    1.011496    1.032430    1.185307    1.409093
     1.414578    1.597532    1.808423    1.834108    1.856584    2.487861    2.689961    2.836514    2.931258    3.610648
     3.787306    4.111545    4.414207    4.771363    4.875640    5.505840    5.585650

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    7 emc=    -94.1553284292 demc= 3.8476E-07 wnorm= 4.6837E-04 knorm= 6.3822E-04 apxde= 4.5828E-08    *not conv.*     

               starting mcscf iteration...   8
 !timer:                                 cpu_time=     0.478 walltime=     2.271

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58618 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985293699     -132.9872336181        0.0000000000        0.0000010000
    2       -94.0535619709     -132.6422662191        0.0000000000        0.0000010000
    3       -94.0138941404     -132.6025983886        0.0000000000        0.0000010000
    4       -93.7713599740     -132.3600642221        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.396440214013159E-005
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    5

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999990 pnorm= 4.1515E-04 rznorm= 5.4183E-07 rpnorm= 1.7244E-07 noldr=  5 nnewr=  5 nolds=  2 nnews=  2
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790499  -23.169715   -3.089914   -2.348550   -2.033886   -1.932256

 qvv(*) eigenvalues. symmetry block  1
    -0.035881    0.040171    0.084477    0.266518    0.404479    0.997543    1.011500    1.032421    1.185284    1.409127
     1.414603    1.597512    1.808412    1.834136    1.856597    2.487881    2.689967    2.836520    2.931259    3.610650
     3.787301    4.111535    4.414206    4.771368    4.875643    5.505851    5.585661

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    8 emc=    -94.1553284937 demc= 6.4550E-08 wnorm= 1.9172E-04 knorm= 4.3945E-04 apxde= 1.2908E-08    *not conv.*     

               starting mcscf iteration...   9
 !timer:                                 cpu_time=     0.541 walltime=     2.616

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58618 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   5 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985323829     -132.9872366311        0.0000000000        0.0000010000
    2       -94.0535597848     -132.6422640330        0.0000000000        0.0000010000
    3       -94.0138933522     -132.6025976003        0.0000000000        0.0000010000
    4       -93.7714030991     -132.3601073473        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  6.898808342772562E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 4.7724E-07 rpnorm= 1.8371E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790417  -23.169776   -3.089897   -2.348550   -2.033817   -1.932260

 qvv(*) eigenvalues. symmetry block  1
    -0.035872    0.040176    0.084468    0.266514    0.404473    0.997535    1.011503    1.032415    1.185269    1.409150
     1.414620    1.597498    1.808405    1.834155    1.856607    2.487894    2.689971    2.836524    2.931261    3.610651
     3.787297    4.111529    4.414205    4.771371    4.875645    5.505858    5.585668

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    9 emc=    -94.1553285066 demc= 1.2906E-08 wnorm= 5.5190E-07 knorm= 2.0291E-08 apxde= 5.0162E-15    *not conv.*     

               starting mcscf iteration...  10
 !timer:                                 cpu_time=     0.593 walltime=     2.903

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58618 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   2 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985323833     -132.9872366315        0.0000000000        0.0000010000
    2       -94.0535597855     -132.6422640337        0.0000000000        0.0000010000
    3       -94.0138933511     -132.6025975993        0.0000000000        0.0000010000
    4       -93.7714030978     -132.3601073460        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  5.965241754420209E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 3.4313E-07 rpnorm= 8.2571E-10 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790417  -23.169776   -3.089897   -2.348550   -2.033817   -1.932260

 qvv(*) eigenvalues. symmetry block  1
    -0.035872    0.040176    0.084468    0.266514    0.404473    0.997535    1.011503    1.032415    1.185269    1.409150
     1.414620    1.597498    1.808405    1.834155    1.856607    2.487894    2.689971    2.836524    2.931261    3.610651
     3.787297    4.111529    4.414205    4.771371    4.875645    5.505858    5.585668

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=   10 emc=    -94.1553285066 demc=-2.8422E-14 wnorm= 4.7722E-07 knorm= 1.2821E-08 apxde= 2.0624E-15    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.333 total energy=      -94.398532383, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.333 total energy=      -94.053559786, rel. (eV)=   9.387186
   DRT #1 state # 3 wt 0.333 total energy=      -94.013893351, rel. (eV)=  10.466565
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.66575125     1.63269253
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.70155623     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        105 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
    N1_ s       1.998472  -0.000030   1.398408   0.213111   0.000000  -0.000604
    N1_ p       0.000000  -0.000292   0.041269   0.494174   1.315011   0.644431
    N1_ d       0.000000  -0.000105   0.004674   0.005740   0.003805   0.004039
    C1_ s       0.000539   1.999452   0.324541   0.807084   0.000000   0.160542
    C1_ p       0.000476   0.000003   0.098838   0.034165   0.144495   0.703404
    C1_ d      -0.000081   0.000001   0.010188   0.003498   0.003358   0.014052
    H1_ s       0.000301   0.000003   0.052064   0.096151   0.240417   0.066513
    H2_ s       0.000301   0.000003   0.052065   0.096146   0.240418   0.066518
    H3_ s      -0.000004   0.000483   0.008976   0.124965   0.026250   0.170553
    H4_ s      -0.000004   0.000483   0.008977   0.124964   0.026246   0.170552
 
   ao class       7A         8A         9A        10A        11A        12A  
    N1_ s       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
    N1_ p       0.074543   1.424815   0.073979   0.000000   0.000000   0.000000
    N1_ d       0.011182   0.002547   0.003589   0.000000   0.000000   0.000000
    C1_ s       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
    C1_ p       0.904777   0.192226   0.621845   0.000000   0.000000   0.000000
    C1_ d       0.014492   0.013105   0.002143   0.000000   0.000000   0.000000
    H1_ s       0.041375   0.000000   0.000000   0.000000   0.000000   0.000000
    H2_ s       0.041375   0.000000   0.000000   0.000000   0.000000   0.000000
    H3_ s       0.289000   0.000000   0.000000   0.000000   0.000000   0.000000
    H4_ s       0.289006   0.000000   0.000000   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.609356   3.292158   0.496824   0.496826   0.620224   0.620225
      p         4.067930   2.700229   0.000000   0.000000   0.000000   0.000000
      d         0.035473   0.060756   0.000000   0.000000   0.000000   0.000000
    total       7.712759   6.053142   0.496824   0.496826   0.620224   0.620225
 

 Total number of electrons:   16.00000000

 !timer: mcscf                           cpu_time=     0.674 walltime=     3.263
 *** cpu_time / walltime =      0.207
 bummer (warning):timer: cpu_time << walltime.  If possible, increase core memory. event =1
