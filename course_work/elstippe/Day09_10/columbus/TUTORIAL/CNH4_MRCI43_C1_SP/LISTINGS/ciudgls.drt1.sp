1
 program ciudg      
 multireference single and double excitation configuration
 interaction based on the graphical unitary group approach.


 references:  h. lischka, r. shepard, f. b. brown, and i. shavitt,
                  int. j. quantum chem. s 15, 91 (1981).
              r. shepard, r. a. bair, r. a. eades, a. f. wagner,
                  m. j. davis, l. b. harding, and t. h. dunning,
                  int j. quantum chem. s 17, 613 (1983).
              r. ahlrichs, h.-j. boehm, c. ehrhardt, p. scharf,
                  h. schiffer, h. lischka, and m. schindler,
                  j. comp. chem. 6, 200 (1985).
              r. shepard, i. shavitt, r. m. pitzer, d. c. comeau, m. pepper
                  h. lischka, p. g. szalay, r. ahlrichs, f. b. brown, and
                  j.-g. zhao, int. j. quantum chem. symp. 22, 149 (1988).

 This Version of Program CIUDG is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              CIUDG       **
     **    PROGRAM VERSION:      2009-03.    **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


================ Computing sorting integral file structure ================

                    -----z----- -----y----- -----x----- -----w----- ---total---

                CSFs        28        2484       24921       37044       64477
      internal walks        28          92          71          98         289
valid internal walks        28          92          71          98         289
 getinfoarray: info=                     6 :                     1
                  8192                  6552                  8192
                  5460                     0
 icd(3)=                  1191 ci%nnlev=                   595  l2rec=
                  8192  n2max=                  5460
 lcore1,lcore2=             222785122             222761806
 lencor,maxblo             222822400                 60000
========================================
 current settings:
 minbl3        1131
 minbl4        1131
 locmaxbl3    20362
 locmaxbuf    10181
 maxbl3       60000
 maxbl3       60000
 maxbl4       60000
 maxbuf       30006
========================================

 sorted 4-external integrals:     3 records of integral w-combinations 
                                  3 records of integral x-combinations of length= 30006
                                        have been written (wstat,xstat=  100.00  100.00)

 sorted 3-external integrals:     4 records of integral w-combinations 
                                  3 records of integral x-combinations of length= 30006
                                        have been written (wstat,xstat=   75.00  100.00)
 Orig.  diagonal integrals:  1electron:        34
                             0ext.    :        56
                             2ext.    :       378
                             4ext.    :       756


 Orig. off-diag. integrals:  4ext.    :     81081
                             3ext.    :     76167
                             2ext.    :     31752
                             1ext.    :      6804
                             0ext.    :       546
                             2ext. SO :         0
                             1ext. SO :         0
                             0ext. SO :         0
                             1electron:       378


 Sorted integrals            3ext.  w :     71253 x :     66339
                             4ext.  w :     71253 x :     61425


Cycle #  1 sortfile size=   1572864(       3 records of   524288) #buckets=   2
distributed memory consumption per node=         0 available core 222785122
Cycle #  2 sortfile size=   1572864(       3 records of   524288) #buckets=   2
distributed memory consumption per node=         0 available core 222785122
 minimum size of srtscr:    524288 WP (     1 records)
 maximum size of srtscr:   1572864 WP (     3 records)
diagi   file:      4 records  of   6144 WP each=>      24576 WP total
ofdgi   file:     10 records  of   6144 WP each=>      61440 WP total
fil3w   file:      4 records  of  30006 WP each=>     120024 WP total
fil3x   file:      3 records  of  30006 WP each=>      90018 WP total
fil4w   file:      3 records  of  30006 WP each=>      90018 WP total
fil4x   file:      3 records  of  30006 WP each=>      90018 WP total
 compressed index vector length=                     1
 echo of the input for program ciudg:
 ------------------------------------------------------------------------
  &input
  NTYPE = 0,
  GSET = 0,
   DAVCOR =10,
  NCOREL = 12
  NROOT = 3
  IVMODE = 3
  NBKITR = 1
  NVBKMN = 3
  RTOLBK = 1e-3,1e-3,1e-3,
  NITER = 60
  NVCIMN = 5
  RTOLCI = 1e-3,1e-3,1e-3,
  NVCIMX = 8
  NVRFMX = 8
  NVBKMX = 8
  IDEN  = 1
  CSFPRN = 10,
 /&end
 ------------------------------------------------------------------------
lodens (list->root)=  3
invlodens (root->list)= -1 -1  1
 bummer (warning):resetting fileloc for seriel operation0
 USING SEGMENTS OF EQUAL SIZE

****************  list of control variables  ****************
 lvlprt =    0      nroot  =    3      noldv  =   0      noldhv =   0
 nunitv =    5      nbkitr =    1      niter  =  60      davcor =  10
 csfprn =   10      ivmode =    3      istrt  =   0      vout   =   0
 iortls =    0      nvbkmx =    8      ibktv  =  -1      ibkthv =  -1
 nvcimx =    8      icitv  =   -1      icithv =  -1      frcsub =   0
 nvbkmn =    3      nvcimn =    5      maxseg = 300      nrfitr =  30
 ncorel =   12      nvrfmx =    8      nvrfmn =   5      iden   =   1
 itran  =    0      froot  =    0      rtmode =   0      ncouple=   1
 skipso =    F      dalton2=    0      molcas =   0      finalv =   0
 finalw =    0      cosmocalc=   0    with_tsklst=   0
 nsegwx =    1     1     1     1
 nseg0x =    1     1     1     1
 nseg1x =    1     1     1     1
 nseg2x =    1     1     1     1
 nseg3x =    1     1     1     1
 nseg4x =    1     1     1     1
 no0ex  =      0    no1ex  =      0    no2ex  =     0    no3ex  =     0
 no4ex  =      0    nodiag =      0
 cdg4ex =    1      c3ex1ex=    1      c2ex0ex=   1
 fileloc=    0     0     0     0     0     0     0     1     1     1
 directhd=   1      noaqccshift_zyxw=      0
 critical_crit=-1.00000    critical_delta= 0.05000

 ctol   = 0.010000    lrtshift=1.000000    smalld =0.001000


 convergence tolerances of bk and full diagonalization steps
 root #       rtolbk        rtol
 ------      --------      ------
    1        1.000E-03    1.000E-03
    2        1.000E-03    1.000E-03
    3        1.000E-03    1.000E-03
 Computing density:                    .drt1.state3
 Main memory management:
 global                1 DP per process
 vdisk                 0 DP per process
 stack                 0 DP per process
 core          222822399 DP per process
 gapointer%node_offset(*)=                     0                     0
 gapointer%node_width(*)=                     0                     0

********** Integral sort section *************


 workspace allocation information: lencor= 222822399

 echo of the input for program cisrt:
 ------------------------------------------------------------------------
  &input
  maxbl3=60000
  maxbl4=60000
  &end
 ------------------------------------------------------------------------
 
 ( 6) listing file:                    ciudgls             
 ( 5) input file:                      cisrtin   
 (17) cidrt file:                      cidrtfl             
 (11) transformed integrals file:      moints    
 (12) diagonal integral file:          diagint             
 (13) off-diagonal integral file:      ofdgint             
 (31) 4-external w integrals file:     fil4w               
 (32) 4-external x integrals file:     fil4x               
 (33) 3-external w integrals file:     fil3w               
 (34) 3-external x integrals file:     fil3x               
 (21) scratch da sorting file:         srtscr              
 (12) 2-e integral file [fsplit=2]:    moints2   

 input integral file header information:
 Hermit Integral Program : SIFS version  srv-p22-13.cbls.c 12:43:48.612 24-Jun-21
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =   38.588704248                                          
 MCSCF energy =     -94.155328507                                                
 SIFS file created by program tran.      srv-p22-13.cbls.c 12:43:50.919 24-Jun-21

 input energy(*) values:
 energy( 1)=  3.858870424817E+01, ietype=   -1,    core energy of type: Nuc.Rep.

 total core energy =   3.858870424817E+01

 nsym = 1 nmot=  34

 symmetry  =    1
 slabel(*) =  A  
 nmpsy(*)  =   34

 info(*) =          1      8192      6552      8192      5460         0

 orbital labels, i:molab(i)=
   1:tout:001   2:tout:002   3:tout:003   4:tout:004   5:tout:005   6:tout:006   7:tout:007   8:tout:008   9:tout:009  10:tout:010
  11:tout:011  12:tout:012  13:tout:013  14:tout:014  15:tout:015  16:tout:016  17:tout:017  18:tout:018  19:tout:019  20:tout:020
  21:tout:021  22:tout:022  23:tout:023  24:tout:024  25:tout:025  26:tout:026  27:tout:027  28:tout:028  29:tout:029  30:tout:030
  31:tout:031  32:tout:032  33:tout:033  34:tout:034

 input parameters:
 prnopt=  0
 ldamin=   32767 ldamax=  524288 ldainc=    4096
 maxbuf=   30006 maxbl3=   60000 maxbl4=   60000 intmxo=    6144
  Using 32 bit compression 

 drt information:
  cidrt_title                                                                    
 nmotd =  36 nfctd =   2 nfvtc =   0 nmot  =  34
 nlevel =  34 niot  =   7 lowinl=  28
 orbital-to-level map(*)
   -1  -1  28  29  30  31  32  33  34   1   2   3   4   5   6   7   8   9  10  11
   12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
 compressed map(*)
   28  29  30  31  32  33  34   1   2   3   4   5   6   7   8   9  10  11  12  13
   14  15  16  17  18  19  20  21  22  23  24  25  26  27
 levsym(*)
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
    1   1   1   1   1   1   1   1   1   1   1   1   1   1
 repartitioning mu(*)=
   2.  2.  2.  2.  0.  0.  0.

 new core energy added to the energy(*) list.
 from the integral file: h1_core= -9.226324016262E+01

 indxdg: diagonal integral statistics.
 total number of integrals contributing to diagonal matrix elements:      1190
 number with all external indices:       756
 number with half external - half internal indices:       378
 number with all internal indices:        56

 indxof: off-diagonal integral statistics.
    4-external integrals: num=      81081 strt=          1
    3-external integrals: num=      76167 strt=      81082
    2-external integrals: num=      31752 strt=     157249
    1-external integrals: num=       6804 strt=     189001
    0-external integrals: num=        546 strt=     195805

 total number of off-diagonal integrals:      196350


 indxof(2nd)  ittp=   3 numx(ittp)=       31752
 indxof(2nd)  ittp=   4 numx(ittp)=        6804
 indxof(2nd)  ittp=   5 numx(ittp)=         546

 intermediate da file sorting parameters:
 nbuk=   2 lendar=  524288 nipbk=  349524 nipsg= 221827977
 pro2e        1     596    1191    1786    2381    2409    2437    3032  702080 1401128
  1925416 1933608 1939068 1960907

 pro2e:    177306 integrals read in    33 records.

 pro2e:         0 integrals 34-ext integrals skipped.
 pro1e        1     596    1191    1786    2381    2409    2437    3032  702080 1401128
  1925416 1933608 1939068 1960907
 pro1e: eref =   -3.498201546360222E+01
 total size of srtscr:                     3  records of                 524288 
 WP =              12582912 Bytes

 new core energy added to the energy(*) list.
 from the hamiltonian repartitioning, eref= -3.498201546360E+01
 putdg        1     596    1191    1786    7930  532218  881743    3032  702080 1401128
  1925416 1933608 1939068 1960907

 putf:       4 buffers of length    6144 written to file 12
 diagonal integral file completed.

 putd34:     3 records of integral w-combinations and
             3 records of integral x-combinations of length= 30006 have been written.
 wstat,xstat=  100.00  100.00
 prep4e:     3 blocks of linear combinations of 4-external integrals processed.
 number of sorted 4-external integrals     132678
 number of original 4-external integrals    81081


 putf34: external integral file complete. nfilw=    31 nfilx=    32 nrecw=     3 nrecx=     3 lbufp= 30006

 putd34:     4 records of integral w-combinations and
             3 records of integral x-combinations of length= 30006 have been written.
 wstat,xstat=   75.00  100.00
 prep3e:     7 blocks of linear combinations of 3-external integrals processed.
 number of sorted 3-external integrals     137592
 number of original 3-external integrals    76167


 putf34: external integral file complete. nfilw=    33 nfilx=    34 nrecw=     4 nrecx=     3 lbufp= 30006
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd     31752         3    157249    157249     31752    196350
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd      6804         4    189001    189001      6804    196350
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd       546         5    195805    195805       546    196350

 putf:      10 buffers of length    6144 written to file 13
 off-diagonal files sort completed.
 !timer: cisrt complete                  cpu_time=     0.148 walltime=     0.589
 executing brd_struct for cisrtinfo
cisrtinfo:
bufszi  6144
 diagfile 4ext:     756 2ext:     378 0ext:      56
 fil4w,fil4x  :   81081 fil3w,fil3x :   76167
 ofdgint  2ext:   31752 1ext:    6804 0ext:     546so0ext:       0so1ext:       0so2ext:       0
buffer minbl4    1131 minbl3    1131 maxbl2    1134nbas:  27   0   0   0   0   0   0   0 maxbuf 30006
 CIUDG version 5.9.7 ( 5-Oct-2004)

 workspace allocation information: lcore= 222822399

 core energy values from the integral file:
 energy( 1)=  3.858870424817E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -9.226324016262E+01, ietype=    6,   fcore energy of type: H1(*)   
 energy( 3)= -3.498201546360E+01, ietype=    5,   fcore energy of type: Vref(*) 

 total core repulsion energy = -8.865655137805E+01
 nmot  =    36 niot  =     7 nfct  =     2 nfvt  =     0
 nrow  =    42 nsym  =     1 ssym  =     1 lenbuf=  1600
 nwalk,xbar:        289       28       92       71       98
 nvalwt,nvalw:      289       28       92       71       98
 ncsft:           64477
 total number of valid internal walks:     289
 nvalz,nvaly,nvalx,nvalw =       28      92      71      98

 cisrt info file parameters:
 file number  12 blocksize   6144
 mxbld   6144
 nd4ext,nd2ext,nd0ext   756   378    56
 n4ext,n3ext,n2ext,n1ext,n0ext,n2int,n1int,n0int    81081    76167    31752     6804      546        0        0        0
 minbl4,minbl3,maxbl2  1131  1131  1134
 maxbuf 30006
 number of external orbitals per symmetry block:  27
 nmsym   1 number of internal orbitals   7
 executing brd_struct for drt
 executing brd_struct for orbinf
 executing brd_struct for momap
 calcthrxt: niot,maxw1=                     7                    95
 block size     0
 pthz,pthy,pthx,pthw:    28    92    71    98 total internal walks:     289
 maxlp3,n2lp,n1lp,n0lp    95     0     0     0
 orbsym(*)= 1 1 1 1 1 1 1
setref: retained number of references =     6
 setref: total/valid number of walks=                    28
                    28
 nmb.of records onel     1
 nmb.of records 2-ext     6
 nmb.of records 1-ext     2
 nmb.of records 0-ext     1
 nmb.of records 2-int     0
 nmb.of records 1-int     0
 nmb.of records 0-int     0
 ---------memory usage in DP -----------------
 < n-ex core usage >
     routines:
    fourex            67263
    threx             60583
    twoex             13947
    onex               7711
    allin              6144
    diagon             7228
               =======
   maximum            67263
 
  __ static summary __ 
   reflst                28
   hrfspc                28
               -------
   static->              28
 
  __ core required  __ 
   totstc                28
   max n-ex           67263
               -------
   totnec->           67291
 
  __ core available __ 
   totspc         222822399
   totnec -           67291
               -------
   totvec->       222755108

 number of external paths / symmetry
 vertex x     351
 vertex w     378
segment: free space=   222755108
 reducing frespc by                   976 to              222754132 
  for index/conft/indsym storage .
 resegmenting ...



                   segmentation summary for type all-internal
 -------------------------------------------------------------------------------
 seg.      no. of|    no. of|  starting|  internal|  starting|  starting|
  no.    internal|        ci|       csf|     walks|      walk|       DRT|
            paths|  elements|    number|     /seg.|    number|    record|
 -------------------------------------------------------------------------------
  Z 1          28|        28|         0|        28|         0|         1|
 -------------------------------------------------------------------------------
  Y 2          92|      2484|        28|        92|        28|         2|
 -------------------------------------------------------------------------------
  X 3          71|     24921|      2512|        71|       120|         3|
 -------------------------------------------------------------------------------
  W 4          98|     37044|     27433|        98|       191|         4|
 -------------------------------------------------------------------------------
max. additional memory requirements:index=           4DP  conft+indsym=         392DP  drtbuffer=         580 DP

dimension of the ci-matrix ->>>     64477

 executing brd_struct for civct
 gentasklist: ntask=                    20
                    TASKLIST
----------------------------------------------------------------------------------------------------
TASK# BRA# KET#  T-TYPE    DESCR.   SEGMENTTYPE    SEGEL              SEGCI          VWALKS   
----------------------------------------------------------------------------------------------------
     1  3   1    24      two-ext xz   2X  3 1      71      28      24921         28      71      28
     2  4   1    25      two-ext wz   2X  4 1      98      28      37044         28      98      28
     3  4   3    26      two-ext wx*  WX  4 3      98      71      37044      24921      98      71
     4  4   3    27      two-ext wx+  WX  4 3      98      71      37044      24921      98      71
     5  2   1    11      one-ext yz   1X  2 1      92      28       2484         28      92      28
     6  3   2    15      1ex3ex yx    3X  3 2      71      92      24921       2484      71      92
     7  4   2    16      1ex3ex yw    3X  4 2      98      92      37044       2484      98      92
     8  1   1     1      allint zz    OX  1 1      28      28         28         28      28      28
     9  2   2     5      0ex2ex yy    OX  2 2      92      92       2484       2484      92      92
    10  3   3     6      0ex2ex xx*   OX  3 3      71      71      24921      24921      71      71
    11  3   3    18      0ex2ex xx+   OX  3 3      71      71      24921      24921      71      71
    12  4   4     7      0ex2ex ww*   OX  4 4      98      98      37044      37044      98      98
    13  4   4    19      0ex2ex ww+   OX  4 4      98      98      37044      37044      98      98
    14  2   2    42      four-ext y   4X  2 2      92      92       2484       2484      92      92
    15  3   3    43      four-ext x   4X  3 3      71      71      24921      24921      71      71
    16  4   4    44      four-ext w   4X  4 4      98      98      37044      37044      98      98
    17  1   1    75      dg-024ext z  OX  1 1      28      28         28         28      28      28
    18  2   2    76      dg-024ext y  OX  2 2      92      92       2484       2484      92      92
    19  3   3    77      dg-024ext x  OX  3 3      71      71      24921      24921      71      71
    20  4   4    78      dg-024ext w  OX  4 4      98      98      37044      37044      98      98
----------------------------------------------------------------------------------------------------
REDTASK #   1 TIME=  19.000 N=  1 (task/type/sgbra)=(   1/24/0) (
REDTASK #   2 TIME=  18.000 N=  1 (task/type/sgbra)=(   2/25/0) (
REDTASK #   3 TIME=  17.000 N=  1 (task/type/sgbra)=(   3/26/1) (
REDTASK #   4 TIME=  16.000 N=  1 (task/type/sgbra)=(   4/27/2) (
REDTASK #   5 TIME=  15.000 N=  1 (task/type/sgbra)=(   5/11/0) (
REDTASK #   6 TIME=  14.000 N=  1 (task/type/sgbra)=(   6/15/0) (
REDTASK #   7 TIME=  13.000 N=  1 (task/type/sgbra)=(   7/16/0) (
REDTASK #   8 TIME=  12.000 N=  1 (task/type/sgbra)=(   8/ 1/0) (
REDTASK #   9 TIME=  11.000 N=  1 (task/type/sgbra)=(   9/ 5/0) (
REDTASK #  10 TIME=  10.000 N=  1 (task/type/sgbra)=(  10/ 6/1) (
REDTASK #  11 TIME=   9.000 N=  1 (task/type/sgbra)=(  11/18/2) (
REDTASK #  12 TIME=   8.000 N=  1 (task/type/sgbra)=(  12/ 7/1) (
REDTASK #  13 TIME=   7.000 N=  1 (task/type/sgbra)=(  13/19/2) (
REDTASK #  14 TIME=   6.000 N=  1 (task/type/sgbra)=(  14/42/1) (
REDTASK #  15 TIME=   5.000 N=  1 (task/type/sgbra)=(  15/43/1) (
REDTASK #  16 TIME=   4.000 N=  1 (task/type/sgbra)=(  16/44/1) (
REDTASK #  17 TIME=   3.000 N=  1 (task/type/sgbra)=(  17/75/1) (
REDTASK #  18 TIME=   2.000 N=  1 (task/type/sgbra)=(  18/76/1) (
REDTASK #  19 TIME=   1.000 N=  1 (task/type/sgbra)=(  19/77/1) (
REDTASK #  20 TIME=   0.000 N=  1 (task/type/sgbra)=(  20/78/1) (
 initializing v-file: 1:                 64477

    ---------trial vector generation----------

    trial vectors will be created by: 

    (ivmode= 3) diagonalizing h in the reference space.                     

      5 vectors will be written to unit 11 beginning with logical record   1

            5 vectors will be created
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         277 2x:           0 4x:           0
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:           0    task #     2:           0    task #     3:           0    task #     4:           0
task #     5:           0    task #     6:           0    task #     7:           0    task #     8:         726
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         251    task #    18:           0    task #    19:           0    task #    20:           0
 reference space has dimension       6
 dsyevx: computed roots 1 to    6(converged:   6)

    root           eigenvalues
    ----           ------------
       1         -94.3985323836
       2         -94.0535597857
       3         -94.0138933507
       4         -93.7714030976
       5         -93.6147885703
       6         -93.5309041066

 strefv generated    5 initial ci vector(s).
    ---------end of vector generation---------

 ufvoutnew: ... writing  recamt=                    28

         vector  1 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                    28

         vector  2 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                    28

         vector  3 from unit 11 written to unit 49 filename cirefv              

 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=    28)...
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:             64477
 number of initial trial vectors:                         5
 number of initial matrix-vector products:                0
 maximum dimension of the subspace vectors:               8
 number of roots to converge:                             3
 number of iterations:                                    1
 residual norm convergence criteria:               0.001000  0.001000  0.001000

          starting bk iteration   1

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2237 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        1945    task #     6:           0    task #     7:           0    task #     8:         726
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2237 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        1945    task #     6:           0    task #     7:           0    task #     8:         726
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2237 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        1945    task #     6:           0    task #     7:           0    task #     8:         726
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2237 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        1945    task #     6:           0    task #     7:           0    task #     8:         726
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2237 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        1945    task #     6:           0    task #     7:           0    task #     8:         726
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -5.74198101
   ht   2    -0.00000000    -5.39700841
   ht   3    -0.00000000     0.00000000    -5.35734197
   ht   4    -0.00000000     0.00000000    -0.00000000    -5.11485172
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000    -4.95823719

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1    1.00000      -1.725485E-14   7.271467E-15  -1.907488E-14  -6.020396E-16
 ref    2   1.730936E-14    1.00000      -1.870907E-13   1.533582E-14  -9.023433E-16
 ref    3   7.235507E-15  -1.871881E-13   -1.00000      -1.965537E-14   1.894209E-15

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1    1.00000        1.00000        1.00000       9.853716E-28   4.764702E-30

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     1.00000000    -0.00000000     0.00000000    -0.00000000    -0.00000000
 ref:   2     0.00000000     1.00000000    -0.00000000     0.00000000    -0.00000000
 ref:   3     0.00000000    -0.00000000    -1.00000000    -0.00000000     0.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3985323836 -1.5099E-14  3.1591E-01  9.7425E-01  1.0000E-03   
 mr-sdci #  1  2    -94.0535597857 -2.4869E-14  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  1  3    -94.0138933507 -1.2434E-14  0.0000E+00  9.9686E-01  1.0000E-03   
 mr-sdci #  1  4    -93.7714030976  1.7764E-15  0.0000E+00  1.0037E+00  1.0000E-04   
 mr-sdci #  1  5    -93.6147885703 -2.1316E-14  0.0000E+00  1.0213E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000996
time for cinew                         0.005898
time for eigenvalue solver             0.000000
time for vector access                 0.000004

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3985323836 -1.5099E-14  3.1591E-01  9.7425E-01  1.0000E-03   
 mr-sdci #  1  2    -94.0535597857 -2.4869E-14  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  1  3    -94.0138933507 -1.2434E-14  0.0000E+00  9.9686E-01  1.0000E-03   
 mr-sdci #  1  4    -93.7714030976  1.7764E-15  0.0000E+00  1.0037E+00  1.0000E-04   
 mr-sdci #  1  5    -93.6147885703 -2.1316E-14  0.0000E+00  1.0213E+00  1.0000E-04   
 
    3 of the   6 expansion vectors are transformed.
    3 of the   5 matrix-vector products are transformed.

    3 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:             64477
 number of initial trial vectors:                         3
 number of initial matrix-vector products:                3
 maximum dimension of the subspace vectors:               8
 number of roots to converge:                             3
 number of iterations:                                   60
 residual norm convergence criteria:               0.001000  0.001000  0.001000

          starting ci iteration   1

 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1    -5.74198101
   ht   2     0.00000000    -5.39700841
   ht   3    -0.00000000    -0.00000000    -5.35734197

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3
 ref    1   -1.00000      -1.627367E-14   7.173407E-15
 ref    2  -1.632179E-14    1.00000      -1.871881E-13
 ref    3  -7.235507E-15  -1.871881E-13   -1.00000    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3
 ref    1    1.00000        1.00000        1.00000    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1    -1.00000000    -0.00000000     0.00000000
 ref:   2    -0.00000000     1.00000000    -0.00000000
 ref:   3    -0.00000000    -0.00000000    -1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3985323836  0.0000E+00  3.1591E-01  9.7425E-01  1.0000E-03   
 mr-sdci #  1  2    -94.0535597857 -8.8818E-16  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  1  3    -94.0138933507 -1.7764E-15  0.0000E+00  9.9686E-01  1.0000E-03   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000709
time for cinew                         0.000999
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration   2

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -5.74198101
   ht   2     0.00000000    -5.39700841
   ht   3    -0.00000000    -0.00000000    -5.35734197
   ht   4     0.31590748    -0.00000055    -0.00584293    -0.34416837

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1   0.961890       8.303030E-07   8.424064E-03   0.273307    
 ref    2  -7.110895E-07    1.00000      -1.416827E-07  -5.309764E-07
 ref    3   7.123604E-03  -1.335579E-07  -0.999958       5.750260E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1   0.925283        1.00000       0.999987       7.472989E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1     0.96188991     0.00000083     0.00842406     0.27330721
 ref:   2    -0.00000071     1.00000000    -0.00000014    -0.00000053
 ref:   3     0.00712360    -0.00000013    -0.99995809     0.00575026

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  2  1    -94.6553953056  2.5686E-01  1.1691E-02  1.8501E-01  1.0000E-03   
 mr-sdci #  2  2    -94.0535597857  4.9560E-13  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  2  3    -94.0139532741  5.9923E-05  0.0000E+00  9.9643E-01  1.0000E-03   
 mr-sdci #  2  4    -91.2172653766 -2.5541E+00  0.0000E+00  1.2100E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000316
time for cinew                         0.004655
time for eigenvalue solver             0.000031
time for vector access                 0.000000

          starting ci iteration   3

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -5.74198101
   ht   2     0.00000000    -5.39700841
   ht   3    -0.00000000    -0.00000000    -5.35734197
   ht   4     0.31590748    -0.00000055    -0.00584293    -0.34416837
   ht   5    -0.02964207    -0.00000041    -0.01182182    -0.00431761    -0.01518393

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.959822      -2.448684E-07  -7.957520E-04  -0.103397      -0.260866    
 ref    2  -2.427896E-07   -1.00000       3.892505E-06  -1.972904E-06   8.154711E-07
 ref    3   1.168346E-03   3.828096E-06   0.999570       2.764275E-02  -9.706839E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.921259        1.00000       0.999141       1.145505E-02   6.814515E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.95982161    -0.00000024    -0.00079575    -0.10339693    -0.26086572
 ref:   2    -0.00000024    -1.00000000     0.00000389    -0.00000197     0.00000082
 ref:   3     0.00116835     0.00000383     0.99957005     0.02764275    -0.00970684

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  3  1    -94.6652633153  9.8680E-03  1.0204E-03  5.3718E-02  1.0000E-03   
 mr-sdci #  3  2    -94.0535597857  1.0804E-11  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  3  3    -94.0158696146  1.9163E-03  0.0000E+00  9.9335E-01  1.0000E-03   
 mr-sdci #  3  4    -91.7756517871  5.5839E-01  0.0000E+00  1.4558E+00  1.0000E-04   
 mr-sdci #  3  5    -91.1996542844 -2.4151E+00  0.0000E+00  1.1847E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000868
time for cinew                         0.005859
time for eigenvalue solver             0.000114
time for vector access                 0.000001

          starting ci iteration   4

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -5.74198101
   ht   2     0.00000000    -5.39700841
   ht   3    -0.00000000    -0.00000000    -5.35734197
   ht   4     0.31590748    -0.00000055    -0.00584293    -0.34416837
   ht   5    -0.02964207    -0.00000041    -0.01182182    -0.00431761    -0.01518393
   ht   6     0.01025825    -0.00000002    -0.00059241     0.00239158     0.00009089    -0.00122399

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.959182       7.917074E-07   8.487621E-03  -0.121627       0.147130      -0.208466    
 ref    2   8.277630E-08   -1.00000       5.361439E-05  -5.319447E-06   5.649053E-07   2.268238E-06
 ref    3  -3.821148E-03   5.283262E-05   0.995578       8.503526E-02  -1.458122E-02  -3.695081E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.920044        1.00000       0.991248       2.202409E-02   2.185995E-02   4.482354E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95918160     0.00000079     0.00848762    -0.12162688     0.14713033    -0.20846626
 ref:   2     0.00000008    -1.00000000     0.00005361    -0.00000532     0.00000056     0.00000227
 ref:   3    -0.00382115     0.00005283     0.99557844     0.08503526    -0.01458122    -0.03695081

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  4  1    -94.6660894366  8.2612E-04  1.0049E-04  1.6896E-02  1.0000E-03   
 mr-sdci #  4  2    -94.0535597858  1.2052E-10  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  4  3    -94.0322127859  1.6343E-02  0.0000E+00  9.5901E-01  1.0000E-03   
 mr-sdci #  4  4    -92.1425071169  3.6686E-01  0.0000E+00  1.1994E+00  1.0000E-04   
 mr-sdci #  4  5    -91.3136379248  1.1398E-01  0.0000E+00  1.4866E+00  1.0000E-04   
 mr-sdci #  4  6    -91.0394269940  2.3829E+00  0.0000E+00  1.3388E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000308
time for cinew                         0.007906
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration   5

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -5.74198101
   ht   2     0.00000000    -5.39700841
   ht   3    -0.00000000    -0.00000000    -5.35734197
   ht   4     0.31590748    -0.00000055    -0.00584293    -0.34416837
   ht   5    -0.02964207    -0.00000041    -0.01182182    -0.00431761    -0.01518393
   ht   6     0.01025825    -0.00000002    -0.00059241     0.00239158     0.00009089    -0.00122399
   ht   7    -0.00485045    -0.00000015    -0.00216262     0.00202895    -0.00038809     0.00001517    -0.00016958

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958798      -1.865814E-02   1.036741E-04   0.123512      -4.113203E-02   9.502585E-02   0.233196    
 ref    2  -2.297721E-07   5.613405E-03   0.999984       1.119197E-05   5.452308E-07   8.603204E-07  -2.567371E-06
 ref    3   5.652947E-03  -0.985588       5.534538E-03  -0.162730      -1.468131E-02  -2.209650E-02   3.698724E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919326       0.971764       0.999999       4.173647E-02   1.907385E-03   9.518168E-03   5.574854E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95879838    -0.01865814     0.00010367     0.12351227    -0.04113203     0.09502585     0.23319624
 ref:   2    -0.00000023     0.00561340     0.99998424     0.00001119     0.00000055     0.00000086    -0.00000257
 ref:   3     0.00565295    -0.98558834     0.00553454    -0.16273041    -0.01468131    -0.02209650     0.03698724

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  5  1    -94.6661859581  9.6521E-05  7.9529E-06  4.9159E-03  1.0000E-03   
 mr-sdci #  5  2    -94.0540489591  4.8917E-04  0.0000E+00  9.0996E-01  1.0000E-03   
 mr-sdci #  5  3    -94.0535597704  2.1347E-02  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  5  4    -92.7694199934  6.2691E-01  0.0000E+00  1.1113E+00  1.0000E-04   
 mr-sdci #  5  5    -91.5595773463  2.4594E-01  0.0000E+00  1.2502E+00  1.0000E-04   
 mr-sdci #  5  6    -91.2671823123  2.2776E-01  0.0000E+00  1.7364E+00  1.0000E-04   
 mr-sdci #  5  7    -90.9413908837  2.2848E+00  0.0000E+00  1.1970E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000328
time for cinew                         0.006499
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   6

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -5.74198101
   ht   2     0.00000000    -5.39700841
   ht   3    -0.00000000    -0.00000000    -5.35734197
   ht   4     0.31590748    -0.00000055    -0.00584293    -0.34416837
   ht   5    -0.02964207    -0.00000041    -0.01182182    -0.00431761    -0.01518393
   ht   6     0.01025825    -0.00000002    -0.00059241     0.00239158     0.00009089    -0.00122399
   ht   7    -0.00485045    -0.00000015    -0.00216262     0.00202895    -0.00038809     0.00001517    -0.00016958
   ht   8     0.00019641    -0.00000002     0.00046952    -0.00005237     0.00009883    -0.00000628     0.00000969    -0.00001034

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958737      -2.519972E-02   1.388295E-06   0.103741       7.923179E-02  -8.119734E-02  -1.065352E-02  -0.237574    
 ref    2  -3.052585E-07   1.134221E-04    1.00000       1.875897E-05  -3.160586E-06   1.380832E-06  -4.028979E-06   1.890852E-06
 ref    3   6.621153E-03  -0.959920       1.143600E-04  -0.268506       3.393674E-02  -2.163653E-02   6.404716E-02  -2.630771E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.922082        1.00000       8.285778E-02   7.429379E-03   7.061148E-03   4.215536E-03   5.713347E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873728    -0.02519972     0.00000139     0.10374127     0.07923179    -0.08119734    -0.01065352    -0.23757393
 ref:   2    -0.00000031     0.00011342     0.99999999     0.00001876    -0.00000316     0.00000138    -0.00000403     0.00000189
 ref:   3     0.00662115    -0.95992012     0.00011436    -0.26850610     0.03393674    -0.02163653     0.06404716    -0.02630771

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  6  1    -94.6661940993  8.1412E-06  7.2805E-07  1.4348E-03  1.0000E-03   
 mr-sdci #  6  2    -94.1151365456  6.1088E-02  0.0000E+00  7.7027E-01  1.0000E-03   
 mr-sdci #  6  3    -94.0535597853  1.4906E-08  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  6  4    -92.9724489347  2.0303E-01  0.0000E+00  1.1107E+00  1.0000E-04   
 mr-sdci #  6  5    -91.8558139384  2.9624E-01  0.0000E+00  1.2513E+00  1.0000E-04   
 mr-sdci #  6  6    -91.5225409710  2.5536E-01  0.0000E+00  1.2015E+00  1.0000E-04   
 mr-sdci #  6  7    -90.9772779052  3.5887E-02  0.0000E+00  1.7491E+00  1.0000E-04   
 mr-sdci #  6  8    -90.9402783879  2.2837E+00  0.0000E+00  1.1753E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000004
time for cinew                         0.009628
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration   7

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964272
   ht   2     0.00000000    -5.45858517
   ht   3     0.00000000     0.00000000    -5.39700841
   ht   4     0.00000000    -0.00000000     0.00000000    -4.31589756
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000    -3.19926256
   ht   6    -0.00031528     0.00032981    -0.00000000    -0.00047621     0.00030606    -0.00000131

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958742      -2.309527E-02   1.052478E-06  -6.377640E-02  -0.113599       3.139419E-02
 ref    2   3.276967E-07   9.936487E-05    1.00000      -2.633503E-05   1.699835E-06  -7.421383E-07
 ref    3  -6.862374E-03  -0.934070       1.023458E-04   0.349276       7.197387E-03   1.586257E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919233       0.873020        1.00000       0.126061       1.295657E-02   1.237217E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95874197    -0.02309527     0.00000105    -0.06377640    -0.11359914     0.03139419
 ref:   2     0.00000033     0.00009936     0.99999999    -0.00002634     0.00000170    -0.00000074
 ref:   3    -0.00686237    -0.93406998     0.00010235     0.34927590     0.00719739     0.01586257

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  7  1    -94.6661948535  7.5419E-07  0.0000E+00  5.4190E-04  1.0000E-03   
 mr-sdci #  7  2    -94.1564883270  4.1352E-02  1.5804E-01  6.5973E-01  1.0000E-03   
 mr-sdci #  7  3    -94.0535597853  2.2022E-11  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  7  4    -93.1491192133  1.7667E-01  0.0000E+00  1.2451E+00  1.0000E-04   
 mr-sdci #  7  5    -92.4261603682  5.7035E-01  0.0000E+00  1.0931E+00  1.0000E-04   
 mr-sdci #  7  6    -91.5155980552 -6.9429E-03  0.0000E+00  1.3287E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000081
time for cinew                         0.009837
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   8

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964272
   ht   2     0.00000000    -5.45858517
   ht   3     0.00000000     0.00000000    -5.39700841
   ht   4     0.00000000    -0.00000000     0.00000000    -4.31589756
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000    -3.19926256
   ht   6    -0.00031528     0.00032981    -0.00000000    -0.00047621     0.00030606    -0.00000131
   ht   7     0.00300738     0.04274190     0.00003138    -0.13061010    -0.14924755    -0.00005683    -0.22577641

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958742       1.320679E-02   1.678992E-06   7.284461E-02  -0.103495       7.913286E-03  -4.740115E-02
 ref    2  -3.416167E-07  -1.164463E-04    1.00000       4.540764E-06   1.399204E-05  -1.514592E-05  -2.754991E-05
 ref    3   6.969172E-03   0.934998       1.200063E-04  -0.192535      -0.102435       0.139336       0.231759    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919236       0.874396        1.00000       4.237595E-02   2.120411E-02   1.947720E-02   5.595912E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95874247     0.01320679     0.00000168     0.07284461    -0.10349524     0.00791329    -0.04740115
 ref:   2    -0.00000034    -0.00011645     0.99999999     0.00000454     0.00001399    -0.00001515    -0.00002755
 ref:   3     0.00696917     0.93499842     0.00012001    -0.19253470    -0.10243461     0.13933621     0.23175904

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  8  1    -94.6661948728  1.9265E-08  0.0000E+00  4.7668E-04  1.0000E-03   
 mr-sdci #  8  2    -94.2911760192  1.3469E-01  1.7735E-02  2.1854E-01  1.0000E-03   
 mr-sdci #  8  3    -94.0535597854  6.5751E-11  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  8  4    -93.2630546754  1.1394E-01  0.0000E+00  1.0404E+00  1.0000E-04   
 mr-sdci #  8  5    -92.4788319785  5.2672E-02  0.0000E+00  1.1112E+00  1.0000E-04   
 mr-sdci #  8  6    -91.5751788122  5.9581E-02  0.0000E+00  1.4266E+00  1.0000E-04   
 mr-sdci #  8  7    -91.2371173176  2.5984E-01  0.0000E+00  1.3442E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001683
time for cinew                         0.005352
time for eigenvalue solver             0.000155
time for vector access                 0.000001

          starting ci iteration   9

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964272
   ht   2     0.00000000    -5.45858517
   ht   3     0.00000000     0.00000000    -5.39700841
   ht   4     0.00000000    -0.00000000     0.00000000    -4.31589756
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000    -3.19926256
   ht   6    -0.00031528     0.00032981    -0.00000000    -0.00047621     0.00030606    -0.00000131
   ht   7     0.00300738     0.04274190     0.00003138    -0.13061010    -0.14924755    -0.00005683    -0.22577641
   ht   8    -0.85391580    -0.15668039     0.00003061     0.11065570     0.00741286    -0.00006937    -0.01318365    -0.16254880

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958738      -6.780219E-03   4.573258E-06  -8.532677E-02  -0.103332      -1.574650E-02   3.346766E-02   5.092524E-02
 ref    2   3.330159E-07   1.329466E-04    1.00000       2.374097E-05   1.341449E-05   1.143379E-05   3.543231E-05  -1.844381E-05
 ref    3  -6.939001E-03  -0.941697       1.333522E-04   0.152573      -0.100755       3.139550E-02  -0.267158       5.691353E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919227       0.886838        1.00000       3.055925E-02   2.082906E-02   1.233629E-03   7.249351E-02   5.832530E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873822    -0.00678022     0.00000457    -0.08532677    -0.10333196    -0.01574650     0.03346766     0.05092524
 ref:   2     0.00000033     0.00013295     0.99999999     0.00002374     0.00001341     0.00001143     0.00003543    -0.00001844
 ref:   3    -0.00693900    -0.94169651     0.00013335     0.15257324    -0.10075496     0.03139550    -0.26715804     0.05691353

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1    -94.6661948815  8.7212E-09  0.0000E+00  4.4368E-04  1.0000E-03   
 mr-sdci #  9  2    -94.3032890850  1.2113E-02 -5.8246E-03  1.2556E-01  1.0000E-03   
 mr-sdci #  9  3    -94.0535597869  1.5296E-09  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci #  9  4    -93.3951381226  1.3208E-01  0.0000E+00  8.4713E-01  1.0000E-04   
 mr-sdci #  9  5    -92.4790105699  1.7859E-04  0.0000E+00  1.1031E+00  1.0000E-04   
 mr-sdci #  9  6    -91.8767230944  3.0154E-01  0.0000E+00  1.5561E+00  1.0000E-04   
 mr-sdci #  9  7    -91.2443960071  7.2787E-03  0.0000E+00  1.4229E+00  1.0000E-04   
 mr-sdci #  9  8    -91.1693815501  2.2910E-01  0.0000E+00  1.4620E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000250
time for cinew                         0.008925
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  10

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964350
   ht   2     0.00000000    -5.64673771
   ht   3    -0.00000000    -0.00000000    -5.39700841
   ht   4    -0.00000000     0.00000000     0.00000000    -4.73858674
   ht   5     0.00000000    -0.00000000     0.00000000    -0.00000000    -3.82245919
   ht   6    -2.55273704     0.94944583    -0.00007119     0.03468462     0.25084670    -1.28128974

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958736       1.050386E-02   7.476088E-06   8.009525E-02  -0.106239       0.205707    
 ref    2  -3.327355E-07  -1.299451E-04    1.00000      -2.165147E-05   1.408216E-05  -3.501131E-05
 ref    3   6.941948E-03   0.944540       1.323700E-04  -0.154756      -0.102249       0.127902    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919224       0.892266        1.00000       3.036479E-02   2.174160E-02   5.867451E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873634     0.01050386     0.00000748     0.08009525    -0.10623910     0.20570743
 ref:   2    -0.00000033    -0.00012995     0.99999999    -0.00002165     0.00001408    -0.00003501
 ref:   3     0.00694195     0.94453991     0.00013237    -0.15475640    -0.10224897     0.12790217

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 10  1    -94.6661948818  3.1374E-10  0.0000E+00  4.4263E-04  1.0000E-03   
 mr-sdci # 10  2    -94.3044346637  1.1456E-03  2.7689E-03  1.0433E-01  1.0000E-03   
 mr-sdci # 10  3    -94.0535597876  6.7045E-10  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci # 10  4    -93.3966714946  1.5334E-03  0.0000E+00  8.3788E-01  1.0000E-04   
 mr-sdci # 10  5    -92.4793006014  2.9003E-04  0.0000E+00  1.1014E+00  1.0000E-04   
 mr-sdci # 10  6    -90.9696358215 -9.0709E-01  0.0000E+00  1.2926E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000197
time for cinew                         0.008401
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  11

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964350
   ht   2     0.00000000    -5.64673771
   ht   3    -0.00000000    -0.00000000    -5.39700841
   ht   4    -0.00000000     0.00000000     0.00000000    -4.73858674
   ht   5     0.00000000    -0.00000000     0.00000000    -0.00000000    -3.82245919
   ht   6    -2.55273704     0.94944583    -0.00007119     0.03468462     0.25084670    -1.28128974
   ht   7     0.94812151    -0.54902966     0.00003119    -0.06111893    -0.07554185     0.51153131    -0.21408326

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958738       1.330020E-02  -4.310894E-06   3.749153E-02   0.138096       6.631376E-03  -0.227354    
 ref    2  -3.297188E-07  -1.469569E-04   -1.00000      -1.003668E-04   1.122200E-05   6.409322E-05   1.510004E-05
 ref    3   6.939069E-03   0.943009      -1.515679E-04  -8.460442E-02   1.626865E-02  -0.239841      -5.108420E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919226       0.889442        1.00000       8.563533E-03   1.933521E-02   5.756745E-02   5.429942E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873762     0.01330020    -0.00000431     0.03749153     0.13809612     0.00663138    -0.22735397
 ref:   2    -0.00000033    -0.00014696    -0.99999998    -0.00010037     0.00001122     0.00006409     0.00001510
 ref:   3     0.00693907     0.94300868    -0.00015157    -0.08460442     0.01626865    -0.23984050    -0.05108420

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 11  1    -94.6661948824  6.1244E-10  0.0000E+00  4.4897E-04  1.0000E-03   
 mr-sdci # 11  2    -94.3079589275  3.5243E-03  1.4206E-03  5.8613E-02  1.0000E-03   
 mr-sdci # 11  3    -94.0535597949  7.3090E-09  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci # 11  4    -93.5714957919  1.7482E-01  0.0000E+00  7.6369E-01  1.0000E-04   
 mr-sdci # 11  5    -92.5876160341  1.0832E-01  0.0000E+00  9.8765E-01  1.0000E-04   
 mr-sdci # 11  6    -91.4798367473  5.1020E-01  0.0000E+00  1.4022E+00  1.0000E-04   
 mr-sdci # 11  7    -90.8911985730 -3.5320E-01  0.0000E+00  1.2235E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000422
time for cinew                         0.008124
time for eigenvalue solver             0.000192
time for vector access                 0.000002

          starting ci iteration  12

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964350
   ht   2     0.00000000    -5.64673771
   ht   3    -0.00000000    -0.00000000    -5.39700841
   ht   4    -0.00000000     0.00000000     0.00000000    -4.73858674
   ht   5     0.00000000    -0.00000000     0.00000000    -0.00000000    -3.82245919
   ht   6    -2.55273704     0.94944583    -0.00007119     0.03468462     0.25084670    -1.28128974
   ht   7     0.94812151    -0.54902966     0.00003119    -0.06111893    -0.07554185     0.51153131    -0.21408326
   ht   8     0.14147663     0.08948067     0.00000423    -0.01793379    -0.01099275     0.04714263    -0.01435129    -0.00732527

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958737       1.276061E-02  -9.424009E-06   3.629499E-02   0.128511       5.962161E-02  -9.514516E-02   0.204879    
 ref    2  -3.101760E-07  -1.704202E-04   -1.00000      -2.437434E-04   5.684840E-05  -1.221045E-04   3.523991E-06  -4.434910E-06
 ref    3   6.928504E-03   0.945389      -1.583461E-04  -8.261557E-02   7.300606E-03   0.104445       0.212139       0.117082    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919224       0.893923        1.00000       8.142719E-03   1.656848E-02   1.446352E-02   5.405561E-02   5.568360E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873682     0.01276061    -0.00000942     0.03629499     0.12851139     0.05962161    -0.09514516     0.20487877
 ref:   2    -0.00000031    -0.00017042    -0.99999995    -0.00024374     0.00005685    -0.00012210     0.00000352    -0.00000443
 ref:   3     0.00692850     0.94538904    -0.00015835    -0.08261557     0.00730061     0.10444506     0.21213912     0.11708239

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 12  1    -94.6661948875  5.0693E-09  0.0000E+00  4.3310E-04  1.0000E-03   
 mr-sdci # 12  2    -94.3090580862  1.0992E-03  1.9619E-04  2.2526E-02  1.0000E-03   
 mr-sdci # 12  3    -94.0535598380  4.3055E-08  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci # 12  4    -93.6448130546  7.3317E-02  0.0000E+00  6.8395E-01  1.0000E-04   
 mr-sdci # 12  5    -92.6428972830  5.5281E-02  0.0000E+00  8.9769E-01  1.0000E-04   
 mr-sdci # 12  6    -91.9498213726  4.6998E-01  0.0000E+00  1.4696E+00  1.0000E-04   
 mr-sdci # 12  7    -91.1756674558  2.8447E-01  0.0000E+00  1.5069E+00  1.0000E-04   
 mr-sdci # 12  8    -90.8416008483 -3.2778E-01  0.0000E+00  1.2696E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000285
time for cinew                         0.009462
time for eigenvalue solver             0.000188
time for vector access                 0.000002

          starting ci iteration  13

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964351
   ht   2    -0.00000000    -5.65250671
   ht   3    -0.00000000    -0.00000000    -5.39700846
   ht   4    -0.00000000    -0.00000000    -0.00000000    -4.98826168
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000    -3.98634590
   ht   6     0.00761506     0.00037457    -0.00001778     0.01886861    -0.01036236    -0.00047569

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958736      -1.265455E-02   2.278355E-05  -2.578639E-02  -0.114713      -6.571373E-02
 ref    2   2.097187E-07   2.051860E-04   0.999999       1.163999E-03  -2.833070E-04   2.559321E-04
 ref    3  -6.926675E-03  -0.945231       1.863954E-04   5.067587E-02   4.333467E-02  -8.282300E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919223       0.893622       0.999998       3.234337E-03   1.503714E-02   1.117801E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873638    -0.01265455     0.00002278    -0.02578639    -0.11471342    -0.06571373
 ref:   2     0.00000021     0.00020519     0.99999922     0.00116400    -0.00028331     0.00025593
 ref:   3    -0.00692668    -0.94523116     0.00018640     0.05067587     0.04333467    -0.08282300

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 13  1    -94.6661949007  1.3274E-08  0.0000E+00  4.0711E-04  1.0000E-03   
 mr-sdci # 13  2    -94.3092643737  2.0629E-04  6.9683E-05  1.2280E-02  1.0000E-03   
 mr-sdci # 13  3    -94.0535604205  5.8251E-07  0.0000E+00  9.9729E-01  1.0000E-03   
 mr-sdci # 13  4    -93.7823210138  1.3751E-01  0.0000E+00  5.7618E-01  1.0000E-04   
 mr-sdci # 13  5    -92.8386477062  1.9575E-01  0.0000E+00  8.6452E-01  1.0000E-04   
 mr-sdci # 13  6    -91.8290998134 -1.2072E-01  0.0000E+00  1.4148E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000431
time for cinew                         0.006991
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  14

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964351
   ht   2    -0.00000000    -5.65250671
   ht   3    -0.00000000    -0.00000000    -5.39700846
   ht   4    -0.00000000    -0.00000000    -0.00000000    -4.98826168
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000    -3.98634590
   ht   6     0.00761506     0.00037457    -0.00001778     0.01886861    -0.01036236    -0.00047569
   ht   7    -0.04122763     0.03749664    -0.00000593    -0.00196203     0.00586978     0.00009054    -0.00070070

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958737      -1.273902E-02   1.275983E-05  -1.733282E-02  -0.105897      -7.645888E-02   3.057860E-02
 ref    2   1.140692E-07   2.347492E-04   0.999996       2.666259E-03  -7.377037E-04   5.416246E-04   2.044737E-04
 ref    3  -6.920414E-03  -0.945882       4.002232E-05   6.451059E-02   1.197100E-02   2.073582E-02   0.155882    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919225       0.894855       0.999992       4.469151E-03   1.135807E-02   6.276228E-03   2.523426E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873716    -0.01273902     0.00001276    -0.01733282    -0.10589721    -0.07645888     0.03057860
 ref:   2     0.00000011     0.00023475     0.99999597     0.00266626    -0.00073770     0.00054162     0.00020447
 ref:   3    -0.00692041    -0.94588185     0.00004002     0.06451059     0.01197100     0.02073582     0.15588191

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 14  1    -94.6661949046  3.8916E-09  0.0000E+00  3.8840E-04  1.0000E-03   
 mr-sdci # 14  2    -94.3093100342  4.5660E-05 -1.9243E-05  8.2597E-03  1.0000E-03   
 mr-sdci # 14  3    -94.0535624818  2.0613E-06  0.0000E+00  9.9728E-01  1.0000E-03   
 mr-sdci # 14  4    -93.8374151920  5.5094E-02  0.0000E+00  5.1157E-01  1.0000E-04   
 mr-sdci # 14  5    -92.9869495564  1.4830E-01  0.0000E+00  8.9785E-01  1.0000E-04   
 mr-sdci # 14  6    -92.1247283750  2.9563E-01  0.0000E+00  1.3704E+00  1.0000E-04   
 mr-sdci # 14  7    -91.3738375883  1.9817E-01  0.0000E+00  1.5680E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.007603
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  15

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964351
   ht   2    -0.00000000    -5.65250671
   ht   3    -0.00000000    -0.00000000    -5.39700846
   ht   4    -0.00000000    -0.00000000    -0.00000000    -4.98826168
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000    -3.98634590
   ht   6     0.00761506     0.00037457    -0.00001778     0.01886861    -0.01036236    -0.00047569
   ht   7    -0.04122763     0.03749664    -0.00000593    -0.00196203     0.00586978     0.00009054    -0.00070070
   ht   8    -0.17191153     0.04074421     0.00000293    -0.00843439     0.01896480     0.00032049    -0.00152125    -0.00542121

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958736       1.251478E-02  -7.019517E-05  -2.149745E-02  -0.103624      -8.638902E-02  -5.391571E-02  -0.209083    
 ref    2   1.269519E-07  -2.301330E-04   0.999996       2.572274E-03  -7.264927E-04   5.135076E-04  -2.520685E-04  -3.834772E-04
 ref    3  -6.920467E-03   0.945868       3.360257E-05   6.398071E-02   1.224405E-02   1.950237E-02  -0.157463      -5.259402E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919222       0.894822       0.999992       4.562288E-03   1.088837E-02   7.843669E-03   2.770160E-02   4.374350E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873572     0.01251478    -0.00007020    -0.02149745    -0.10362396    -0.08638902    -0.05391571    -0.20908297
 ref:   2     0.00000013    -0.00023013     0.99999616     0.00257227    -0.00072649     0.00051351    -0.00025207    -0.00038348
 ref:   3    -0.00692047     0.94586754     0.00003360     0.06398071     0.01224405     0.01950237    -0.15746312    -0.00525940

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 15  1    -94.6661949048  1.4956E-10  0.0000E+00  3.8617E-04  1.0000E-03   
 mr-sdci # 15  2    -94.3093132810  3.2468E-06  1.1931E-05  7.4462E-03  1.0000E-03   
 mr-sdci # 15  3    -94.0535628805  3.9865E-07  0.0000E+00  9.9728E-01  1.0000E-03   
 mr-sdci # 15  4    -93.8383874940  9.7230E-04  0.0000E+00  5.0603E-01  1.0000E-04   
 mr-sdci # 15  5    -92.9871387645  1.8921E-04  0.0000E+00  8.9660E-01  1.0000E-04   
 mr-sdci # 15  6    -92.1266460711  1.9177E-03  0.0000E+00  1.3698E+00  1.0000E-04   
 mr-sdci # 15  7    -91.3759289025  2.0913E-03  0.0000E+00  1.5495E+00  1.0000E-04   
 mr-sdci # 15  8    -91.2044592223  3.6286E-01  0.0000E+00  1.4080E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000361
time for cinew                         0.013226
time for eigenvalue solver             0.000187
time for vector access                 0.000000

          starting ci iteration  16

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964353
   ht   2     0.00000000    -5.65276190
   ht   3    -0.00000000    -0.00000000    -5.39701150
   ht   4    -0.00000000     0.00000000     0.00000000    -5.18183612
   ht   5     0.00000000    -0.00000000    -0.00000000     0.00000000    -4.33058739
   ht   6    -0.08699336    -0.00879703     0.00003955    -0.00951305     0.01006593    -0.00136398

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958738       1.269955E-02   1.845518E-04   9.742757E-03   0.129088       0.171673    
 ref    2  -1.116668E-07  -2.345890E-04  -0.999995      -2.989471E-03   9.008818E-04   5.529392E-04
 ref    3   6.920736E-03   0.945840      -4.439676E-05  -6.090879E-02  -2.025136E-02  -4.780541E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919225       0.894774       0.999990       3.813738E-03   1.707456E-02   3.175724E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873750     0.01269955     0.00018455     0.00974276     0.12908770     0.17167288
 ref:   2    -0.00000011    -0.00023459    -0.99999473    -0.00298947     0.00090088     0.00055294
 ref:   3     0.00692074     0.94583955    -0.00004440    -0.06090879    -0.02025136    -0.04780541

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 16  1    -94.6661949051  2.6656E-10  0.0000E+00  3.8670E-04  1.0000E-03   
 mr-sdci # 16  2    -94.3093157874  2.5064E-06 -1.4423E-05  7.0420E-03  1.0000E-03   
 mr-sdci # 16  3    -94.0535637627  8.8225E-07  0.0000E+00  9.9728E-01  1.0000E-03   
 mr-sdci # 16  4    -93.8456954251  7.3079E-03  0.0000E+00  4.9362E-01  1.0000E-04   
 mr-sdci # 16  5    -93.0219598654  3.4821E-02  0.0000E+00  8.4570E-01  1.0000E-04   
 mr-sdci # 16  6    -91.2342938759 -8.9235E-01  0.0000E+00  1.4655E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000231
time for cinew                         0.004804
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  17

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964353
   ht   2     0.00000000    -5.65276190
   ht   3    -0.00000000    -0.00000000    -5.39701150
   ht   4    -0.00000000     0.00000000     0.00000000    -5.18183612
   ht   5     0.00000000    -0.00000000    -0.00000000     0.00000000    -4.33058739
   ht   6    -0.08699336    -0.00879703     0.00003955    -0.00951305     0.01006593    -0.00136398
   ht   7    -0.14807183    -0.03625283     0.00005385    -0.01077359     0.01343609    -0.00231136    -0.00403362

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735       1.246848E-02   2.962200E-05   1.898952E-02   8.827121E-02   0.237562      -3.382038E-02
 ref    2  -2.749704E-08  -2.581606E-04  -0.999984      -5.366831E-03   1.407309E-03   1.792200E-04   6.699157E-04
 ref    3   6.924681E-03   0.945407      -3.071515E-04  -2.940173E-02  -8.102046E-02   8.866988E-02  -0.207412    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.893950       0.999968       1.253866E-03   1.435810E-02   6.429799E-02   4.416417E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873523     0.01246848     0.00002962     0.01898952     0.08827121     0.23756180    -0.03382038
 ref:   2    -0.00000003    -0.00025816    -0.99998413    -0.00536683     0.00140731     0.00017922     0.00066992
 ref:   3     0.00692468     0.94540718    -0.00030715    -0.02940173    -0.08102046     0.08866988    -0.20741241

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 17  1    -94.6661949063  1.2963E-09  0.0000E+00  3.7568E-04  1.0000E-03   
 mr-sdci # 17  2    -94.3093286935  1.2906E-05  8.1255E-06  4.0139E-03  1.0000E-03   
 mr-sdci # 17  3    -94.0535687791  5.0164E-06  0.0000E+00  9.9725E-01  1.0000E-03   
 mr-sdci # 17  4    -93.8819507514  3.6255E-02  0.0000E+00  3.8799E-01  1.0000E-04   
 mr-sdci # 17  5    -93.1467726353  1.2481E-01  0.0000E+00  6.6239E-01  1.0000E-04   
 mr-sdci # 17  6    -91.3331023301  9.8808E-02  0.0000E+00  1.3533E+00  1.0000E-04   
 mr-sdci # 17  7    -91.0406037232 -3.3533E-01  0.0000E+00  1.5204E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000204
time for cinew                         0.006754
time for eigenvalue solver             0.000155
time for vector access                 0.000000

          starting ci iteration  18

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964353
   ht   2     0.00000000    -5.65276190
   ht   3    -0.00000000    -0.00000000    -5.39701150
   ht   4    -0.00000000     0.00000000     0.00000000    -5.18183612
   ht   5     0.00000000    -0.00000000    -0.00000000     0.00000000    -4.33058739
   ht   6    -0.08699336    -0.00879703     0.00003955    -0.00951305     0.01006593    -0.00136398
   ht   7    -0.14807183    -0.03625283     0.00005385    -0.01077359     0.01343609    -0.00231136    -0.00403362
   ht   8    -0.01265660    -0.01567404    -0.00000643    -0.00040261     0.00224396    -0.00021956    -0.00043166    -0.00009223

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958736      -1.235578E-02  -2.567364E-04  -2.320460E-02  -6.003346E-02  -9.142111E-02   0.243716      -2.168000E-03
 ref    2  -6.811651E-08   3.018041E-04  -0.999933       1.056518E-02  -3.700670E-03   2.784595E-03   2.645790E-05  -1.447644E-04
 ref    3   6.925228E-03  -0.945611       1.664875E-04   3.767239E-02   4.392880E-02   0.100508       9.335164E-02   0.234801    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919222       0.894334       0.999867       2.069286E-03   5.547450E-03   1.846748E-02   6.811217E-02   5.513614E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873551    -0.01235578    -0.00025674    -0.02320460    -0.06003346    -0.09142111     0.24371631    -0.00216800
 ref:   2    -0.00000007     0.00030180    -0.99993326     0.01056518    -0.00370067     0.00278459     0.00002646    -0.00014476
 ref:   3     0.00692523    -0.94561144     0.00016649     0.03767239     0.04392880     0.10050822     0.09335164     0.23480080

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 18  1    -94.6661949064  4.7290E-11  0.0000E+00  3.7636E-04  1.0000E-03   
 mr-sdci # 18  2    -94.3093354235  6.7300E-06  1.4285E-06  1.8565E-03  1.0000E-03   
 mr-sdci # 18  3    -94.0535988006  3.0021E-05  0.0000E+00  9.9714E-01  1.0000E-03   
 mr-sdci # 18  4    -93.9021553455  2.0205E-02  0.0000E+00  3.4659E-01  1.0000E-04   
 mr-sdci # 18  5    -93.2637088822  1.1694E-01  0.0000E+00  7.0614E-01  1.0000E-04   
 mr-sdci # 18  6    -92.2419385730  9.0884E-01  0.0000E+00  1.3117E+00  1.0000E-04   
 mr-sdci # 18  7    -91.3296226891  2.8902E-01  0.0000E+00  1.3491E+00  1.0000E-04   
 mr-sdci # 18  8    -90.9465688106 -2.5789E-01  0.0000E+00  1.4863E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001189
time for cinew                         0.007981
time for eigenvalue solver             0.000154
time for vector access                 0.000002

          starting ci iteration  19

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964353
   ht   2    -0.00000000    -5.65278405
   ht   3    -0.00000000     0.00000000    -5.39704742
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.24560397
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000    -4.60715750
   ht   6    -0.00066259    -0.00010807     0.00001861     0.00050963    -0.00091002    -0.00000285

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735      -1.237089E-02   7.719861E-05   2.041455E-02   5.926587E-02   1.469365E-02
 ref    2   3.344615E-07   3.568476E-04   0.999718      -1.962736E-02   1.112331E-02  -6.902772E-03
 ref    3  -6.924772E-03  -0.945570       3.531138E-04  -3.231877E-02  -5.777740E-02   3.298694E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.894255       0.999435       1.846490E-03   6.974399E-03   1.351690E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873536    -0.01237089     0.00007720     0.02041455     0.05926587     0.01469365
 ref:   2     0.00000033     0.00035685     0.99971756    -0.01962736     0.01112331    -0.00690277
 ref:   3    -0.00692477    -0.94556975     0.00035311    -0.03231877    -0.05777740     0.03298694

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 19  1    -94.6661949066  2.2050E-10  0.0000E+00  3.7456E-04  1.0000E-03   
 mr-sdci # 19  2    -94.3093367130  1.2896E-06  0.0000E+00  8.1567E-04  1.0000E-03   
 mr-sdci # 19  3    -94.0538083787  2.0958E-04  3.6868E-01  9.9650E-01  1.0000E-03   
 mr-sdci # 19  4    -93.9095868071  7.4315E-03  0.0000E+00  3.3843E-01  1.0000E-04   
 mr-sdci # 19  5    -93.3881306072  1.2442E-01  0.0000E+00  5.8950E-01  1.0000E-04   
 mr-sdci # 19  6    -91.9064824397 -3.3546E-01  0.0000E+00  1.4919E+00  1.0000E-04   
 
 root number  3 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000138
time for cinew                         0.014882
time for eigenvalue solver             0.000135
time for vector access                 0.000000

          starting ci iteration  20

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964353
   ht   2    -0.00000000    -5.65278405
   ht   3    -0.00000000     0.00000000    -5.39704742
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.24560397
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000    -4.60715750
   ht   6    -0.00066259    -0.00010807     0.00001861     0.00050963    -0.00091002    -0.00000285
   ht   7     0.10408084    -0.33358793     0.36313912    -0.02425050     0.01360188     0.00002321    -0.57736774

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958735      -1.396335E-03  -1.237232E-02   2.039803E-02  -5.928813E-02   1.473099E-02   3.795918E-03
 ref    2   9.231034E-08  -0.951985      -6.951951E-04  -4.425329E-03  -5.064938E-03  -1.066814E-02  -0.305871    
 ref    3  -6.924770E-03   5.011930E-03  -0.945565      -3.226858E-02   5.785201E-02   3.284666E-02  -1.437043E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.906302       0.894247       1.476925E-03   6.887591E-03   1.409714E-03   9.377792E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873536    -0.00139633    -0.01237232     0.02039803    -0.05928813     0.01473099     0.00379592
 ref:   2     0.00000009    -0.95198499    -0.00069520    -0.00442533    -0.00506494    -0.01066814    -0.30587089
 ref:   3    -0.00692477     0.00501193    -0.94556515    -0.03226858     0.05785201     0.03284666    -0.01437043

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 20  1    -94.6661949066  5.5955E-14  0.0000E+00  3.7456E-04  1.0000E-03   
 mr-sdci # 20  2    -94.3290675798  1.9731E-02  1.8976E-02  2.2466E-01  1.0000E-03   
 mr-sdci # 20  3    -94.3093366908  2.5553E-01  0.0000E+00  8.7025E-04  1.0000E-03   
 mr-sdci # 20  4    -93.9096348636  4.8056E-05  0.0000E+00  3.3756E-01  1.0000E-04   
 mr-sdci # 20  5    -93.3881932600  6.2653E-05  0.0000E+00  5.8911E-01  1.0000E-04   
 mr-sdci # 20  6    -91.9065346911  5.2251E-05  0.0000E+00  1.4919E+00  1.0000E-04   
 mr-sdci # 20  7    -91.3859361345  5.6313E-02  0.0000E+00  1.5069E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000023
time for cinew                         0.013921
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  21

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964353
   ht   2    -0.00000000    -5.65278405
   ht   3    -0.00000000     0.00000000    -5.39704742
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.24560397
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000    -4.60715750
   ht   6    -0.00066259    -0.00010807     0.00001861     0.00050963    -0.00091002    -0.00000285
   ht   7     0.10408084    -0.33358793     0.36313912    -0.02425050     0.01360188     0.00002321    -0.57736774
   ht   8     0.20422255    -0.18707093    -0.01364099    -0.01190362    -0.00395270     0.00002202    -0.04372186    -0.04087918

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735       1.775902E-03  -1.236780E-02  -2.085514E-02  -5.788500E-02  -2.753485E-02   3.676854E-02  -1.614957E-02
 ref    2   1.062780E-07  -0.950416       6.759979E-05  -6.425041E-04   5.634256E-04   6.692987E-03  -7.541793E-02  -0.301609    
 ref    3  -6.924768E-03  -1.809371E-03  -0.945575       3.317270E-02   5.531975E-02  -6.120324E-03  -8.845749E-02   2.808118E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.903297       0.894266       1.535778E-03   6.411266E-03   8.404225E-04   1.486452E-02   9.201732E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873536     0.00177590    -0.01236780    -0.02085514    -0.05788500    -0.02753485     0.03676854    -0.01614957
 ref:   2     0.00000011    -0.95041616     0.00006760    -0.00064250     0.00056343     0.00669299    -0.07541793    -0.30160896
 ref:   3    -0.00692477    -0.00180937    -0.94557535     0.03317270     0.05531975    -0.00612032    -0.08845749     0.02808118

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 21  1    -94.6661949066  2.6645E-15  0.0000E+00  3.7456E-04  1.0000E-03   
 mr-sdci # 21  2    -94.3439016427  1.4834E-02  5.1778E-04  6.1645E-02  1.0000E-03   
 mr-sdci # 21  3    -94.3093367202  2.9388E-08  0.0000E+00  8.0341E-04  1.0000E-03   
 mr-sdci # 21  4    -93.9099507578  3.1589E-04  0.0000E+00  3.3661E-01  1.0000E-04   
 mr-sdci # 21  5    -93.3898880185  1.6948E-03  0.0000E+00  5.9041E-01  1.0000E-04   
 mr-sdci # 21  6    -91.9352864399  2.8752E-02  0.0000E+00  1.4664E+00  1.0000E-04   
 mr-sdci # 21  7    -91.6127980298  2.2686E-01  0.0000E+00  1.3834E+00  1.0000E-04   
 mr-sdci # 21  8    -91.3252201096  3.7865E-01  0.0000E+00  1.4843E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000293
time for cinew                         0.008922
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  22

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964353
   ht   2    -0.00000000    -5.68735026
   ht   3    -0.00000000     0.00000000    -5.65278534
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.25339938
   ht   5     0.00000000     0.00000000    -0.00000000     0.00000000    -4.73333664
   ht   6    -1.54065936    -0.00372874     0.61577549    -0.05819877     0.03885525    -0.48129220

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735      -2.237494E-03   1.236151E-02   2.162777E-02   5.399537E-02   0.224112    
 ref    2  -1.187064E-07   0.950391      -7.281355E-05   7.295920E-04  -8.624567E-04   1.290392E-02
 ref    3   6.925041E-03   1.468538E-03   0.945571      -3.257856E-02  -5.818713E-02   0.163610    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.903251       0.894257       1.529655E-03   6.301986E-03   7.716083E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873499    -0.00223749     0.01236151     0.02162777     0.05399537     0.22411151
 ref:   2    -0.00000012     0.95039122    -0.00007281     0.00072959    -0.00086246     0.01290392
 ref:   3     0.00692504     0.00146854     0.94557062    -0.03257856    -0.05818713     0.16361038

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 22  1    -94.6661949066  9.0887E-12  0.0000E+00  3.7398E-04  1.0000E-03   
 mr-sdci # 22  2    -94.3439147069  1.3064E-05 -7.3769E-04  6.1474E-02  1.0000E-03   
 mr-sdci # 22  3    -94.3093367226  2.4053E-09  0.0000E+00  8.0486E-04  1.0000E-03   
 mr-sdci # 22  4    -93.9099827247  3.1967E-05  0.0000E+00  3.3629E-01  1.0000E-04   
 mr-sdci # 22  5    -93.3905116210  6.2360E-04  0.0000E+00  5.9253E-01  1.0000E-04   
 mr-sdci # 22  6    -91.3345172813 -6.0077E-01  0.0000E+00  1.3778E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000541
time for cinew                         0.007360
time for eigenvalue solver             0.000142
time for vector access                 0.000001

          starting ci iteration  23

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964353
   ht   2    -0.00000000    -5.68735026
   ht   3    -0.00000000     0.00000000    -5.65278534
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.25339938
   ht   5     0.00000000     0.00000000    -0.00000000     0.00000000    -4.73333664
   ht   6    -1.54065936    -0.00372874     0.61577549    -0.05819877     0.03885525    -0.48129220
   ht   7     2.29730234     0.00374505    -0.84621477     0.09189790    -0.05748720     0.70851652    -1.04443526

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958736      -3.100673E-04  -1.235130E-02   2.260319E-02   4.991586E-02  -0.105188      -0.238708    
 ref    2  -7.012510E-07  -0.949906      -4.408759E-05  -3.191682E-04   1.318680E-03   3.483185E-02  -8.985321E-03
 ref    3   6.925560E-03   1.215563E-03  -0.945581      -3.355926E-02  -5.401902E-02   0.147927      -0.146679    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919222       0.902324       0.894276       1.637230E-03   5.411387E-03   3.416005E-02   7.857699E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873551    -0.00031007    -0.01235130     0.02260319     0.04991586    -0.10518809    -0.23870806
 ref:   2    -0.00000070    -0.94990640    -0.00004409    -0.00031917     0.00131868     0.03483185    -0.00898532
 ref:   3     0.00692556     0.00121556    -0.94558081    -0.03355926    -0.05401902     0.14792651    -0.14667897

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 23  1    -94.6661949067  4.8916E-11  0.0000E+00  3.7426E-04  1.0000E-03   
 mr-sdci # 23  2    -94.3449949387  1.0802E-03  3.4775E-04  2.6938E-02  1.0000E-03   
 mr-sdci # 23  3    -94.3093367393  1.6714E-08  0.0000E+00  7.8464E-04  1.0000E-03   
 mr-sdci # 23  4    -93.9101237129  1.4099E-04  0.0000E+00  3.3570E-01  1.0000E-04   
 mr-sdci # 23  5    -93.3922695377  1.7579E-03  0.0000E+00  5.9432E-01  1.0000E-04   
 mr-sdci # 23  6    -91.6763469264  3.4183E-01  0.0000E+00  1.4805E+00  1.0000E-04   
 mr-sdci # 23  7    -91.3293450698 -2.8345E-01  0.0000E+00  1.3636E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001317
time for cinew                         0.005828
time for eigenvalue solver             0.000174
time for vector access                 0.000000

          starting ci iteration  24

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964353
   ht   2    -0.00000000    -5.68735026
   ht   3    -0.00000000     0.00000000    -5.65278534
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.25339938
   ht   5     0.00000000     0.00000000    -0.00000000     0.00000000    -4.73333664
   ht   6    -1.54065936    -0.00372874     0.61577549    -0.05819877     0.03885525    -0.48129220
   ht   7     2.29730234     0.00374505    -0.84621477     0.09189790    -0.05748720     0.70851652    -1.04443526
   ht   8    -0.06610599     0.02250972    -0.04903588    -0.00991885     0.00007665    -0.01244196     0.01915858    -0.00206083

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958736      -2.043752E-06   1.234867E-02   2.276633E-02   4.816801E-02   1.265169E-02   0.247806       8.834258E-02
 ref    2  -8.214377E-07   0.948524      -3.054354E-05   9.415006E-04  -8.163685E-03   0.113596       1.281777E-02  -2.308386E-02
 ref    3   6.925922E-03   1.500505E-04   0.945592      -3.430084E-02  -4.794066E-02  -3.638801E-02  -3.365182E-02   0.255922    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919222       0.899697       0.894296       1.695740E-03   4.685110E-03   1.438818E-02   6.270437E-02   7.383325E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873560    -0.00000204     0.01234867     0.02276633     0.04816801     0.01265169     0.24780563     0.08834258
 ref:   2    -0.00000082     0.94852361    -0.00003054     0.00094150    -0.00816368     0.11359590     0.01281777    -0.02308386
 ref:   3     0.00692592     0.00015005     0.94559177    -0.03430084    -0.04794066    -0.03638801    -0.03365182     0.25592181

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 24  1    -94.6661949067  2.3212E-11  0.0000E+00  3.7443E-04  1.0000E-03   
 mr-sdci # 24  2    -94.3452792369  2.8430E-04  4.2413E-05  1.0120E-02  1.0000E-03   
 mr-sdci # 24  3    -94.3093367592  1.9828E-08  0.0000E+00  7.4564E-04  1.0000E-03   
 mr-sdci # 24  4    -93.9102136606  8.9948E-05  0.0000E+00  3.3519E-01  1.0000E-04   
 mr-sdci # 24  5    -93.3966236331  4.3541E-03  0.0000E+00  5.8991E-01  1.0000E-04   
 mr-sdci # 24  6    -92.6774058704  1.0011E+00  0.0000E+00  1.1957E+00  1.0000E-04   
 mr-sdci # 24  7    -91.4182867459  8.8942E-02  0.0000E+00  1.3936E+00  1.0000E-04   
 mr-sdci # 24  8    -91.1725782483 -1.5264E-01  0.0000E+00  1.4444E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.008970
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  25

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964353
   ht   2     0.00000000    -5.68872786
   ht   3     0.00000000    -0.00000000    -5.65278538
   ht   4     0.00000000    -0.00000000    -0.00000000    -5.25366228
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000    -4.74007226
   ht   6     0.00037361     0.00177819     0.00600905     0.00253935     0.00075128    -0.00008637

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958736       2.737272E-05   1.234777E-02  -2.281219E-02  -4.776883E-02   1.183280E-02
 ref    2   1.022470E-06   0.948286       1.384616E-06  -1.274251E-04   5.281103E-03  -6.894137E-02
 ref    3  -6.925611E-03  -8.849850E-05   0.945598       3.483931E-02   4.558499E-02  -6.262841E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919222       0.899246       0.894309       1.734190E-03   4.387743E-03   8.815246E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873555     0.00002737     0.01234777    -0.02281219    -0.04776883     0.01183280
 ref:   2     0.00000102     0.94828594     0.00000138    -0.00012743     0.00528110    -0.06894137
 ref:   3    -0.00692561    -0.00008850     0.94559829     0.03483931     0.04558499    -0.06262841

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 25  1    -94.6661949068  7.1948E-11  0.0000E+00  3.7356E-04  1.0000E-03   
 mr-sdci # 25  2    -94.3453102358  3.0999E-05  8.3243E-06  3.9601E-03  1.0000E-03   
 mr-sdci # 25  3    -94.3093367874  2.8284E-08  0.0000E+00  6.7992E-04  1.0000E-03   
 mr-sdci # 25  4    -93.9104035122  1.8985E-04  0.0000E+00  3.3392E-01  1.0000E-04   
 mr-sdci # 25  5    -93.3988660287  2.2424E-03  0.0000E+00  5.8613E-01  1.0000E-04   
 mr-sdci # 25  6    -91.9149072958 -7.6250E-01  0.0000E+00  1.6885E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000308
time for cinew                         0.004509
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  26

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964353
   ht   2     0.00000000    -5.68872786
   ht   3     0.00000000    -0.00000000    -5.65278538
   ht   4     0.00000000    -0.00000000    -0.00000000    -5.25366228
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000    -4.74007226
   ht   6     0.00037361     0.00177819     0.00600905     0.00253935     0.00075128    -0.00008637
   ht   7    -0.00089737     0.00154016    -0.00730891    -0.00169514    -0.00072440    -0.00000157    -0.00003678

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958736       1.587450E-05  -1.234727E-02  -2.277229E-02   4.493904E-02  -2.072430E-02  -2.326114E-02
 ref    2   8.771428E-07  -0.947985       4.478976E-06   2.133705E-04  -2.414236E-02  -0.120911       3.293270E-03
 ref    3  -6.925499E-03  -1.123727E-04  -0.945601       3.464769E-02  -3.571119E-02   6.148267E-02   0.110813    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919222       0.898676       0.894313       1.719085E-03   3.877660E-03   1.882909E-02   1.283144E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873553     0.00001587    -0.01234727    -0.02277229     0.04493904    -0.02072430    -0.02326114
 ref:   2     0.00000088    -0.94798516     0.00000448     0.00021337    -0.02414236    -0.12091101     0.00329327
 ref:   3    -0.00692550    -0.00011237    -0.94560058     0.03464769    -0.03571119     0.06148267     0.11081296

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 26  1    -94.6661949068  2.9114E-12  0.0000E+00  3.7339E-04  1.0000E-03   
 mr-sdci # 26  2    -94.3453180534  7.8176E-06  2.2194E-06  2.2254E-03  1.0000E-03   
 mr-sdci # 26  3    -94.3093367885  1.0239E-09  0.0000E+00  6.8179E-04  1.0000E-03   
 mr-sdci # 26  4    -93.9104088261  5.3139E-06  0.0000E+00  3.3386E-01  1.0000E-04   
 mr-sdci # 26  5    -93.4028279054  3.9619E-03  0.0000E+00  6.1650E-01  1.0000E-04   
 mr-sdci # 26  6    -93.2395167235  1.3246E+00  0.0000E+00  1.1716E+00  1.0000E-04   
 mr-sdci # 26  7    -91.3564350034 -6.1852E-02  0.0000E+00  1.6415E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000380
time for cinew                         0.004765
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  27

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1    -6.00964353
   ht   2     0.00000000    -5.68872786
   ht   3     0.00000000    -0.00000000    -5.65278538
   ht   4     0.00000000    -0.00000000    -0.00000000    -5.25366228
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000    -4.74007226
   ht   6     0.00037361     0.00177819     0.00600905     0.00253935     0.00075128    -0.00008637
   ht   7    -0.00089737     0.00154016    -0.00730891    -0.00169514    -0.00072440    -0.00000157    -0.00003678
   ht   8    -0.02730426     0.00140871     0.00284065    -0.00184286     0.00060995     0.00000719     0.00000304    -0.00013752

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735      -8.870795E-05  -1.235178E-02  -2.247664E-02   6.981242E-02  -1.769452E-02  -0.205829      -1.624692E-02
 ref    2   1.016469E-06  -0.947957       6.427416E-06   6.647208E-05  -0.102652      -4.931755E-02  -4.910117E-02   4.916324E-03
 ref    3  -6.925456E-03  -1.104876E-04  -0.945600       3.466515E-02   1.572564E-02   5.947004E-02   5.255713E-02   0.109107    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919220       0.898623       0.894313       1.706876E-03   1.565846E-02   6.282003E-03   4.753896E-02   1.219252E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873473    -0.00008871    -0.01235178    -0.02247664     0.06981242    -0.01769452    -0.20582949    -0.01624692
 ref:   2     0.00000102    -0.94795721     0.00000643     0.00006647    -0.10265180    -0.04931755    -0.04910117     0.00491632
 ref:   3    -0.00692546    -0.00011049    -0.94560049     0.03466515     0.01572564     0.05947004     0.05255713     0.10910721

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 27  1    -94.6661949068  4.1493E-11  0.0000E+00  3.7235E-04  1.0000E-03   
 mr-sdci # 27  2    -94.3453186707  6.1727E-07 -1.0171E-05  2.0112E-03  1.0000E-03   
 mr-sdci # 27  3    -94.3093367896  1.1244E-09  0.0000E+00  6.7250E-04  1.0000E-03   
 mr-sdci # 27  4    -93.9104125193  3.6932E-06  0.0000E+00  3.3390E-01  1.0000E-04   
 mr-sdci # 27  5    -93.4535011011  5.0673E-02  0.0000E+00  9.2166E-01  1.0000E-04   
 mr-sdci # 27  6    -93.3833477994  1.4383E-01  0.0000E+00  7.2350E-01  1.0000E-04   
 mr-sdci # 27  7    -91.4386368828  8.2202E-02  0.0000E+00  1.4260E+00  1.0000E-04   
 mr-sdci # 27  8    -91.3563394085  1.8376E-01  0.0000E+00  1.6363E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001337
time for cinew                         0.009989
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  28

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.00964353
   ht   2     0.00000000    -5.68876729
   ht   3    -0.00000000    -0.00000000    -5.65278541
   ht   4    -0.00000000     0.00000000    -0.00000000    -5.25386114
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000    -4.79694972
   ht   6     0.18076498    -0.00083687    -0.05345823     0.00821695     0.00793408    -0.00617930

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735      -6.124590E-06   1.234809E-02   2.269747E-02  -1.816040E-02  -0.260736    
 ref    2  -1.162673E-06   0.947925      -8.254827E-06   7.993977E-05   0.106055      -5.378894E-03
 ref    3   6.925242E-03   7.466511E-05   0.945599      -3.457418E-02   2.908376E-03  -9.283297E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.898561       0.894310       1.710556E-03   1.158597E-02   7.663001E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873530    -0.00000612     0.01234809     0.02269747    -0.01816040    -0.26073573
 ref:   2    -0.00000116     0.94792462    -0.00000825     0.00007994     0.10605521    -0.00537889
 ref:   3     0.00692524     0.00007467     0.94559907    -0.03457418     0.00290838    -0.09283297

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 28  1    -94.6661949068  1.9109E-11  0.0000E+00  3.7313E-04  1.0000E-03   
 mr-sdci # 28  2    -94.3453191758  5.0517E-07  1.2148E-06  1.6887E-03  1.0000E-03   
 mr-sdci # 28  3    -94.3093367904  7.6159E-10  0.0000E+00  6.7635E-04  1.0000E-03   
 mr-sdci # 28  4    -93.9104153933  2.8740E-06  0.0000E+00  3.3389E-01  1.0000E-04   
 mr-sdci # 28  5    -93.5399977718  8.6497E-02  0.0000E+00  8.0196E-01  1.0000E-04   
 mr-sdci # 28  6    -91.3689152886 -2.0144E+00  0.0000E+00  1.3515E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000415
time for cinew                         0.009611
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  29

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        2840 2x:        1009 4x:         261
All internal counts: zz :         753 yy:        5352 xx:        5876 ww:        8796
One-external counts: yz :        2237 yx:        3877 yw:        4707
Two-external counts: yy :        1416 ww:        2108 xx:        1554 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:        1034
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.00964353
   ht   2     0.00000000    -5.68876729
   ht   3    -0.00000000    -0.00000000    -5.65278541
   ht   4    -0.00000000     0.00000000    -0.00000000    -5.25386114
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000    -4.79694972
   ht   6     0.18076498    -0.00083687    -0.05345823     0.00821695     0.00793408    -0.00617930
   ht   7    -0.01974831     0.00108673     0.00205286    -0.00143714    -0.00213701     0.00064256    -0.00007192

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735       1.295945E-05   1.234876E-02   2.259949E-02   2.297308E-02   8.540722E-02  -0.246218    
 ref    2  -1.265977E-06   0.947899      -1.002387E-05   3.475593E-04  -9.539000E-02   5.418419E-02   1.079894E-02
 ref    3   6.924226E-03  -2.592250E-05   0.945596      -3.422374E-02  -4.416436E-02  -0.119501      -0.136189    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.898512       0.894303       1.682122E-03   1.157751E-02   2.451092E-02   7.928748E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 29  1    -94.6661949069  1.0783E-10  0.0000E+00  3.7264E-04  1.0000E-03   
 mr-sdci # 29  2    -94.3453200788  9.0298E-07  2.1078E-07  6.5006E-04  1.0000E-03   
 mr-sdci # 29  3    -94.3093367914  1.0534E-09  0.0000E+00  6.7046E-04  1.0000E-03   
 mr-sdci # 29  4    -93.9104220150  6.6217E-06  0.0000E+00  3.3389E-01  1.0000E-04   
 mr-sdci # 29  5    -93.6815083412  1.4151E-01  0.0000E+00  6.3315E-01  1.0000E-04   
 mr-sdci # 29  6    -91.6164258733  2.4751E-01  0.0000E+00  1.6025E+00  1.0000E-04   
 mr-sdci # 29  7    -91.3433001637 -9.5337E-02  0.0000E+00  1.3578E+00  1.0000E-04   
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000448
time for cinew                         0.008037
time for eigenvalue solver             0.000201
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after 29 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 29  1    -94.6661949069  1.0783E-10  0.0000E+00  3.7264E-04  1.0000E-03   
 mr-sdci # 29  2    -94.3453200788  9.0298E-07  2.1078E-07  6.5006E-04  1.0000E-03   
 mr-sdci # 29  3    -94.3093367914  1.0534E-09  0.0000E+00  6.7046E-04  1.0000E-03   
 mr-sdci # 29  4    -93.9104220150  6.6217E-06  0.0000E+00  3.3389E-01  1.0000E-04   
 mr-sdci # 29  5    -93.6815083412  1.4151E-01  0.0000E+00  6.3315E-01  1.0000E-04   
 mr-sdci # 29  6    -91.6164258733  2.4751E-01  0.0000E+00  1.6025E+00  1.0000E-04   
 mr-sdci # 29  7    -91.3433001637 -9.5337E-02  0.0000E+00  1.3578E+00  1.0000E-04   

####################CIUDGINFO####################

   ci vector at position   1 energy=  -94.666194906939
   ci vector at position   2 energy=  -94.345320078822
   ci vector at position   3 energy=  -94.309336791412

################END OF CIUDGINFO################

 
    3 of the   8 expansion vectors are transformed.
    3 of the   7 matrix-vector products are transformed.

    3 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =       -94.6661949069

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1 -0.906930                        +-   +-   +-   +-   +-   +-      
 z*  1  1       2 -0.297229                        +-   +-   +-   +-   +-   +     - 
 z*  1  1       3  0.092699                        +-   +-   +-   +-   +-        +- 
 z*  1  1       6  0.011528                        +-   +-   +-   +-        +-   +- 
 z   1  1      10  0.012282                        +-   +-   +-        +-   +-   +- 
 z   1  1      19 -0.013638                        +-   +    +-    -   +-   +-   +- 
 z   1  1      21  0.024954                        +-        +-   +-   +-   +-   +- 
 y   1  1      34 -0.019048              6( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      65 -0.039739             10( a  )   +-   +-   +-   +-   +-         - 
 y   1  1      73  0.016610             18( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     117  0.032545              8( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     131  0.010128             22( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     165  0.012058              2( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     276 -0.036166              5( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     278 -0.032409              7( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     282 -0.026649             11( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     286  0.022947             15( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     409  0.012181              3( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     599  0.024566              4( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     609  0.019636             14( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1    1082  0.021557              1( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1084 -0.011188              3( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1088 -0.012772              7( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1092 -0.010991             11( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1094  0.014933             13( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1379  0.015586              1( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1381 -0.013478              3( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1383 -0.014879              5( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1387  0.020435              9( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1391  0.012606             13( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1399  0.016026             21( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1424  0.012416             19( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1425 -0.017204             20( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1477  0.021897             18( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1559 -0.014838             19( a  )   +-   +     -   +-   +-   +-    - 
 y   1  1    1730 -0.014343              1( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1732 -0.015818              3( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1734 -0.012981              5( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1736 -0.015623              7( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1740 -0.014660             11( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1745 -0.014121             16( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    2108 -0.013279              1( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2110 -0.016844              3( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2114 -0.012585              7( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2118 -0.022997             11( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2131 -0.012671             24( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2213 -0.010038             25( a  )   +    +-   +-    -   +-   +-    - 
 y   1  1    2289  0.012400             20( a  )   +    +-    -   +-   +-   +-    - 
 x   1  1    2526  0.013885    4( a  )   6( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2573  0.010762    6( a  )  12( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2692 -0.010426    9( a  )  20( a  )   +-   +-   +-   +-    -    -      
 x   1  1    3586 -0.015097    6( a  )   7( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3606 -0.011037    5( a  )  10( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3666 -0.012648   10( a  )  15( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3710  0.010258    9( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    4275  0.011954    2( a  )   5( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4286  0.010173    4( a  )   7( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4293  0.013989    5( a  )   8( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4295  0.013712    7( a  )   8( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4366  0.010617    8( a  )  15( a  )   +-   +-   +-    -    -   +-      
 x   1  1    6385  0.011384    2( a  )   6( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6411  0.013286    2( a  )  10( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6413  0.011862    4( a  )  10( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6461 -0.018063   10( a  )  14( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6537 -0.012434   11( a  )  19( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6555 -0.010352   11( a  )  20( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6582  0.012559   19( a  )  21( a  )   +-   +-    -   +-   +-    -      
 x   1  1    7080 -0.010342    2( a  )   4( a  )   +-   +-    -   +-    -   +-      
 x   1  1    9198 -0.014528    2( a  )   7( a  )   +-   +-    -    -   +-   +-      
 x   1  1    9243  0.011428    7( a  )  12( a  )   +-   +-    -    -   +-   +-      
 x   1  1   11430  0.010804    7( a  )  18( a  )   +-    -   +-   +-   +-    -      
 x   1  1   18422  0.010811   10( a  )  16( a  )    -   +-   +-   +-   +-    -      
 w   1  1   27454  0.010574    6( a  )   6( a  )   +-   +-   +-   +-   +-           
 w   1  1   27484  0.010128    6( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27488  0.024439   10( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27604  0.012624   18( a  )  18( a  )   +-   +-   +-   +-   +-           
 w   1  1   27623  0.011105   19( a  )  19( a  )   +-   +-   +-   +-   +-           
 w   1  1   27830 -0.011076    4( a  )   6( a  )   +-   +-   +-   +-   +     -      
 w   1  1   28570  0.014946    2( a  )   2( a  )   +-   +-   +-   +-        +-      
 w   1  1   28573  0.016968    3( a  )   3( a  )   +-   +-   +-   +-        +-      
 w   1  1   28575 -0.015276    2( a  )   4( a  )   +-   +-   +-   +-        +-      
 w   1  1   28577  0.016398    4( a  )   4( a  )   +-   +-   +-   +-        +-      
 w   1  1   28603  0.016459    8( a  )   8( a  )   +-   +-   +-   +-        +-      
 w   1  1   28612  0.014765    9( a  )   9( a  )   +-   +-   +-   +-        +-      
 w   1  1   28635  0.017480    2( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   28637 -0.015274    4( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   28645  0.016730   12( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   28648 -0.012707    3( a  )  13( a  )   +-   +-   +-   +-        +-      
 w   1  1   28654  0.010766    9( a  )  13( a  )   +-   +-   +-   +-        +-      
 w   1  1   28662  0.010757    4( a  )  14( a  )   +-   +-   +-   +-        +-      
 w   1  1   28670 -0.010903   12( a  )  14( a  )   +-   +-   +-   +-        +-      
 w   1  1   29728 -0.010101    6( a  )   7( a  )   +-   +-   +-   +    +-    -      
 w   1  1   29751  0.010679    5( a  )  10( a  )   +-   +-   +-   +    +-    -      
 w   1  1   29766  0.010106   10( a  )  11( a  )   +-   +-   +-   +    +-    -      
 w   1  1   29816 -0.011941   10( a  )  15( a  )   +-   +-   +-   +    +-    -      
 w   1  1   30462  0.016747    2( a  )   3( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30466 -0.013552    3( a  )   4( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30501  0.013310    8( a  )   9( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30520 -0.010047    8( a  )  11( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30526  0.013442    3( a  )  12( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30547 -0.011495   12( a  )  13( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30561  0.010666   13( a  )  14( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30564  0.011713    2( a  )  15( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30574  0.012291   12( a  )  15( a  )   +-   +-   +-   +     -   +-      
 w   1  1   32726  0.010276    1( a  )   1( a  )   +-   +-   +-        +-   +-      
 w   1  1   32731  0.010874    3( a  )   3( a  )   +-   +-   +-        +-   +-      
 w   1  1   32740  0.016088    5( a  )   5( a  )   +-   +-   +-        +-   +-      
 w   1  1   32753  0.016451    7( a  )   7( a  )   +-   +-   +-        +-   +-      
 w   1  1   32835 -0.011881    5( a  )  15( a  )   +-   +-   +-        +-   +-      
 w   1  1   32845  0.013273   15( a  )  15( a  )   +-   +-   +-        +-   +-      
 w   1  1   32896  0.010486   18( a  )  18( a  )   +-   +-   +-        +-   +-      
 w   1  1   35040 -0.012688    2( a  )  10( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35042 -0.011914    4( a  )  10( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35094 -0.018419   10( a  )  14( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35757  0.016028    2( a  )   4( a  )   +-   +-   +    +-    -   +-      
 w   1  1   35842  0.016294    2( a  )  14( a  )   +-   +-   +    +-    -   +-      
 w   1  1   38019  0.013271    1( a  )   2( a  )   +-   +-   +     -   +-   +-      
 w   1  1   38113 -0.010189    5( a  )  14( a  )   +-   +-   +     -   +-   +-      
 w   1  1   38136  0.014808   14( a  )  15( a  )   +-   +-   +     -   +-   +-      
 w   1  1   40286  0.012406    1( a  )   1( a  )   +-   +-        +-   +-   +-      
 w   1  1   40288  0.020010    2( a  )   2( a  )   +-   +-        +-   +-   +-      
 w   1  1   40293  0.011795    2( a  )   4( a  )   +-   +-        +-   +-   +-      
 w   1  1   40351  0.012209   11( a  )  11( a  )   +-   +-        +-   +-   +-      
 w   1  1   40353  0.011245    2( a  )  12( a  )   +-   +-        +-   +-   +-      
 w   1  1   40378  0.021548    2( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1   40380  0.014460    4( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1   40390  0.021031   14( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1   42571 -0.012034    3( a  )   6( a  )   +-   +    +-   +-   +-    -      
 w   1  1   42599 -0.011025    1( a  )  10( a  )   +-   +    +-   +-   +-    -      
 w   1  1   43314 -0.010031    2( a  )   3( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43318  0.017129    3( a  )   4( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43353 -0.015033    8( a  )   9( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43378 -0.014323    3( a  )  12( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43389  0.012463    2( a  )  13( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43399  0.015794   12( a  )  13( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43403  0.011239    3( a  )  14( a  )   +-   +    +-   +-    -   +-      
 w   1  1   45581  0.013836    1( a  )   3( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45583 -0.010701    3( a  )   3( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45588 -0.012281    1( a  )   5( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45601  0.010822    3( a  )   7( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45658  0.015346    3( a  )  13( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45683  0.016330    1( a  )  15( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45695  0.019049   13( a  )  15( a  )   +-   +    +-    -   +-   +-      
 w   1  1   47847  0.012252    1( a  )   2( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47850  0.011473    2( a  )   3( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47937  0.013367    1( a  )  14( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47949  0.010637   13( a  )  14( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47964  0.012347   14( a  )  15( a  )   +-   +     -   +-   +-   +-      
 w   1  1   50119  0.010620    3( a  )   3( a  )   +-        +-   +-   +-   +-      
 w   1  1   50204  0.010100   13( a  )  13( a  )   +-        +-   +-   +-   +-      
 w   1  1   50233  0.010284   15( a  )  15( a  )   +-        +-   +-   +-   +-      
 w   1  1   52427  0.013947    1( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52429  0.010500    3( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52431  0.010371    5( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52446  0.011368   10( a  )  11( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52511  0.014716   10( a  )  16( a  )   +    +-   +-   +-   +-    -      
 w   1  1   57765 -0.014438    1( a  )  14( a  )   +    +-    -   +-   +-   +-      
 w   1  1   57775 -0.011989   11( a  )  14( a  )   +    +-    -   +-   +-   +-      
 w   1  1   57807 -0.014396   14( a  )  16( a  )   +    +-    -   +-   +-   +-      

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             157
     0.01> rq > 0.001           3808
    0.001> rq > 0.0001          6448
   0.0001> rq > 0.00001         5167
  0.00001> rq > 0.000001        2554
 0.000001> rq                  46341
           all                 64477
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         277 2x:           0 4x:           0
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:         201
  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.906929890809     85.612507469055
     2     2     -0.297228592220     28.052878805194
     3     3      0.092699428979     -8.766042535747
     4     4     -0.000003263084      0.000308130388
     5     5      0.000004009641     -0.000378851159
     6     6      0.011528160882     -1.092951158148

 number of reference csfs (nref) is     6.  root number (iroot) is  1.
 c0**2 =   0.91959275  c**2 (all zwalks) =   0.92070333

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.398163757370   "relaxed" cnot**2         =   0.919592745530
 eci       =    -94.666194906939   deltae = eci - eref       =  -0.268031149569
 eci+dv1   =    -94.687746555789   dv1 = (1-cnot**2)*deltae  =  -0.021551648849
 eci+dv2   =    -94.689630986611   dv2 = dv1 / cnot**2       =  -0.023436079671
 eci+dv3   =    -94.691876534040   dv3 = dv1 / (2*cnot**2-1) =  -0.025681627100
 eci+pople =    -94.686909115325   ( 12e- scaled deltae )    =  -0.288745357955


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 2) =       -94.3453200788

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       4 -0.917471                        +-   +-   +-   +-   +    +-    - 
 z*  1  1       5 -0.238293                        +-   +-   +-   +-   +     -   +- 
 z   1  1      11 -0.081833                        +-   +-   +    +-   +-   +-    - 
 z   1  1      12  0.066205                        +-   +-   +    +-   +-    -   +- 
 y   1  1      30  0.020989              2( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      88 -0.017984              6( a  )   +-   +-   +-   +-    -   +-      
 y   1  1      92  0.018006             10( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     115  0.027611              6( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     119  0.011454             10( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     127 -0.026436             18( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     146  0.028231             10( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     154 -0.013968             18( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     173  0.014122             10( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     192  0.018228              2( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     194 -0.021750              4( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     198 -0.021790              8( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     202  0.019817             12( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     204 -0.017314             14( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     207 -0.018186             17( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     225  0.012451              8( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     330  0.021339              5( a  )   +-   +-   +-    -   +    +-    - 
 y   1  1     336  0.010664             11( a  )   +-   +-   +-    -   +    +-    - 
 y   1  1     357 -0.023211              5( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     359 -0.024214              7( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     363 -0.017871             11( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     367  0.017269             15( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     434 -0.020459              1( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     436  0.031010              3( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     442 -0.041907              9( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     446 -0.026495             13( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     448  0.016447             15( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     495 -0.010736              8( a  )   +-   +-   +-        +-   +-    - 
 y   1  1     586 -0.012698             18( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     653 -0.013073              4( a  )   +-   +-    -   +-   +    +-    - 
 y   1  1     680  0.018022              4( a  )   +-   +-    -   +-   +     -   +- 
 y   1  1     690  0.014346             14( a  )   +-   +-    -   +-   +     -   +- 
 y   1  1     735  0.014902              5( a  )   +-   +-    -   +    +-   +-    - 
 y   1  1     821 -0.013329             10( a  )   +-   +-   +    +-   +-    -    - 
 y   1  1     840 -0.010720              2( a  )   +-   +-   +    +-    -   +-    - 
 y   1  1     846  0.025387              8( a  )   +-   +-   +    +-    -   +-    - 
 y   1  1     897 -0.014330              5( a  )   +-   +-   +     -   +-   +-    - 
 y   1  1     899 -0.017272              7( a  )   +-   +-   +     -   +-   +-    - 
 y   1  1     903 -0.018620             11( a  )   +-   +-   +     -   +-   +-    - 
 y   1  1     977  0.012864              4( a  )   +-   +-        +-   +-   +-    - 
 y   1  1     981  0.010096              8( a  )   +-   +-        +-   +-   +-    - 
 y   1  1    1136 -0.010471              1( a  )   +-    -   +-   +-   +    +-    - 
 y   1  1    1163  0.014622              1( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1165 -0.010263              3( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1169 -0.010021              7( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1175  0.012613             13( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1298  0.013625              1( a  )   +-    -   +    +-   +-   +-    - 
 y   1  1    1304 -0.010349              7( a  )   +-    -   +    +-   +-   +-    - 
 y   1  1    1397 -0.019924             19( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1408 -0.038918              3( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1410  0.011802              5( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1414  0.041668              9( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1416 -0.022308             11( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1418  0.017835             13( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1420 -0.027439             15( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1426  0.011827             21( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1429 -0.013170             24( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1437  0.012055              5( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    1461 -0.017793              2( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1467  0.022409              8( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1471 -0.017301             12( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1531  0.013935             18( a  )   +-   +    +-    -    -   +-   +- 
 y   1  1    1551 -0.013867             11( a  )   +-   +     -   +-   +-   +-    - 
 y   1  1    1613 -0.011159             19( a  )   +-   +     -   +-    -   +-   +- 
 y   1  1    1784  0.010542              1( a  )    -   +-   +-   +-   +    +-    - 
 y   1  1    1813 -0.011100              3( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1815 -0.010033              5( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1817 -0.012123              7( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1821 -0.010058             11( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1826 -0.010452             16( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    2127  0.010318             20( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2135  0.011903              1( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2143  0.016187              9( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2147  0.011331             13( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2155  0.010429             21( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2162  0.011407              1( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    2172  0.014956             11( a  )   +    +-   +-   +-    -    -   +- 
 x   1  1    4288 -0.013925    6( a  )   7( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4408 -0.013087    5( a  )  18( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4414 -0.010901   11( a  )  18( a  )   +-   +-   +-    -    -   +-      
 x   1  1    5341 -0.012742    6( a  )   7( a  )   +-   +-   +-    -   +     -    - 
 x   1  1    5421 -0.010741   10( a  )  15( a  )   +-   +-   +-    -   +     -    - 
 x   1  1    7514  0.011946   10( a  )  14( a  )   +-   +-    -   +-    -   +     - 
 x   1  1    8166  0.011010    2( a  )  10( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8168  0.010098    4( a  )  10( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8216 -0.014931   10( a  )  14( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8292 -0.010396   11( a  )  19( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8337  0.010460   19( a  )  21( a  )   +-   +-    -   +-   +     -    - 
 x   1  1   10251 -0.014961    2( a  )   7( a  )   +-   +-    -    -   +    +-    - 
 x   1  1   10296  0.011835    7( a  )  12( a  )   +-   +-    -    -   +    +-    - 
 w   1  1   27832  0.023900    6( a  )   6( a  )   +-   +-   +-   +-   +     -      
 w   1  1   27856  0.010757    9( a  )   9( a  )   +-   +-   +-   +-   +     -      
 w   1  1   27862 -0.012060    6( a  )  10( a  )   +-   +-   +-   +-   +     -      
 w   1  1   27866 -0.013647   10( a  )  10( a  )   +-   +-   +-   +-   +     -      
 w   1  1   28020 -0.012098   19( a  )  20( a  )   +-   +-   +-   +-   +     -      
 w   1  1   28210 -0.012564    6( a  )   6( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28240 -0.011521    6( a  )  10( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28244 -0.028671   10( a  )  10( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28255 -0.011356   11( a  )  11( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28360 -0.014009   18( a  )  18( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28379 -0.011975   19( a  )  19( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28584  0.014969    2( a  )   6( a  )   +-   +-   +-   +-        +-      
 w   1  1   28586 -0.017274    4( a  )   6( a  )   +-   +-   +-   +-        +-      
 w   1  1   28601  0.010442    6( a  )   8( a  )   +-   +-   +-   +-        +-      
 w   1  1   28639  0.017216    6( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   30484 -0.012921    6( a  )   7( a  )   +-   +-   +-   +     -   +-      
 w   1  1   31974  0.010802    2( a  )   3( a  )   +-   +-   +-   +         +-    - 
 w   1  1   33860  0.010859    1( a  )   1( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33865  0.012153    3( a  )   3( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33874  0.016378    5( a  )   5( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33887  0.016603    7( a  )   7( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33938  0.010325    1( a  )  13( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33940 -0.010078    3( a  )  13( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33950  0.010361   13( a  )  13( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33967  0.010353    3( a  )  15( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33969 -0.011744    5( a  )  15( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33979  0.013282   15( a  )  15( a  )   +-   +-   +-        +    +-    - 
 w   1  1   34030  0.010313   18( a  )  18( a  )   +-   +-   +-        +    +-    - 
 w   1  1   36228  0.012017   10( a  )  14( a  )   +-   +-   +    +-    -   +     - 
 w   1  1   36930 -0.010256    2( a  )  10( a  )   +-   +-   +    +-   +     -    - 
 w   1  1   36984 -0.014439   10( a  )  14( a  )   +-   +-   +    +-   +     -    - 
 w   1  1   37269  0.012650    2( a  )   4( a  )   +-   +-   +    +-        +-    - 
 w   1  1   37354  0.011503    2( a  )  14( a  )   +-   +-   +    +-        +-    - 
 w   1  1   39153  0.013656    1( a  )   2( a  )   +-   +-   +     -   +    +-    - 
 w   1  1   39231  0.010062    2( a  )  13( a  )   +-   +-   +     -   +    +-    - 
 w   1  1   39270  0.014057   14( a  )  15( a  )   +-   +-   +     -   +    +-    - 
 w   1  1   41420  0.012540    1( a  )   1( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41422  0.020664    2( a  )   2( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41427  0.010356    2( a  )   4( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41455  0.011323    8( a  )   8( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41485  0.012068   11( a  )  11( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41487  0.011817    2( a  )  12( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41512  0.019954    2( a  )  14( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41514  0.013799    4( a  )  14( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41524  0.019692   14( a  )  14( a  )   +-   +-        +-   +    +-    - 
 w   1  1   43327 -0.014344    3( a  )   6( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43393  0.012063    6( a  )  13( a  )   +-   +    +-   +-    -   +-      
 w   1  1   44461 -0.011041    3( a  )   6( a  )   +-   +    +-   +-   +     -    - 
 w   1  1   44830  0.011361    3( a  )   4( a  )   +-   +    +-   +-        +-    - 
 w   1  1   44865 -0.010905    8( a  )   9( a  )   +-   +    +-   +-        +-    - 
 w   1  1   44911  0.010534   12( a  )  13( a  )   +-   +    +-   +-        +-    - 
 w   1  1   46715  0.014925    1( a  )   3( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46717 -0.012316    3( a  )   3( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46722 -0.012706    1( a  )   5( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46735  0.010885    3( a  )   7( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46792  0.016638    3( a  )  13( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46817  0.016877    1( a  )  15( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46829  0.019718   13( a  )  15( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   48981  0.012235    1( a  )   2( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   48984  0.012806    2( a  )   3( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49042 -0.011165    8( a  )  11( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49071  0.012876    1( a  )  14( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49083  0.010421   13( a  )  14( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49098  0.011689   14( a  )  15( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   51248  0.010041    1( a  )   1( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51253  0.012387    3( a  )   3( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51326  0.010186    1( a  )  13( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51338  0.010574   13( a  )  13( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51355  0.010343    3( a  )  15( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51367  0.010923   15( a  )  15( a  )   +-        +-   +-   +    +-    - 
 w   1  1   54317  0.011905    1( a  )  10( a  )   +    +-   +-   +-   +     -    - 
 w   1  1   54401  0.012985   10( a  )  16( a  )   +    +-   +-   +-   +     -    - 
 w   1  1   58899 -0.014341    1( a  )  14( a  )   +    +-    -   +-   +    +-    - 
 w   1  1   58909 -0.011874   11( a  )  14( a  )   +    +-    -   +-   +    +-    - 
 w   1  1   58941 -0.014100   14( a  )  16( a  )   +    +-    -   +-   +    +-    - 
 w   1  1   61181 -0.010172    1( a  )  15( a  )   +     -   +-   +-   +    +-    - 

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             168
     0.01> rq > 0.001           4018
    0.001> rq > 0.0001          6461
   0.0001> rq > 0.00001         3802
  0.00001> rq > 0.000001        3987
 0.000001> rq                  46039
           all                 64477
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         277 2x:           0 4x:           0
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:         201
  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.000013529878      0.001270500751
     2     2      0.000123941690     -0.011649135293
     3     3      0.000085797474     -0.008023078793
     4     4     -0.917471158522     86.291071501211
     5     5     -0.238292856545     22.413636966618
     6     6     -0.000052615016      0.004927996806

 number of reference csfs (nref) is     6.  root number (iroot) is  2.
 c0**2 =   0.89853684  c**2 (all zwalks) =   0.90961660

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.053552136190   "relaxed" cnot**2         =   0.898536837875
 eci       =    -94.345320078822   deltae = eci - eref       =  -0.291767942632
 eci+dv1   =    -94.374923776888   dv1 = (1-cnot**2)*deltae  =  -0.029603698066
 eci+dv2   =    -94.378266639072   dv2 = dv1 / cnot**2       =  -0.032946560250
 eci+dv3   =    -94.382460557760   dv3 = dv1 / (2*cnot**2-1) =  -0.037140478938
 eci+pople =    -94.374949423630   ( 12e- scaled deltae )    =  -0.321397287440


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 3) =       -94.3093367914

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1  0.304317                        +-   +-   +-   +-   +-   +-      
 z*  1  1       2 -0.895546                        +-   +-   +-   +-   +-   +     - 
 z*  1  1       3  0.063689                        +-   +-   +-   +-   +-        +- 
 z*  1  1       6 -0.018465                        +-   +-   +-   +-        +-   +- 
 z   1  1      13  0.012545                        +-   +-   +    +-    -   +-   +- 
 z   1  1      19  0.014784                        +-   +    +-    -   +-   +-   +- 
 z   1  1      21 -0.022364                        +-        +-   +-   +-   +-   +- 
 z   1  1      27 -0.037082                        +     -   +-   +-   +-   +-   +- 
 y   1  1      38 -0.018712             10( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      65 -0.026075             10( a  )   +-   +-   +-   +-   +-         - 
 y   1  1      73 -0.011778             18( a  )   +-   +-   +-   +-   +-         - 
 y   1  1      86 -0.010438              4( a  )   +-   +-   +-   +-    -   +-      
 y   1  1      90  0.046390              8( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     104  0.011782             22( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     111 -0.012563              2( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     117  0.027006              8( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     121 -0.012009             12( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     123  0.010587             14( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     144  0.018729              8( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     165  0.010050              2( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     245 -0.016845              1( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     249 -0.061420              5( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     251 -0.051880              7( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     255 -0.035204             11( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     259  0.032324             15( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     272 -0.010254              1( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     276 -0.032472              5( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     278 -0.012574              7( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     280  0.011145              9( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     282 -0.016811             11( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     286  0.013575             15( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     303 -0.020580              5( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     305 -0.020784              7( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     309 -0.016476             11( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     313  0.015643             15( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     407 -0.010705              1( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     409  0.016098              3( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     415 -0.012589              9( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     572  0.044916              4( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     576  0.013192              8( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     580 -0.011165             12( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     582  0.021031             14( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     597  0.017329              2( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     599  0.013400              4( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     609  0.016538             14( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     626  0.016242              4( a  )   +-   +-    -   +-   +-        +- 
 y   1  1     636  0.013169             14( a  )   +-   +-    -   +-   +-        +- 
 y   1  1    1055  0.038883              1( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1057 -0.020886              3( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1059 -0.021206              5( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1061 -0.016458              7( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1065 -0.013945             11( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1067  0.021538             13( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1069  0.010161             15( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1082  0.019436              1( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1109  0.014731              1( a  )   +-    -   +-   +-   +-        +- 
 y   1  1    1121  0.011295             13( a  )   +-    -   +-   +-   +-        +- 
 y   1  1    1379  0.015131              1( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1381 -0.018807              3( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1387  0.014296              9( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1391  0.010868             13( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1424 -0.013851             19( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1451  0.010643             19( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    1452 -0.012069             20( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    1477 -0.010185             18( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1504  0.014428             18( a  )   +-   +    +-    -   +-    -   +- 
 y   1  1    1703 -0.022921              1( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1705 -0.019459              3( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1707 -0.014504              5( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1709 -0.018772              7( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1713 -0.013498             11( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1718 -0.017650             16( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1757 -0.010381              1( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    1759 -0.010693              3( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    1763 -0.011024              7( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    1767 -0.010275             11( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    2110 -0.015523              3( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2114 -0.017262              7( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2118 -0.023240             11( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2131 -0.011360             24( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2154  0.012385             20( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2206 -0.021111             18( a  )   +    +-   +-    -   +-   +-    - 
 y   1  1    2288  0.018686             19( a  )   +    +-    -   +-   +-   +-    - 
 y   1  1    2368 -0.011178             18( a  )   +     -   +-   +-   +-   +-    - 
 x   1  1    2526  0.010974    4( a  )   6( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2552 -0.010361    4( a  )  10( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2656  0.010286    8( a  )  18( a  )   +-   +-   +-   +-    -    -      
 x   1  1    3706 -0.012440    5( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3708 -0.013601    7( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3712 -0.014153   11( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    4059  0.010084    7( a  )  18( a  )   +-   +-   +-    -   +-         - 
 x   1  1    4626  0.011817    2( a  )   5( a  )   +-   +-   +-    -    -   +     - 
 x   1  1    4644  0.011097    5( a  )   8( a  )   +-   +-   +-    -    -   +     - 
 x   1  1    4646  0.010958    7( a  )   8( a  )   +-   +-   +-    -    -   +     - 
 x   1  1    6411 -0.014323    2( a  )  10( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6461  0.013398   10( a  )  14( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6537  0.010624   11( a  )  19( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6762  0.011066    2( a  )  10( a  )   +-   +-    -   +-   +-         - 
 x   1  1    6812 -0.014779   10( a  )  14( a  )   +-   +-    -   +-   +-         - 
 x   1  1    6888 -0.010049   11( a  )  19( a  )   +-   +-    -   +-   +-         - 
 x   1  1    9549 -0.014259    2( a  )   7( a  )   +-   +-    -    -   +-   +     - 
 x   1  1    9594  0.010833    7( a  )  12( a  )   +-   +-    -    -   +-   +     - 
 w   1  1   27434 -0.010513    1( a  )   1( a  )   +-   +-   +-   +-   +-           
 w   1  1   27454  0.011928    6( a  )   6( a  )   +-   +-   +-   +-   +-           
 w   1  1   27484 -0.016140    6( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27488 -0.034359   10( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27489 -0.012148    1( a  )  11( a  )   +-   +-   +-   +-   +-           
 w   1  1   27499 -0.011735   11( a  )  11( a  )   +-   +-   +-   +-   +-           
 w   1  1   27642 -0.015839   19( a  )  20( a  )   +-   +-   +-   +-   +-           
 w   1  1   27654 -0.012609   11( a  )  21( a  )   +-   +-   +-   +-   +-           
 w   1  1   27751 -0.011567   18( a  )  25( a  )   +-   +-   +-   +-   +-           
 w   1  1   28948  0.015142    2( a  )   2( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28951  0.016494    3( a  )   3( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28953 -0.014645    2( a  )   4( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28955  0.015987    4( a  )   4( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28981  0.016006    8( a  )   8( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28990  0.014270    9( a  )   9( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29013  0.017194    2( a  )  12( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29015 -0.014486    4( a  )  12( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29023  0.016258   12( a  )  12( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29026 -0.011674    3( a  )  13( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29032  0.010137    9( a  )  13( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29040  0.010606    4( a  )  14( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29048 -0.010288   12( a  )  14( a  )   +-   +-   +-   +-        +     - 
 w   1  1   30840  0.016873    2( a  )   3( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30844 -0.013187    3( a  )   4( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30879  0.013524    8( a  )   9( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30904  0.013133    3( a  )  12( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30925 -0.010922   12( a  )  13( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30939  0.010445   13( a  )  14( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30942  0.011254    2( a  )  15( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30952  0.011558   12( a  )  15( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   33104  0.010324    1( a  )   1( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33109  0.010701    3( a  )   3( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33118  0.013488    5( a  )   5( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33131  0.014015    7( a  )   7( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33213 -0.010340    5( a  )  15( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33223  0.012695   15( a  )  15( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33274  0.010083   18( a  )  18( a  )   +-   +-   +-        +-   +     - 
 w   1  1   35040  0.012707    2( a  )  10( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35094  0.013477   10( a  )  14( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35418 -0.011076    2( a  )  10( a  )   +-   +-   +    +-   +-         - 
 w   1  1   35472 -0.014848   10( a  )  14( a  )   +-   +-   +    +-   +-         - 
 w   1  1   36135  0.015510    2( a  )   4( a  )   +-   +-   +    +-    -   +     - 
 w   1  1   36220  0.016221    2( a  )  14( a  )   +-   +-   +    +-    -   +     - 
 w   1  1   38397  0.012607    1( a  )   2( a  )   +-   +-   +     -   +-   +     - 
 w   1  1   38514  0.014788   14( a  )  15( a  )   +-   +-   +     -   +-   +     - 
 w   1  1   40664  0.011138    1( a  )   1( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40666  0.018543    2( a  )   2( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40671  0.010439    2( a  )   4( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40729  0.011412   11( a  )  11( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40731  0.011069    2( a  )  12( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40756  0.021065    2( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40758  0.013533    4( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40768  0.021176   14( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1   42571 -0.011012    3( a  )   6( a  )   +-   +    +-   +-   +-    -      
 w   1  1   42707 -0.010499    1( a  )  18( a  )   +-   +    +-   +-   +-    -      
 w   1  1   43692 -0.010313    2( a  )   3( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43696  0.017120    3( a  )   4( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43731 -0.015315    8( a  )   9( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43756 -0.013982    3( a  )  12( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43767  0.012383    2( a  )  13( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43777  0.015350   12( a  )  13( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43781  0.011209    3( a  )  14( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   45959  0.013522    1( a  )   3( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   45961 -0.010841    3( a  )   3( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   45966 -0.011259    1( a  )   5( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   46036  0.014853    3( a  )  13( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   46061  0.015483    1( a  )  15( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   46073  0.018067   13( a  )  15( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   48225  0.010728    1( a  )   2( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48228  0.011027    2( a  )   3( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48315  0.012290    1( a  )  14( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48327  0.010254   13( a  )  14( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48342  0.012138   14( a  )  15( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   50497  0.010720    3( a  )   3( a  )   +-        +-   +-   +-   +     - 
 w   1  1   52427 -0.011994    1( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52511 -0.011394   10( a  )  16( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52805  0.010795    1( a  )  10( a  )   +    +-   +-   +-   +-         - 
 w   1  1   52889  0.011071   10( a  )  16( a  )   +    +-   +-   +-   +-         - 
 w   1  1   58143 -0.013447    1( a  )  14( a  )   +    +-    -   +-   +-   +     - 
 w   1  1   58153 -0.011437   11( a  )  14( a  )   +    +-    -   +-   +-   +     - 
 w   1  1   58185 -0.014226   14( a  )  16( a  )   +    +-    -   +-   +-   +     - 

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             181
     0.01> rq > 0.001           4078
    0.001> rq > 0.0001          6688
   0.0001> rq > 0.00001         5233
  0.00001> rq > 0.000001        3633
 0.000001> rq                  44662
           all                 64477
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         277 2x:           0 4x:           0
All internal counts: zz :         753 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        1945    task #     6:        3465    task #     7:        4159    task #     8:         726
task #     9:        5724    task #    10:        1011    task #    11:        1011    task #    12:        1341
task #    13:        1341    task #    14:         670    task #    15:         672    task #    16:         672
task #    17:         251    task #    18:        1296    task #    19:         773    task #    20:         201
  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.304317415355    -28.616578141806
     2     2     -0.895545556514     84.186729242487
     3     3      0.063688862169     -6.019082384734
     4     4     -0.000064239554      0.006040513655
     5     5     -0.000038064523      0.003570792887
     6     6     -0.018465001458      1.740653224030

 number of reference csfs (nref) is     6.  root number (iroot) is  3.
 c0**2 =   0.89900817  c**2 (all zwalks) =   0.90132415

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.011453700561   "relaxed" cnot**2         =   0.899008166099
 eci       =    -94.309336791412   deltae = eci - eref       =  -0.297883090852
 eci+dv1   =    -94.339420551046   dv1 = (1-cnot**2)*deltae  =  -0.030083759633
 eci+dv2   =    -94.342800068797   dv2 = dv1 / cnot**2       =  -0.033463277385
 eci+dv3   =    -94.347034966775   dv3 = dv1 / (2*cnot**2-1) =  -0.037698175362
 eci+pople =    -94.339418863805   ( 12e- scaled deltae )    =  -0.327965163245
maximum overlap with reference    1(overlap= 0.95874)
weight of reference states=  0.9192

 information on vector: 1 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    2(overlap= 0.94790)
weight of reference states=  0.8985

 information on vector: 2 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    3(overlap= 0.94560)
weight of reference states=  0.8943

 information on vector: 3 from unit 11 written to unit 48 filename civout              
 passed aftci ... 
 readint2: molcas,dalton2=                     0                     0
 files%faoints=aoints              
lodens (list->root)=  1  2  3
sifcfg setup: record length 4096 DP
# d1 elements per record  3272
# d2 elements per record  2730
  The MR-CISD density will be calculated.
 item #                     1 suffix=:.drt1.state1:
 read_civout: repnuc=  -88.6565513780480     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 96% root-following 0
 MR-CISD energy:   -94.66619491    -6.00964353
 residuum:     0.00037264
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821    -0.01624692
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894     0.00491632
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907     0.10910721

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821     0.00000000
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894     0.00000000
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 computing MRCISD density
 densi: densityinfo%a4den=   1.00000000000000     
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    1
--------------------------------------------------------------------------------
================================================================================
   DYZ=     215  DYX=     335  DYW=     401
   D0Z=     147  D0Y=     776  D0X=     489  D0W=     695
  DDZI=     109 DDYI=     378 DDXI=     272 DDWI=     359
  DDZE=       0 DDYE=      92 DDXE=      71 DDWE=      98
================================================================================
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

          modens reordered block   1

               a     1        a     2        a     3        a     4        a     5        a     6        a     7        a     8
  a     1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     3    0.00000        0.00000        1.98299       2.819781E-04   3.348756E-08   3.126733E-03   2.372699E-08   3.127524E-07
  a     4    0.00000        0.00000       2.819781E-04    1.97601      -1.389738E-08  -2.245894E-03   7.125979E-09   5.346487E-07
  a     5    0.00000        0.00000       3.348756E-08  -1.389738E-08    1.97353       2.609326E-09  -1.871877E-03   3.387117E-06
  a     6    0.00000        0.00000       3.126733E-03  -2.245894E-03   2.609326E-09    1.96964       8.148115E-09  -1.241083E-06
  a     7    0.00000        0.00000       2.372699E-08   7.125979E-09  -1.871877E-03   8.148115E-09    1.97270       1.768303E-06
  a     8    0.00000        0.00000       3.127524E-07   5.346487E-07   3.387117E-06  -1.241083E-06   1.768303E-06    1.85280    
  a     9    0.00000        0.00000      -1.344135E-06  -5.317879E-07  -1.089726E-05   3.430693E-06  -2.838191E-06   0.360287    
  a    10    0.00000        0.00000      -4.541249E-04   3.685993E-03  -2.835951E-07  -8.577966E-04   9.578521E-08   5.529282E-08
  a    11    0.00000        0.00000      -7.926906E-08   9.594700E-08  -5.781095E-04   3.693187E-07   6.614530E-03   3.141970E-06
  a    12    0.00000        0.00000       3.687425E-04  -2.745416E-03  -3.184425E-08  -3.776705E-03  -3.809909E-08   3.680584E-07
  a    13    0.00000        0.00000       3.467944E-08   6.548618E-08  -2.234830E-03   1.133784E-07   4.083639E-03   9.284485E-07
  a    14    0.00000        0.00000      -3.652258E-03  -1.882249E-03  -3.707260E-08  -4.420496E-03   3.859846E-07  -5.527527E-07
  a    15    0.00000        0.00000       1.012366E-06  -2.233063E-06   8.040973E-07  -3.003158E-06   1.157227E-06  -2.231005E-02
  a    16    0.00000        0.00000       3.611646E-03  -6.976284E-03  -1.719013E-06  -8.980055E-03  -2.706498E-06   7.109019E-06
  a    17    0.00000        0.00000      -7.254955E-07   1.489058E-06  -7.723118E-03   2.134945E-06  -1.205870E-02  -6.831623E-07
  a    18    0.00000        0.00000       3.693223E-03  -3.527754E-03   1.708860E-07   2.410525E-03   4.436648E-07  -1.616677E-07
  a    19    0.00000        0.00000      -1.117464E-06   2.239116E-06  -4.879100E-07   5.536188E-06  -1.778972E-06  -6.849752E-03
  a    20    0.00000        0.00000       1.164423E-03  -2.302006E-03  -1.206437E-07  -5.953564E-03   1.122327E-07  -6.805583E-06
  a    21    0.00000        0.00000      -2.128754E-09  -1.610617E-09  -5.957158E-03   3.311950E-07   1.073101E-02   9.724334E-07
  a    22    0.00000        0.00000       8.559871E-04  -2.764242E-04  -1.707495E-07   1.077472E-02  -2.429336E-07  -3.314800E-07
  a    23    0.00000        0.00000      -2.784919E-09  -9.244833E-08  -9.874412E-03  -2.062642E-07  -7.628771E-03   4.175091E-07
  a    24    0.00000        0.00000       1.005996E-03   6.477784E-03  -2.022476E-07   2.504141E-03  -3.120514E-07   1.482172E-07
  a    25    0.00000        0.00000       3.406720E-03   6.701952E-03   3.364064E-07  -5.780643E-05   7.438911E-08   2.233545E-07
  a    26    0.00000        0.00000      -1.283495E-07  -2.859841E-07   7.288095E-03  -1.755206E-07   1.706030E-03  -1.689052E-07
  a    27    0.00000        0.00000       2.523552E-07  -5.698950E-07   8.092757E-07  -3.299425E-07   3.754759E-07   3.478473E-03
  a    28    0.00000        0.00000       1.244972E-08   1.928587E-07  -3.904148E-06   7.444577E-08  -1.347252E-06   9.396832E-08
  a    29    0.00000        0.00000      -6.612335E-07   2.399148E-07  -1.366678E-06   1.019236E-07  -4.326835E-07  -2.466417E-08
  a    30    0.00000        0.00000       6.141487E-03  -2.312562E-03  -1.026256E-07  -3.503705E-04  -1.046908E-08   1.739699E-07
  a    31    0.00000        0.00000      -1.933677E-07   2.296748E-07  -6.903485E-04   4.335409E-07  -2.536684E-04   2.124918E-08
  a    32    0.00000        0.00000      -5.004970E-03   4.231009E-03   3.541829E-08   9.857056E-03   5.784195E-08  -2.533116E-07
  a    33    0.00000        0.00000       1.821964E-03   2.890489E-03  -1.269536E-07  -1.366678E-03  -9.589028E-08   7.590741E-07
  a    34    0.00000        0.00000      -3.355745E-07  -8.413513E-07   6.641386E-08  -2.902631E-07  -1.895181E-07   5.552790E-03
  a    35    0.00000        0.00000       2.360776E-04   8.742271E-04   1.176131E-07   7.125759E-04   1.058760E-07   1.739631E-07
  a    36    0.00000        0.00000      -1.599007E-08  -5.335765E-08   1.124177E-03  -8.230607E-08   1.100026E-03   1.551432E-07

               a     9        a    10        a    11        a    12        a    13        a    14        a    15        a    16
  a     3  -1.344135E-06  -4.541249E-04  -7.926906E-08   3.687425E-04   3.467944E-08  -3.652258E-03   1.012366E-06   3.611646E-03
  a     4  -5.317879E-07   3.685993E-03   9.594700E-08  -2.745416E-03   6.548618E-08  -1.882249E-03  -2.233063E-06  -6.976284E-03
  a     5  -1.089726E-05  -2.835951E-07  -5.781095E-04  -3.184425E-08  -2.234830E-03  -3.707260E-08   8.040973E-07  -1.719013E-06
  a     6   3.430693E-06  -8.577966E-04   3.693187E-07  -3.776705E-03   1.133784E-07  -4.420496E-03  -3.003158E-06  -8.980055E-03
  a     7  -2.838191E-06   9.578521E-08   6.614530E-03  -3.809909E-08   4.083639E-03   3.859846E-07   1.157227E-06  -2.706498E-06
  a     8   0.360287       5.529282E-08   3.141970E-06   3.680584E-07   9.284485E-07  -5.527527E-07  -2.231005E-02   7.109019E-06
  a     9   0.132126       1.826004E-08   3.276202E-07   4.582312E-08   1.321617E-07  -3.099888E-07  -7.887187E-03   2.615266E-06
  a    10   1.826004E-08   8.669736E-03   5.239415E-08   1.182337E-04  -2.480224E-09  -2.141927E-04   7.572797E-07   2.352432E-03
  a    11   3.276202E-07   5.239415E-08   1.018259E-02   2.948157E-08  -1.297721E-04  -4.645523E-08  -2.251402E-08  -6.473973E-07
  a    12   4.582312E-08   1.182337E-04   2.948157E-08   8.820292E-03  -5.141570E-09   9.994843E-04   1.568768E-07   6.102564E-04
  a    13   1.321617E-07  -2.480224E-09  -1.297721E-04  -5.141570E-09   7.483523E-03  -8.143918E-08   3.877172E-08  -3.265257E-07
  a    14  -3.099888E-07  -2.141927E-04  -4.645523E-08   9.994843E-04  -8.143918E-08   8.058547E-03   1.155847E-06   3.668755E-03
  a    15  -7.887187E-03   7.572797E-07  -2.251402E-08   1.568768E-07   3.877172E-08   1.155847E-06   3.863180E-03   1.512608E-06
  a    16   2.615266E-06   2.352432E-03  -6.473973E-07   6.102564E-04  -3.265257E-07   3.668755E-03   1.512608E-06   8.435387E-03
  a    17   8.670989E-08  -5.266754E-07  -2.131933E-03  -2.108932E-07  -1.125072E-03  -8.629979E-07   1.835655E-08  -4.125780E-07
  a    18  -2.229350E-07   1.057413E-03   6.319534E-08  -8.599893E-04   7.427995E-08  -1.118161E-03   3.891761E-08   1.470601E-04
  a    19   1.367309E-02  -3.077552E-06   5.240823E-08  -3.244795E-06  -3.295814E-09  -3.733305E-06   1.429221E-03  -3.880672E-06
  a    20   1.310246E-05   3.010018E-03  -1.859407E-07   3.318166E-03  -7.983001E-08   3.758387E-03   2.408084E-06   3.460048E-03
  a    21   6.174583E-08   3.013990E-07   6.824798E-03   1.927845E-07  -3.104286E-03  -3.270726E-08  -6.094351E-08  -3.062043E-07
  a    22   2.563853E-07   5.517797E-03   1.009994E-08  -4.010956E-03   2.140507E-07  -2.458322E-03  -5.937417E-08  -2.657469E-04
  a    23   7.630664E-08  -1.680591E-08   4.218728E-03   1.773305E-07   6.580798E-03   1.478009E-07   3.980104E-08  -1.910315E-07
  a    24   2.693109E-07   2.585727E-03  -1.713193E-07   4.205425E-03   1.710647E-07  -5.218061E-03  -7.231906E-07  -2.371361E-03
  a    25  -2.023511E-08   2.983800E-03  -2.618537E-07   2.352807E-03   3.639362E-08   2.138006E-03   4.343563E-07   1.179710E-03
  a    26   4.555165E-07  -1.591862E-07  -1.014334E-03  -1.339491E-07   3.375616E-04  -8.179817E-08  -1.583182E-07  -2.289757E-07
  a    27  -3.443565E-03  -3.598033E-08  -8.945399E-08  -3.035014E-08   4.253915E-08  -2.352077E-08   1.072610E-03  -3.615312E-07
  a    28   5.863135E-08   2.790028E-08   4.500494E-07   4.808503E-08   2.380721E-09   1.191812E-08  -1.419859E-08   1.186357E-08
  a    29  -4.617736E-09   8.407862E-10   1.840968E-07  -4.649689E-08   1.380912E-08  -1.421833E-08  -8.516389E-10   1.117632E-08
  a    30   3.105203E-07   1.494056E-06  -9.316769E-09   5.982306E-04   6.394169E-11   1.338006E-04  -4.317617E-08  -6.888407E-06
  a    31   1.426077E-09   2.464993E-08   8.638099E-04   1.149582E-08  -5.285881E-05  -7.392878E-09   3.541861E-08   3.864570E-07
  a    32  -5.083171E-08  -3.357425E-04  -3.386148E-09  -2.586445E-04   4.031847E-09   2.751846E-05  -7.952170E-07  -2.428256E-03
  a    33   3.906537E-07   3.582678E-04   4.152513E-11  -4.696673E-04  -2.597722E-08   2.087324E-05   1.448010E-07   4.605732E-04
  a    34   1.223910E-03  -6.251834E-08   1.607554E-08   7.262480E-08  -1.212120E-08   7.229094E-10   1.267741E-05  -2.870634E-08
  a    35  -1.058650E-07  -2.203106E-04  -1.840640E-08   3.628150E-05   7.731549E-08   4.874837E-04  -3.628118E-08  -5.427400E-05
  a    36   2.003462E-08   1.916916E-08  -2.038714E-04  -2.985630E-08   6.126730E-04  -6.849710E-08   8.783606E-09   5.251299E-08

               a    17        a    18        a    19        a    20        a    21        a    22        a    23        a    24
  a     3  -7.254955E-07   3.693223E-03  -1.117464E-06   1.164423E-03  -2.128754E-09   8.559871E-04  -2.784919E-09   1.005996E-03
  a     4   1.489058E-06  -3.527754E-03   2.239116E-06  -2.302006E-03  -1.610617E-09  -2.764242E-04  -9.244833E-08   6.477784E-03
  a     5  -7.723118E-03   1.708860E-07  -4.879100E-07  -1.206437E-07  -5.957158E-03  -1.707495E-07  -9.874412E-03  -2.022476E-07
  a     6   2.134945E-06   2.410525E-03   5.536188E-06  -5.953564E-03   3.311950E-07   1.077472E-02  -2.062642E-07   2.504141E-03
  a     7  -1.205870E-02   4.436648E-07  -1.778972E-06   1.122327E-07   1.073101E-02  -2.429336E-07  -7.628771E-03  -3.120514E-07
  a     8  -6.831623E-07  -1.616677E-07  -6.849752E-03  -6.805583E-06   9.724334E-07  -3.314800E-07   4.175091E-07   1.482172E-07
  a     9   8.670989E-08  -2.229350E-07   1.367309E-02   1.310246E-05   6.174583E-08   2.563853E-07   7.630664E-08   2.693109E-07
  a    10  -5.266754E-07   1.057413E-03  -3.077552E-06   3.010018E-03   3.013990E-07   5.517797E-03  -1.680591E-08   2.585727E-03
  a    11  -2.131933E-03   6.319534E-08   5.240823E-08  -1.859407E-07   6.824798E-03   1.009994E-08   4.218728E-03  -1.713193E-07
  a    12  -2.108932E-07  -8.599893E-04  -3.244795E-06   3.318166E-03   1.927845E-07  -4.010956E-03   1.773305E-07   4.205425E-03
  a    13  -1.125072E-03   7.427995E-08  -3.295814E-09  -7.983001E-08  -3.104286E-03   2.140507E-07   6.580798E-03   1.710647E-07
  a    14  -8.629979E-07  -1.118161E-03  -3.733305E-06   3.758387E-03  -3.270726E-08  -2.458322E-03   1.478009E-07  -5.218061E-03
  a    15   1.835655E-08   3.891761E-08   1.429221E-03   2.408084E-06  -6.094351E-08  -5.937417E-08   3.980104E-08  -7.231906E-07
  a    16  -4.125780E-07   1.470601E-04  -3.880672E-06   3.460048E-03  -3.062043E-07  -2.657469E-04  -1.910315E-07  -2.371361E-03
  a    17   6.403829E-03  -8.082779E-08   6.354404E-08  -7.809523E-07  -1.303621E-03   1.078096E-07  -1.046081E-03   5.316675E-07
  a    18  -8.082779E-08   4.582446E-03   1.040942E-07  -1.768884E-04   6.118340E-08   1.879809E-03   1.949250E-09   5.070852E-04
  a    19   6.354404E-08   1.040942E-07   7.963864E-03   2.019585E-07  -2.501356E-08   5.736150E-07   6.410469E-08   5.485320E-07
  a    20  -7.809523E-07  -1.768884E-04   2.019585E-07   7.586834E-03   1.310941E-08  -5.736983E-04   3.583003E-08  -4.785185E-04
  a    21  -1.303621E-03   6.118340E-08  -2.501356E-08   1.310941E-08   6.442969E-03   2.297992E-08  -1.438432E-04   1.182400E-08
  a    22   1.078096E-07   1.879809E-03   5.736150E-07  -5.736983E-04   2.297992E-08   6.682386E-03   4.631258E-08   1.212946E-03
  a    23  -1.046081E-03   1.949250E-09   6.410469E-08   3.583003E-08  -1.438432E-04   4.631258E-08   9.139659E-03  -3.582586E-08
  a    24   5.316675E-07   5.070852E-04   5.485320E-07  -4.785185E-04   1.182400E-08   1.212946E-03  -3.582586E-08   7.693127E-03
  a    25  -3.039231E-07  -6.407464E-04  -1.807518E-06   1.830766E-03  -6.883125E-08   4.621275E-04   3.245304E-08   3.472789E-04
  a    26  -9.828522E-04   6.188909E-08   9.358568E-08  -4.476300E-08  -1.041227E-03  -4.664188E-09   5.000959E-04  -3.927334E-08
  a    27  -1.765781E-07   5.864976E-08  -7.577512E-04  -7.711921E-07  -1.198436E-07  -5.023312E-08   6.387422E-08  -3.320563E-08
  a    28   5.185112E-07  -3.449823E-08   4.005831E-09   4.374684E-08   4.504868E-07   1.039190E-08  -1.439781E-07   7.947494E-09
  a    29   4.144801E-07  -1.069910E-07  -1.708534E-09  -1.171133E-07   1.314833E-07   6.605095E-08   1.168391E-07  -9.627485E-09
  a    30   5.604306E-08   1.140139E-03  -9.888081E-07   1.170044E-03  -1.315104E-08  -5.402703E-04   2.555594E-08   1.327406E-04
  a    31   2.110361E-03  -5.627365E-08   4.037962E-09  -9.363623E-09   4.522353E-04   3.059238E-08   7.094900E-04   4.146098E-09
  a    32   4.603633E-07  -4.904263E-05   4.766097E-07  -4.320955E-04  -1.514454E-08   2.200372E-04  -4.883731E-08   1.780103E-04
  a    33  -1.174652E-07  -2.257979E-04  -9.167201E-07   1.201146E-03   3.141561E-09   3.595913E-04  -1.762086E-08  -3.518005E-04
  a    34   2.898846E-08   7.109872E-09   7.752192E-04   4.534435E-07  -3.588874E-09  -8.421696E-08  -3.773830E-09   3.168538E-08
  a    35   1.914970E-08  -4.648957E-04  -9.659244E-07   9.359113E-04  -6.869952E-08  -4.381463E-04   6.249268E-08  -9.408372E-04
  a    36   1.710713E-04   5.912238E-08  -5.063357E-09  -1.145367E-07  -4.552488E-04   7.306625E-08   2.038031E-04   1.251277E-07

               a    25        a    26        a    27        a    28        a    29        a    30        a    31        a    32
  a     3   3.406720E-03  -1.283495E-07   2.523552E-07   1.244972E-08  -6.612335E-07   6.141487E-03  -1.933677E-07  -5.004970E-03
  a     4   6.701952E-03  -2.859841E-07  -5.698950E-07   1.928587E-07   2.399148E-07  -2.312562E-03   2.296748E-07   4.231009E-03
  a     5   3.364064E-07   7.288095E-03   8.092757E-07  -3.904148E-06  -1.366678E-06  -1.026256E-07  -6.903485E-04   3.541829E-08
  a     6  -5.780643E-05  -1.755206E-07  -3.299425E-07   7.444577E-08   1.019236E-07  -3.503705E-04   4.335409E-07   9.857056E-03
  a     7   7.438911E-08   1.706030E-03   3.754759E-07  -1.347252E-06  -4.326835E-07  -1.046908E-08  -2.536684E-04   5.784195E-08
  a     8   2.233545E-07  -1.689052E-07   3.478473E-03   9.396832E-08  -2.466417E-08   1.739699E-07   2.124918E-08  -2.533116E-07
  a     9  -2.023511E-08   4.555165E-07  -3.443565E-03   5.863135E-08  -4.617736E-09   3.105203E-07   1.426077E-09  -5.083171E-08
  a    10   2.983800E-03  -1.591862E-07  -3.598033E-08   2.790028E-08   8.407862E-10   1.494056E-06   2.464993E-08  -3.357425E-04
  a    11  -2.618537E-07  -1.014334E-03  -8.945399E-08   4.500494E-07   1.840968E-07  -9.316769E-09   8.638099E-04  -3.386148E-09
  a    12   2.352807E-03  -1.339491E-07  -3.035014E-08   4.808503E-08  -4.649689E-08   5.982306E-04   1.149582E-08  -2.586445E-04
  a    13   3.639362E-08   3.375616E-04   4.253915E-08   2.380721E-09   1.380912E-08   6.394169E-11  -5.285881E-05   4.031847E-09
  a    14   2.138006E-03  -8.179817E-08  -2.352077E-08   1.191812E-08  -1.421833E-08   1.338006E-04  -7.392878E-09   2.751846E-05
  a    15   4.343563E-07  -1.583182E-07   1.072610E-03  -1.419859E-08  -8.516389E-10  -4.317617E-08   3.541861E-08  -7.952170E-07
  a    16   1.179710E-03  -2.289757E-07  -3.615312E-07   1.186357E-08   1.117632E-08  -6.888407E-06   3.864570E-07  -2.428256E-03
  a    17  -3.039231E-07  -9.828522E-04  -1.765781E-07   5.185112E-07   4.144801E-07   5.604306E-08   2.110361E-03   4.603633E-07
  a    18  -6.407464E-04   6.188909E-08   5.864976E-08  -3.449823E-08  -1.069910E-07   1.140139E-03  -5.627365E-08  -4.904263E-05
  a    19  -1.807518E-06   9.358568E-08  -7.577512E-04   4.005831E-09  -1.708534E-09  -9.888081E-07   4.037962E-09   4.766097E-07
  a    20   1.830766E-03  -4.476300E-08  -7.711921E-07   4.374684E-08  -1.171133E-07   1.170044E-03  -9.363623E-09  -4.320955E-04
  a    21  -6.883125E-08  -1.041227E-03  -1.198436E-07   4.504868E-07   1.314833E-07  -1.315104E-08   4.522353E-04  -1.514454E-08
  a    22   4.621275E-04  -4.664188E-09  -5.023312E-08   1.039190E-08   6.605095E-08  -5.402703E-04   3.059238E-08   2.200372E-04
  a    23   3.245304E-08   5.000959E-04   6.387422E-08  -1.439781E-07   1.168391E-07   2.555594E-08   7.094900E-04  -4.883731E-08
  a    24   3.472789E-04  -3.927334E-08  -3.320563E-08   7.947494E-09  -9.627485E-09   1.327406E-04   4.146098E-09   1.780103E-04
  a    25   3.437274E-03  -5.694775E-08   6.873461E-08  -2.172550E-08   4.127784E-08  -3.768346E-04  -4.801578E-08  -2.054768E-04
  a    26  -5.694775E-08   1.896301E-03  -3.735311E-07   9.948100E-07   1.852054E-07   2.053985E-08  -8.475745E-04   2.210415E-08
  a    27   6.873461E-08  -3.735311E-07   4.465614E-03  -1.571660E-09   5.310285E-09   3.774006E-09  -1.153750E-07   1.181506E-08
  a    28  -2.172550E-08   9.948100E-07  -1.571660E-09   3.935825E-03   6.062261E-04   5.271837E-08   2.704346E-07   2.322011E-09
  a    29   4.127784E-08   1.852054E-07   5.310285E-09   6.062261E-04   2.607866E-03   3.932329E-08   1.004991E-07  -1.862453E-08
  a    30  -3.768346E-04   2.053985E-08   3.774006E-09   5.271837E-08   3.932329E-08   2.002487E-03   2.261965E-08   1.227867E-04
  a    31  -4.801578E-08  -8.475745E-04  -1.153750E-07   2.704346E-07   1.004991E-07   2.261965E-08   2.438619E-03  -1.194613E-08
  a    32  -2.054768E-04   2.210415E-08   1.181506E-08   2.322011E-09  -1.862453E-08   1.227867E-04  -1.194613E-08   2.174504E-03
  a    33  -3.723273E-04   4.859548E-08   1.200324E-07  -5.824392E-09  -3.969631E-08   3.467135E-04  -1.190237E-08  -2.578034E-04
  a    34   1.152224E-07  -7.165781E-08   4.889434E-04  -6.340617E-09   3.001874E-09  -7.728426E-08   5.663348E-09   2.860913E-08
  a    35   6.378614E-04  -1.039972E-07  -9.147256E-09   2.286217E-09   6.898227E-09  -8.573069E-05   2.962688E-08   9.532890E-05
  a    36  -1.241357E-07  -6.242976E-04  -8.013510E-08   2.453202E-07   5.431370E-08   1.543169E-08   1.625343E-04  -2.243684E-08

               a    33        a    34        a    35        a    36
  a     3   1.821964E-03  -3.355745E-07   2.360776E-04  -1.599007E-08
  a     4   2.890489E-03  -8.413513E-07   8.742271E-04  -5.335765E-08
  a     5  -1.269536E-07   6.641386E-08   1.176131E-07   1.124177E-03
  a     6  -1.366678E-03  -2.902631E-07   7.125759E-04  -8.230607E-08
  a     7  -9.589028E-08  -1.895181E-07   1.058760E-07   1.100026E-03
  a     8   7.590741E-07   5.552790E-03   1.739631E-07   1.551432E-07
  a     9   3.906537E-07   1.223910E-03  -1.058650E-07   2.003462E-08
  a    10   3.582678E-04  -6.251834E-08  -2.203106E-04   1.916916E-08
  a    11   4.152513E-11   1.607554E-08  -1.840640E-08  -2.038714E-04
  a    12  -4.696673E-04   7.262480E-08   3.628150E-05  -2.985630E-08
  a    13  -2.597722E-08  -1.212120E-08   7.731549E-08   6.126730E-04
  a    14   2.087324E-05   7.229094E-10   4.874837E-04  -6.849710E-08
  a    15   1.448010E-07   1.267741E-05  -3.628118E-08   8.783606E-09
  a    16   4.605732E-04  -2.870634E-08  -5.427400E-05   5.251299E-08
  a    17  -1.174652E-07   2.898846E-08   1.914970E-08   1.710713E-04
  a    18  -2.257979E-04   7.109872E-09  -4.648957E-04   5.912238E-08
  a    19  -9.167201E-07   7.752192E-04  -9.659244E-07  -5.063357E-09
  a    20   1.201146E-03   4.534435E-07   9.359113E-04  -1.145367E-07
  a    21   3.141561E-09  -3.588874E-09  -6.869952E-08  -4.552488E-04
  a    22   3.595913E-04  -8.421696E-08  -4.381463E-04   7.306625E-08
  a    23  -1.762086E-08  -3.773830E-09   6.249268E-08   2.038031E-04
  a    24  -3.518005E-04   3.168538E-08  -9.408372E-04   1.251277E-07
  a    25  -3.723273E-04   1.152224E-07   6.378614E-04  -1.241357E-07
  a    26   4.859548E-08  -7.165781E-08  -1.039972E-07  -6.242976E-04
  a    27   1.200324E-07   4.889434E-04  -9.147256E-09  -8.013510E-08
  a    28  -5.824392E-09  -6.340617E-09   2.286217E-09   2.453202E-07
  a    29  -3.969631E-08   3.001874E-09   6.898227E-09   5.431370E-08
  a    30   3.467135E-04  -7.728426E-08  -8.573069E-05   1.543169E-08
  a    31  -1.190237E-08   5.663348E-09   2.962688E-08   1.625343E-04
  a    32  -2.578034E-04   2.860913E-08   9.532890E-05  -2.243684E-08
  a    33   1.488571E-03  -4.473168E-08   2.746668E-04  -3.467858E-08
  a    34  -4.473168E-08   1.245409E-03  -6.162677E-08  -1.365291E-09
  a    35   2.746668E-04  -6.162677E-08   1.330765E-03  -1.578858E-08
  a    36  -3.467858E-08  -1.365291E-09  -1.578858E-08   1.170578E-03

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98374602     1.97674015     1.97513790     1.97142136     1.96851533     1.92550515
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.06407142     0.01939344     0.01784900     0.01655840     0.01511318     0.01441246     0.00721783     0.00638874
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00594467     0.00507869     0.00417094     0.00393857     0.00374268     0.00237275     0.00219499     0.00177509
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00177503     0.00122674     0.00114741     0.00100321     0.00078188     0.00076173     0.00073293     0.00031236
              MO    33       MO    34       MO    35       MO    36
  occ(*)=     0.00030541     0.00026658     0.00025347     0.00014451


 total number of electrons =   16.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
    N1_ s       1.998472  -0.000030   1.268730   0.292748   0.000000   0.000000
    N1_ p       0.000000  -0.000292  -0.005749   0.188923   0.430987   0.936459
    N1_ d       0.000000  -0.000105   0.003777   0.002703   0.000438   0.017377
    C1_ s       0.000539   1.999452   0.463001   0.791269   0.000000   0.000000
    C1_ p       0.000476   0.000003  -0.003788   0.172430   0.946739   0.262439
    C1_ d      -0.000081   0.000001   0.005694   0.005714   0.001017   0.017846
    H1_ s       0.000301   0.000003   0.080588   0.063329   0.051147   0.244845
    H2_ s       0.000301   0.000003   0.080588   0.063326   0.051146   0.244849
    H3_ s      -0.000004   0.000483   0.045451   0.198153   0.246831   0.123800
    H4_ s      -0.000004   0.000483   0.045453   0.198145   0.246832   0.123808
 
   ao class       7a         8a         9a        10a        11a        12a  
    N1_ s       0.015048   0.000000   0.000000   0.002714   0.000000   0.004095
    N1_ p       0.994746   1.333444   0.022696   0.005070   0.008581   0.003653
    N1_ d       0.005011   0.004910   0.000823   0.000821   0.000927   0.000647
    C1_ s       0.014564   0.000000   0.000000   0.003841   0.000000   0.000183
    C1_ p       0.656507   0.574649   0.040613   0.006253   0.000040   0.000192
    C1_ d       0.012930   0.012502  -0.000061   0.000273   0.000377   0.000026
    H1_ s       0.072158  -0.000000   0.000000   0.000023   0.003516   0.003694
    H2_ s       0.072161   0.000000   0.000000   0.000023   0.003517   0.003694
    H3_ s       0.062695   0.000000   0.000000   0.000188   0.000445   0.000187
    H4_ s       0.062694   0.000000   0.000000   0.000188   0.000445   0.000187
 
   ao class      13a        14a        15a        16a        17a        18a  
    N1_ s       0.000000   0.000244   0.000000   0.000000   0.000625   0.000607
    N1_ p       0.001413   0.000228   0.001477   0.003743   0.001067   0.000342
    N1_ d       0.000340   0.000035   0.001946   0.001513   0.002378   0.002276
    C1_ s       0.000000   0.003343   0.000000   0.000000   0.000335   0.000693
    C1_ p       0.007482   0.003108   0.002056   0.000882   0.000661   0.000080
    C1_ d       0.000195   0.000140   0.001675   0.000251   0.000872   0.000864
    H1_ s       0.000077   0.000176  -0.000014  -0.000000  -0.000003   0.000048
    H2_ s       0.000077   0.000176  -0.000014  -0.000000  -0.000003   0.000048
    H3_ s       0.002765   0.003482   0.000046   0.000000   0.000005   0.000061
    H4_ s       0.002765   0.003482   0.000046  -0.000000   0.000005   0.000061
 
   ao class      19a        20a        21a        22a        23a        24a  
    N1_ s       0.000000   0.000000   0.000065  -0.000000   0.000000   0.000076
    N1_ p       0.000000   0.000897   0.000406  -0.000000   0.000085   0.000095
    N1_ d       0.003472   0.001494   0.000767   0.000397   0.000882   0.000358
    C1_ s       0.000000   0.000000   0.000381   0.000000   0.000000   0.000278
    C1_ p       0.000000   0.000188   0.000721   0.000000   0.000440   0.000085
    C1_ d       0.000699   0.001360   0.001054   0.001975   0.000347   0.000723
    H1_ s       0.000000   0.000000  -0.000011  -0.000000   0.000037  -0.000004
    H2_ s      -0.000000  -0.000000  -0.000011   0.000000   0.000037  -0.000004
    H3_ s       0.000000   0.000000   0.000185   0.000000   0.000184   0.000083
    H4_ s      -0.000000   0.000000   0.000185   0.000000   0.000184   0.000083
 
   ao class      25a        26a        27a        28a        29a        30a  
    N1_ s       0.000000   0.000373   0.000000   0.000046   0.000000   0.000022
    N1_ p       0.000056  -0.000002   0.000206   0.000170   0.000084   0.000178
    N1_ d       0.000344   0.000184   0.000041   0.000234   0.000137   0.000163
    C1_ s       0.000001  -0.000052  -0.000000   0.000020   0.000000   0.000012
    C1_ p       0.001062  -0.000016   0.000091   0.000181   0.000220   0.000008
    C1_ d       0.000311   0.000280   0.000480   0.000173   0.000341   0.000183
    H1_ s      -0.000000   0.000180   0.000090   0.000016   0.000000   0.000047
    H2_ s      -0.000000   0.000180   0.000090   0.000016   0.000000   0.000047
    H3_ s       0.000000   0.000049   0.000075   0.000073   0.000000   0.000051
    H4_ s       0.000000   0.000049   0.000075   0.000073   0.000000   0.000051
 
   ao class      31a        32a        33a        34a        35a        36a  
    N1_ s       0.000000   0.000000   0.000022   0.000000   0.000002   0.000008
    N1_ p       0.000072   0.000029  -0.000001   0.000032   0.000017   0.000004
    N1_ d       0.000146   0.000000   0.000004   0.000002   0.000006   0.000003
    C1_ s      -0.000000   0.000000   0.000024   0.000000  -0.000001   0.000057
    C1_ p      -0.000001   0.000000   0.000087   0.000087   0.000010   0.000023
    C1_ d       0.000090   0.000004   0.000028   0.000003   0.000015   0.000024
    H1_ s       0.000094   0.000122   0.000028   0.000009   0.000084   0.000000
    H2_ s       0.000094   0.000122   0.000028   0.000009   0.000084   0.000000
    H3_ s       0.000120   0.000017   0.000043   0.000062   0.000018   0.000012
    H4_ s       0.000120   0.000017   0.000043   0.000062   0.000018   0.000012


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.583869   3.277942   0.520582   0.520584   0.685560   0.685562
      p         3.929117   2.674008   0.000000   0.000000   0.000000   0.000000
      d         0.054451   0.068326   0.000000   0.000000   0.000000   0.000000
    total       7.567436   6.020276   0.520582   0.520584   0.685560   0.685562
 

 Total number of electrons:   16.00000000

 item #                     2 suffix=:.drt1.state2:
 read_civout: repnuc=  -88.6565513780480     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 96% root-following 0
 MR-CISD energy:   -94.66619491    -6.00964353
 residuum:     0.00037264
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.34532008    -5.68876870
 residuum:     0.00065006
 deltae:     0.00000090
 apxde:     0.00000021

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821    -0.01624692
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894     0.00491632
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907     0.10910721

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821     0.00000000
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894     0.00000000
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 computing MRCISD density
 densi: densityinfo%a4den=   1.00000000000000     
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    2
--------------------------------------------------------------------------------
================================================================================
   DYZ=     215  DYX=     335  DYW=     401
   D0Z=     147  D0Y=     776  D0X=     489  D0W=     695
  DDZI=     109 DDYI=     378 DDXI=     272 DDWI=     359
  DDZE=       0 DDYE=      92 DDXE=      71 DDWE=      98
================================================================================
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

          modens reordered block   1

               a     1        a     2        a     3        a     4        a     5        a     6        a     7        a     8
  a     1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     3    0.00000        0.00000        1.98322       1.002739E-03   7.560863E-08   2.301568E-03   9.739854E-08   3.570807E-07
  a     4    0.00000        0.00000       1.002739E-03    1.97206      -1.069016E-09   1.394967E-03  -2.972909E-07   1.061914E-07
  a     5    0.00000        0.00000       7.560863E-08  -1.069016E-09    1.95906       1.413092E-07  -5.956572E-02  -2.664365E-07
  a     6    0.00000        0.00000       2.301568E-03   1.394967E-03   1.413092E-07    1.96546       1.239629E-06  -4.704762E-08
  a     7    0.00000        0.00000       9.739854E-08  -2.972909E-07  -5.956572E-02   1.239629E-06    1.00842      -1.160431E-04
  a     8    0.00000        0.00000       3.570807E-07   1.061914E-07  -2.664365E-07  -4.704762E-08  -1.160431E-04    1.90338    
  a     9    0.00000        0.00000      -2.052204E-06  -4.240549E-07  -1.618405E-05   2.270850E-08  -7.076204E-05   0.226463    
  a    10    0.00000        0.00000      -9.096724E-03  -6.712773E-03   4.956633E-07  -4.587551E-03  -7.522877E-07  -6.614130E-07
  a    11    0.00000        0.00000       4.314946E-07   5.216386E-07   7.168408E-03  -1.355336E-07  -1.685473E-02   8.211097E-05
  a    12    0.00000        0.00000      -4.246351E-03  -6.264283E-03   1.942142E-07   6.492710E-03  -5.233688E-07  -1.279531E-07
  a    13    0.00000        0.00000      -7.842506E-08   1.349685E-07   9.166638E-03  -1.797849E-07   1.870408E-02   3.167077E-07
  a    14    0.00000        0.00000      -7.021937E-03  -1.884400E-03   1.510874E-07   1.630932E-02  -3.520259E-07   5.162317E-07
  a    15    0.00000        0.00000       8.406901E-07  -2.557154E-06  -3.879047E-07  -2.295541E-06   7.619098E-06  -1.811162E-02
  a    16    0.00000        0.00000       2.885562E-03  -6.663491E-03   1.221025E-06  -5.835182E-03   3.024809E-06   6.931532E-06
  a    17    0.00000        0.00000      -6.520527E-07   1.477633E-06   5.735992E-03   1.114131E-06   1.428745E-02  -1.288148E-05
  a    18    0.00000        0.00000       3.427560E-03  -1.994034E-05  -1.492973E-07  -9.003126E-03  -6.658007E-07  -1.373182E-07
  a    19    0.00000        0.00000       3.884409E-06   3.129226E-06  -9.153721E-07  -5.202905E-06   4.612780E-06  -2.472729E-03
  a    20    0.00000        0.00000      -4.507929E-03  -3.242739E-03   1.790160E-07   5.259996E-03  -2.720821E-07  -3.226365E-06
  a    21    0.00000        0.00000       3.490607E-08   8.259420E-08   1.401914E-04  -2.034256E-07  -1.930310E-02   1.840571E-05
  a    22    0.00000        0.00000      -1.509633E-03   3.277089E-03   2.761206E-08  -4.596668E-03   4.469015E-07  -3.541954E-07
  a    23    0.00000        0.00000      -4.610228E-08  -8.335385E-08   1.512165E-03   2.048344E-07   1.603725E-02   9.199902E-06
  a    24    0.00000        0.00000       7.133592E-04  -2.029562E-03   6.382871E-09   7.695614E-05   2.542089E-07  -3.024204E-07
  a    25    0.00000        0.00000      -1.042177E-04   2.764664E-03   1.728519E-07   2.025384E-03   1.038849E-06   6.709282E-07
  a    26    0.00000        0.00000       1.214309E-08  -9.787806E-08   2.655191E-03  -2.515736E-07   2.010840E-02  -4.355759E-06
  a    27    0.00000        0.00000       1.268697E-07   2.671290E-07   1.032815E-06   2.965963E-07  -5.008763E-06   1.871048E-02
  a    28    0.00000        0.00000       6.081769E-07  -2.805772E-07  -1.070537E-06   1.006349E-06  -7.860756E-06  -3.111037E-07
  a    29    0.00000        0.00000       1.692304E-06   8.682604E-07   8.414728E-07  -1.183471E-06  -1.895040E-06  -2.246780E-08
  a    30    0.00000        0.00000      -4.274733E-04  -1.509494E-03   1.191945E-07   7.244504E-03   9.931614E-10   1.581444E-07
  a    31    0.00000        0.00000      -2.742795E-07   2.663958E-07   6.012079E-03   5.624708E-07  -5.520713E-03   3.316312E-06
  a    32    0.00000        0.00000      -5.750289E-03   3.549434E-03  -1.672781E-07   1.379574E-02   2.361998E-07  -3.799566E-07
  a    33    0.00000        0.00000       7.400501E-04   4.504633E-03  -2.125554E-07  -3.758154E-03  -1.816606E-07   1.650936E-06
  a    34    0.00000        0.00000      -5.924247E-08  -9.129967E-07  -1.765362E-07   5.553232E-07  -8.516277E-07   5.506909E-03
  a    35    0.00000        0.00000       3.486177E-05   1.987005E-03   2.835815E-07   4.984705E-04   6.347045E-07  -1.446255E-07
  a    36    0.00000        0.00000       7.064396E-09  -1.438909E-07   2.497168E-03  -1.177273E-07   4.459501E-03   1.079015E-06

               a     9        a    10        a    11        a    12        a    13        a    14        a    15        a    16
  a     3  -2.052204E-06  -9.096724E-03   4.314946E-07  -4.246351E-03  -7.842506E-08  -7.021937E-03   8.406901E-07   2.885562E-03
  a     4  -4.240549E-07  -6.712773E-03   5.216386E-07  -6.264283E-03   1.349685E-07  -1.884400E-03  -2.557154E-06  -6.663491E-03
  a     5  -1.618405E-05   4.956633E-07   7.168408E-03   1.942142E-07   9.166638E-03   1.510874E-07  -3.879047E-07   1.221025E-06
  a     6   2.270850E-08  -4.587551E-03  -1.355336E-07   6.492710E-03  -1.797849E-07   1.630932E-02  -2.295541E-06  -5.835182E-03
  a     7  -7.076204E-05  -7.522877E-07  -1.685473E-02  -5.233688E-07   1.870408E-02  -3.520259E-07   7.619098E-06   3.024809E-06
  a     8   0.226463      -6.614130E-07   8.211097E-05  -1.279531E-07   3.167077E-07   5.162317E-07  -1.811162E-02   6.931532E-06
  a     9    1.06038       2.641139E-07  -1.202987E-05   7.973289E-07  -1.218740E-05  -8.263493E-07  -1.438382E-02   4.241036E-06
  a    10   2.641139E-07   9.111071E-03   2.458993E-08  -3.223011E-04   2.359106E-08   2.716228E-05   7.616036E-07   2.217869E-03
  a    11  -1.202987E-05   2.458993E-08   9.346098E-03  -3.637371E-08   4.747809E-04  -4.316966E-08   1.629000E-06  -8.965952E-07
  a    12   7.973289E-07  -3.223011E-04  -3.637371E-08   1.008897E-02   4.259972E-08   9.196870E-04   3.010502E-07   1.025598E-03
  a    13  -1.218740E-05   2.359106E-08   4.747809E-04   4.259972E-08   5.814372E-03  -8.832194E-08  -3.564900E-07   7.049937E-09
  a    14  -8.263493E-07   2.716228E-05  -4.316966E-08   9.196870E-04  -8.832194E-08   7.972977E-03   1.046400E-06   3.324171E-03
  a    15  -1.438382E-02   7.616036E-07   1.629000E-06   3.010502E-07  -3.564900E-07   1.046400E-06   8.316169E-03  -1.780055E-09
  a    16   4.241036E-06   2.217869E-03  -8.965952E-07   1.025598E-03   7.049937E-09   3.324171E-03  -1.780055E-09   8.356744E-03
  a    17   1.095626E-05  -5.169480E-07  -3.161572E-03  -4.169577E-07   2.818339E-04  -7.640034E-07   1.391266E-07  -5.149250E-07
  a    18  -6.895641E-07   2.120058E-03   1.597107E-07  -3.505116E-03  -1.011918E-08  -7.415956E-04   1.030197E-08  -5.889810E-05
  a    19   2.956098E-02  -3.127780E-06   8.667920E-07  -3.999352E-06  -2.698351E-08  -3.628687E-06   8.667621E-04  -3.657761E-06
  a    20   2.827639E-05   3.042044E-03  -2.188772E-07   4.095089E-03  -5.111539E-08   3.683078E-03   1.967726E-06   3.456468E-03
  a    21   1.950110E-06   2.577244E-07   5.689368E-03   1.205245E-07  -1.838450E-03  -2.027987E-08   5.870687E-07  -5.968009E-07
  a    22   2.537571E-07   5.810148E-03   2.990130E-08  -4.717327E-03   1.554544E-07  -1.981727E-03  -6.803271E-08  -3.870937E-04
  a    23   5.645281E-06  -2.251418E-08   4.107195E-03   2.106919E-07   5.243140E-03   1.334665E-07   1.931775E-07  -5.124824E-08
  a    24   8.064343E-07   2.077560E-03  -1.932015E-07   5.105074E-03   1.614215E-07  -4.932377E-03  -5.258388E-07  -1.928068E-03
  a    25  -1.902620E-07   2.994923E-03  -2.173319E-07   1.625617E-03   1.774556E-08   2.047850E-03   2.921829E-07   9.968503E-04
  a    26  -4.610681E-06  -1.414197E-07  -6.588270E-04  -6.401151E-08   1.870822E-04  -8.002945E-08  -5.928642E-08  -1.110739E-07
  a    27  -4.430782E-03  -6.665488E-08   6.700205E-07  -4.586145E-08   4.291287E-07   1.664841E-09  -6.261744E-04   2.115159E-07
  a    28   2.006527E-08  -7.040723E-07   2.692533E-07   4.267787E-07   1.103501E-07  -3.632065E-08   2.686502E-09  -1.193227E-07
  a    29   2.479427E-08  -7.952885E-08   4.279304E-08  -5.795160E-07   1.066843E-07   2.221165E-07  -1.705209E-09   1.453666E-07
  a    30   1.683269E-07   6.511346E-04   3.737290E-09  -6.010387E-04  -2.551589E-09   1.195990E-05  -8.305739E-08  -1.920600E-04
  a    31   4.691285E-06  -1.429666E-08   1.407370E-04   1.438346E-08   3.875359E-04   1.121542E-08   3.807919E-07   2.788328E-07
  a    32   1.014278E-07  -3.082594E-04   2.036564E-08  -1.839305E-04  -1.534348E-08   3.212746E-04  -7.352213E-07  -2.240131E-03
  a    33  -6.193541E-07   2.257324E-04  -4.051278E-08   8.273889E-04  -1.234160E-08   6.440214E-05   1.773594E-07   5.488014E-04
  a    34  -4.899101E-03  -6.759656E-08   3.377459E-07  -1.686753E-07   8.333752E-08  -8.770237E-09  -1.698625E-04   1.970780E-08
  a    35  -2.956392E-07  -2.808890E-04  -1.571000E-08  -5.479109E-05   7.050199E-08   3.449333E-04  -8.201378E-08  -1.487947E-04
  a    36  -2.245978E-06   2.579793E-08  -2.295422E-04  -1.062622E-08   5.660066E-04  -5.028012E-08  -1.654334E-10   9.500643E-08

               a    17        a    18        a    19        a    20        a    21        a    22        a    23        a    24
  a     3  -6.520527E-07   3.427560E-03   3.884409E-06  -4.507929E-03   3.490607E-08  -1.509633E-03  -4.610228E-08   7.133592E-04
  a     4   1.477633E-06  -1.994034E-05   3.129226E-06  -3.242739E-03   8.259420E-08   3.277089E-03  -8.335385E-08  -2.029562E-03
  a     5   5.735992E-03  -1.492973E-07  -9.153721E-07   1.790160E-07   1.401914E-04   2.761206E-08   1.512165E-03   6.382871E-09
  a     6   1.114131E-06  -9.003126E-03  -5.202905E-06   5.259996E-03  -2.034256E-07  -4.596668E-03   2.048344E-07   7.695614E-05
  a     7   1.428745E-02  -6.658007E-07   4.612780E-06  -2.720821E-07  -1.930310E-02   4.469015E-07   1.603725E-02   2.542089E-07
  a     8  -1.288148E-05  -1.373182E-07  -2.472729E-03  -3.226365E-06   1.840571E-05  -3.541954E-07   9.199902E-06  -3.024204E-07
  a     9   1.095626E-05  -6.895641E-07   2.956098E-02   2.827639E-05   1.950110E-06   2.537571E-07   5.645281E-06   8.064343E-07
  a    10  -5.169480E-07   2.120058E-03  -3.127780E-06   3.042044E-03   2.577244E-07   5.810148E-03  -2.251418E-08   2.077560E-03
  a    11  -3.161572E-03   1.597107E-07   8.667920E-07  -2.188772E-07   5.689368E-03   2.990130E-08   4.107195E-03  -1.932015E-07
  a    12  -4.169577E-07  -3.505116E-03  -3.999352E-06   4.095089E-03   1.205245E-07  -4.717327E-03   2.106919E-07   5.105074E-03
  a    13   2.818339E-04  -1.011918E-08  -2.698351E-08  -5.111539E-08  -1.838450E-03   1.554544E-07   5.243140E-03   1.614215E-07
  a    14  -7.640034E-07  -7.415956E-04  -3.628687E-06   3.683078E-03  -2.027987E-08  -1.981727E-03   1.334665E-07  -4.932377E-03
  a    15   1.391266E-07   1.030197E-08   8.667621E-04   1.967726E-06   5.870687E-07  -6.803271E-08   1.931775E-07  -5.258388E-07
  a    16  -5.149250E-07  -5.889810E-05  -3.657761E-06   3.456468E-03  -5.968009E-07  -3.870937E-04  -5.124824E-08  -1.928068E-03
  a    17   5.911346E-03   5.471300E-08   1.177651E-06  -8.153572E-07  -2.541752E-03   1.903174E-07  -4.544524E-04   3.972713E-07
  a    18   5.471300E-08   7.504015E-03   8.924927E-07  -1.059691E-03   1.371836E-07   3.208718E-03  -6.591312E-08  -9.062028E-04
  a    19   1.177651E-06   8.924927E-07   1.048635E-02   2.005787E-06   2.918177E-07   6.766782E-07   3.925095E-07  -8.575634E-08
  a    20  -8.153572E-07  -1.059691E-03   2.005787E-06   8.201963E-03  -1.362726E-08  -6.671238E-04   5.680082E-08   1.689522E-04
  a    21  -2.541752E-03   1.371836E-07   2.918177E-07  -1.362726E-08   4.834555E-03   5.412277E-08   4.439949E-04  -1.427710E-09
  a    22   1.903174E-07   3.208718E-03   6.766782E-07  -6.671238E-04   5.412277E-08   6.803613E-03  -4.886972E-09   4.363982E-04
  a    23  -4.544524E-04  -6.591312E-08   3.925095E-07   5.680082E-08   4.439949E-04  -4.886972E-09   7.529949E-03  -2.553379E-08
  a    24   3.972713E-07  -9.062028E-04  -8.575634E-08   1.689522E-04  -1.427710E-09   4.363982E-04  -2.553379E-08   7.713110E-03
  a    25  -2.254084E-07  -1.147065E-04  -1.297095E-06   1.316700E-03  -3.738351E-08   7.357336E-04   2.019067E-08   2.838027E-05
  a    26  -5.162576E-04   1.074772E-08  -2.514662E-08  -4.284207E-09  -6.382913E-04  -2.441621E-08   5.627158E-04  -7.464951E-09
  a    27  -1.040873E-06   4.192162E-08  -1.260197E-03  -1.318050E-06  -2.511186E-07  -7.642009E-08   7.767162E-08  -6.324579E-08
  a    28   2.957290E-07  -2.290238E-07   7.133924E-09  -1.380004E-07   2.509080E-07  -3.387113E-07  -1.758497E-07   1.837410E-08
  a    29   2.782250E-07  -9.493213E-08   1.694190E-10  -1.829988E-07  -2.412779E-08   3.122906E-08   1.057780E-07  -2.914750E-07
  a    30   1.117266E-07   1.940352E-03  -6.727801E-07   8.509538E-04  -7.596592E-09   4.322085E-04  -4.431099E-09  -2.371015E-04
  a    31   1.599288E-03  -7.205009E-08   2.945184E-07   5.426192E-09  -3.047296E-04   1.450196E-08   7.017821E-04   1.055650E-08
  a    32   4.360837E-07  -2.997855E-04   2.550495E-07  -1.967791E-04   1.159011E-08   1.766564E-04  -4.037689E-08   8.567100E-05
  a    33  -1.572769E-07  -9.211895E-04  -1.494631E-06   1.807438E-03  -1.195548E-08  -2.411576E-04   1.740867E-09   3.145736E-04
  a    34  -1.615504E-07   1.271746E-07   6.719828E-04   2.200286E-07   2.936365E-08   5.640847E-09  -5.851384E-09  -9.843643E-08
  a    35   6.014801E-08  -3.688835E-04  -7.673650E-07   7.188188E-04  -5.708482E-08  -4.153080E-04   4.448254E-08  -8.938465E-04
  a    36   3.080533E-04   3.673748E-08   3.069218E-09  -8.348458E-08  -3.948962E-04   6.287394E-08   1.124515E-04   1.211818E-07

               a    25        a    26        a    27        a    28        a    29        a    30        a    31        a    32
  a     3  -1.042177E-04   1.214309E-08   1.268697E-07   6.081769E-07   1.692304E-06  -4.274733E-04  -2.742795E-07  -5.750289E-03
  a     4   2.764664E-03  -9.787806E-08   2.671290E-07  -2.805772E-07   8.682604E-07  -1.509494E-03   2.663958E-07   3.549434E-03
  a     5   1.728519E-07   2.655191E-03   1.032815E-06  -1.070537E-06   8.414728E-07   1.191945E-07   6.012079E-03  -1.672781E-07
  a     6   2.025384E-03  -2.515736E-07   2.965963E-07   1.006349E-06  -1.183471E-06   7.244504E-03   5.624708E-07   1.379574E-02
  a     7   1.038849E-06   2.010840E-02  -5.008763E-06  -7.860756E-06  -1.895040E-06   9.931614E-10  -5.520713E-03   2.361998E-07
  a     8   6.709282E-07  -4.355759E-06   1.871048E-02  -3.111037E-07  -2.246780E-08   1.581444E-07   3.316312E-06  -3.799566E-07
  a     9  -1.902620E-07  -4.610681E-06  -4.430782E-03   2.006527E-08   2.479427E-08   1.683269E-07   4.691285E-06   1.014278E-07
  a    10   2.994923E-03  -1.414197E-07  -6.665488E-08  -7.040723E-07  -7.952885E-08   6.511346E-04  -1.429666E-08  -3.082594E-04
  a    11  -2.173319E-07  -6.588270E-04   6.700205E-07   2.692533E-07   4.279304E-08   3.737290E-09   1.407370E-04   2.036564E-08
  a    12   1.625617E-03  -6.401151E-08  -4.586145E-08   4.267787E-07  -5.795160E-07  -6.010387E-04   1.438346E-08  -1.839305E-04
  a    13   1.774556E-08   1.870822E-04   4.291287E-07   1.103501E-07   1.066843E-07  -2.551589E-09   3.875359E-04  -1.534348E-08
  a    14   2.047850E-03  -8.002945E-08   1.664841E-09  -3.632065E-08   2.221165E-07   1.195990E-05   1.121542E-08   3.212746E-04
  a    15   2.921829E-07  -5.928642E-08  -6.261744E-04   2.686502E-09  -1.705209E-09  -8.305739E-08   3.807919E-07  -7.352213E-07
  a    16   9.968503E-04  -1.110739E-07   2.115159E-07  -1.193227E-07   1.453666E-07  -1.920600E-04   2.788328E-07  -2.240131E-03
  a    17  -2.254084E-07  -5.162576E-04  -1.040873E-06   2.957290E-07   2.782250E-07   1.117266E-07   1.599288E-03   4.360837E-07
  a    18  -1.147065E-04   1.074772E-08   4.192162E-08  -2.290238E-07  -9.493213E-08   1.940352E-03  -7.205009E-08  -2.997855E-04
  a    19  -1.297095E-06  -2.514662E-08  -1.260197E-03   7.133924E-09   1.694190E-10  -6.727801E-07   2.945184E-07   2.550495E-07
  a    20   1.316700E-03  -4.284207E-09  -1.318050E-06  -1.380004E-07  -1.829988E-07   8.509538E-04   5.426192E-09  -1.967791E-04
  a    21  -3.738351E-08  -6.382913E-04  -2.511186E-07   2.509080E-07  -2.412779E-08  -7.596592E-09  -3.047296E-04   1.159011E-08
  a    22   7.357336E-04  -2.441621E-08  -7.642009E-08  -3.387113E-07   3.122906E-08   4.322085E-04   1.450196E-08   1.766564E-04
  a    23   2.019067E-08   5.627158E-04   7.767162E-08  -1.758497E-07   1.057780E-07  -4.431099E-09   7.017821E-04  -4.037689E-08
  a    24   2.838027E-05  -7.464951E-09  -6.324579E-08   1.837410E-08  -2.914750E-07  -2.371015E-04   1.055650E-08   8.567100E-05
  a    25   3.388431E-03  -6.404250E-08   1.359737E-07  -1.384133E-07   9.853563E-09  -3.069543E-04  -3.621417E-08  -2.035182E-04
  a    26  -6.404250E-08   1.775497E-03  -5.027207E-07   1.177664E-06   3.099011E-07   1.986987E-08  -5.752898E-04   1.620660E-08
  a    27   1.359737E-07  -5.027207E-07   6.220052E-03  -1.730585E-08   6.439256E-09  -8.138761E-09  -2.079321E-07   4.373981E-11
  a    28  -1.384133E-07   1.177664E-06  -1.730585E-08   4.149502E-03   7.128952E-04   1.415738E-07   1.130621E-07   4.925515E-08
  a    29   9.853563E-09   3.099011E-07   6.439256E-09   7.128952E-04   2.748813E-03   3.271257E-08  -3.752826E-08   5.425079E-08
  a    30  -3.069543E-04   1.986987E-08  -8.138761E-09   1.415738E-07   3.271257E-08   2.256796E-03   1.691923E-09  -1.520310E-04
  a    31  -3.621417E-08  -5.752898E-04  -2.079321E-07   1.130621E-07  -3.752826E-08   1.691923E-09   1.880667E-03   1.730229E-08
  a    32  -2.035182E-04   1.620660E-08   4.373981E-11   4.925515E-08   5.425079E-08  -1.520310E-04   1.730229E-08   2.169465E-03
  a    33  -4.400819E-04   5.174874E-08   1.023264E-07  -7.932803E-09  -4.479697E-09   3.502133E-04  -6.277626E-10  -4.409393E-05
  a    34   1.293715E-07  -2.270123E-08   5.434603E-04  -5.981428E-09   1.213049E-09  -7.387390E-08  -3.125920E-08  -1.423775E-09
  a    35   6.455380E-04  -9.416834E-08  -1.428396E-08   1.311044E-08   1.796776E-08  -6.035514E-05   3.719066E-08   9.861326E-05
  a    36  -1.199720E-07  -5.133023E-04  -2.465874E-08   1.927133E-07   5.286402E-08   1.321574E-08   2.250447E-04  -2.384930E-08

               a    33        a    34        a    35        a    36
  a     3   7.400501E-04  -5.924247E-08   3.486177E-05   7.064396E-09
  a     4   4.504633E-03  -9.129967E-07   1.987005E-03  -1.438909E-07
  a     5  -2.125554E-07  -1.765362E-07   2.835815E-07   2.497168E-03
  a     6  -3.758154E-03   5.553232E-07   4.984705E-04  -1.177273E-07
  a     7  -1.816606E-07  -8.516277E-07   6.347045E-07   4.459501E-03
  a     8   1.650936E-06   5.506909E-03  -1.446255E-07   1.079015E-06
  a     9  -6.193541E-07  -4.899101E-03  -2.956392E-07  -2.245978E-06
  a    10   2.257324E-04  -6.759656E-08  -2.808890E-04   2.579793E-08
  a    11  -4.051278E-08   3.377459E-07  -1.571000E-08  -2.295422E-04
  a    12   8.273889E-04  -1.686753E-07  -5.479109E-05  -1.062622E-08
  a    13  -1.234160E-08   8.333752E-08   7.050199E-08   5.660066E-04
  a    14   6.440214E-05  -8.770237E-09   3.449333E-04  -5.028012E-08
  a    15   1.773594E-07  -1.698625E-04  -8.201378E-08  -1.654334E-10
  a    16   5.488014E-04   1.970780E-08  -1.487947E-04   9.500643E-08
  a    17  -1.572769E-07  -1.615504E-07   6.014801E-08   3.080533E-04
  a    18  -9.211895E-04   1.271746E-07  -3.688835E-04   3.673748E-08
  a    19  -1.494631E-06   6.719828E-04  -7.673650E-07   3.069218E-09
  a    20   1.807438E-03   2.200286E-07   7.188188E-04  -8.348458E-08
  a    21  -1.195548E-08   2.936365E-08  -5.708482E-08  -3.948962E-04
  a    22  -2.411576E-04   5.640847E-09  -4.153080E-04   6.287394E-08
  a    23   1.740867E-09  -5.851384E-09   4.448254E-08   1.124515E-04
  a    24   3.145736E-04  -9.843643E-08  -8.938465E-04   1.211818E-07
  a    25  -4.400819E-04   1.293715E-07   6.455380E-04  -1.199720E-07
  a    26   5.174874E-08  -2.270123E-08  -9.416834E-08  -5.133023E-04
  a    27   1.023264E-07   5.434603E-04  -1.428396E-08  -2.465874E-08
  a    28  -7.932803E-09  -5.981428E-09   1.311044E-08   1.927133E-07
  a    29  -4.479697E-09   1.213049E-09   1.796776E-08   5.286402E-08
  a    30   3.502133E-04  -7.387390E-08  -6.035514E-05   1.321574E-08
  a    31  -6.277626E-10  -3.125920E-08   3.719066E-08   2.250447E-04
  a    32  -4.409393E-05  -1.423775E-09   9.861326E-05  -2.384930E-08
  a    33   1.528620E-03   1.080547E-08   9.585041E-05  -1.109637E-08
  a    34   1.080547E-08   1.584340E-03  -2.387246E-08   1.714171E-08
  a    35   9.585041E-05  -2.387246E-08   1.309647E-03  -3.690123E-08
  a    36  -1.109637E-08   1.714171E-08  -3.690123E-08   1.010635E-03

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98372511     1.97227640     1.96534534     1.96288415     1.96076338     1.00665218
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     1.00446058     0.02025597     0.01787457     0.01652660     0.01574750     0.01057991     0.01048795     0.00741252
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00602127     0.00577586     0.00568976     0.00501212     0.00453788     0.00444851     0.00244981     0.00199856
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00184889     0.00143285     0.00137278     0.00104998     0.00088004     0.00065939     0.00053539     0.00033660
              MO    33       MO    34       MO    35       MO    36
  occ(*)=     0.00032534     0.00025556     0.00019376     0.00018351


 total number of electrons =   16.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
    N1_ s       1.998472  -0.000030   1.441584   0.083337   0.037896   0.000000
    N1_ p       0.000000  -0.000292  -0.002592   0.720928   0.434091   1.210261
    N1_ d       0.000000  -0.000105   0.004476   0.005722   0.001586   0.002354
    C1_ s       0.000539   1.999452   0.271624   0.782362   0.204885   0.000000
    C1_ p       0.000476   0.000003   0.040906   0.012172   0.793626   0.211922
    C1_ d      -0.000081   0.000001   0.006622   0.002986   0.012795   0.002846
    H1_ s       0.000301   0.000003   0.096105   0.101925   0.027047   0.231511
    H2_ s       0.000301   0.000003   0.096104   0.101919   0.027043   0.231523
    H3_ s      -0.000004   0.000483   0.014448   0.080463   0.213191   0.036233
    H4_ s      -0.000004   0.000483   0.014448   0.080463   0.213184   0.036235
 
   ao class       7a         8a         9a        10a        11a        12a  
    N1_ s       0.000000   0.000000   0.000000   0.001645   0.003239   0.000000
    N1_ p       1.264907   0.068437   0.339874   0.002302   0.000437   0.006641
    N1_ d       0.004354   0.006429   0.006056   0.000399   0.000874   0.001465
    C1_ s       0.000000   0.000000   0.000000   0.006927   0.000018   0.000000
    C1_ p       0.682294   0.450531   0.652548   0.000524   0.006868   0.000133
    C1_ d       0.009207   0.009819   0.005822   0.001967   0.001158   0.000684
    H1_ s       0.000001   0.034387   0.000012   0.000093   0.000966   0.003362
    H2_ s       0.000001   0.034387   0.000012   0.000093   0.000965   0.003362
    H3_ s       0.000000   0.201329   0.000069   0.003153   0.001675   0.000440
    H4_ s       0.000000   0.201333   0.000069   0.003153   0.001675   0.000440
 
   ao class      13a        14a        15a        16a        17a        18a  
    N1_ s       0.001843   0.000000   0.000000   0.000000   0.000516   0.000000
    N1_ p       0.005473   0.010130   0.002727   0.000083   0.001143   0.000114
    N1_ d       0.000713   0.000005   0.001110   0.000004   0.002330   0.003683
    C1_ s       0.000022   0.000000   0.000000   0.000000   0.000343   0.000000
    C1_ p       0.001359   0.000153   0.002782   0.007229   0.000757   0.000013
    C1_ d       0.000266   0.000292   0.000950   0.000097   0.000939   0.001966
    H1_ s       0.002768   0.000000   0.000065   0.000000  -0.000001  -0.000000
    H2_ s       0.002768   0.000000   0.000065   0.000000  -0.000001  -0.000000
    H3_ s       0.000268   0.000000   0.001394   0.000000  -0.000003   0.000000
    H4_ s       0.000268   0.000000   0.001394   0.000000  -0.000003   0.000000
 
   ao class      19a        20a        21a        22a        23a        24a  
    N1_ s       0.000225   0.000000   0.000742   0.000000   0.000000   0.000244
    N1_ p       0.000218   0.000635   0.000418   0.000000   0.000000   0.000259
    N1_ d       0.001178   0.000761   0.001924   0.003792   0.000361   0.000237
    C1_ s       0.002197   0.000000   0.000601   0.000000   0.000000  -0.000013
    C1_ p       0.000685   0.002541   0.000570   0.000000   0.000000   0.000205
    C1_ d       0.001125   0.000539   0.000293   0.000656   0.002088   0.001051
    H1_ s       0.000051   0.000018  -0.000016  -0.000000  -0.000000   0.000006
    H2_ s       0.000051   0.000018  -0.000016  -0.000000   0.000000   0.000006
    H3_ s      -0.000021   0.000250   0.000011   0.000000   0.000000   0.000003
    H4_ s      -0.000021   0.000250   0.000011  -0.000000  -0.000000   0.000003
 
   ao class      25a        26a        27a        28a        29a        30a  
    N1_ s      -0.000000   0.000288   0.000000   0.000000  -0.000002   0.000070
    N1_ p       0.000062   0.000016   0.000039   0.000198   0.000266  -0.000006
    N1_ d       0.000712   0.000435   0.000485   0.000141   0.000240   0.000062
    C1_ s      -0.000000   0.000108   0.000000  -0.000000   0.000035   0.000032
    C1_ p       0.000286   0.000108   0.000014   0.000098   0.000057  -0.000013
    C1_ d       0.000403   0.000131   0.000835   0.000279   0.000052   0.000252
    H1_ s       0.000046   0.000129   0.000000   0.000143   0.000099   0.000019
    H2_ s       0.000046   0.000129   0.000000   0.000143   0.000099   0.000019
    H3_ s       0.000147   0.000045   0.000000   0.000025   0.000017   0.000112
    H4_ s       0.000147   0.000045   0.000000   0.000025   0.000017   0.000112
 
   ao class      31a        32a        33a        34a        35a        36a  
    N1_ s      -0.000000   0.000000   0.000016   0.000000   0.000012   0.000001
    N1_ p       0.000019   0.000035   0.000025   0.000042   0.000006   0.000012
    N1_ d       0.000064  -0.000000   0.000006   0.000002   0.000004   0.000007
    C1_ s       0.000000   0.000000   0.000016   0.000000   0.000038   0.000007
    C1_ p      -0.000001   0.000005   0.000005   0.000103   0.000007   0.000078
    C1_ d       0.000128   0.000006   0.000015   0.000002   0.000016   0.000038
    H1_ s       0.000034   0.000126   0.000107   0.000009   0.000012   0.000004
    H2_ s       0.000034   0.000126   0.000106   0.000009   0.000012   0.000004
    H3_ s       0.000129   0.000019   0.000015   0.000044   0.000043   0.000016
    H4_ s       0.000129   0.000019   0.000015   0.000044   0.000043   0.000016


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.570097   3.269195   0.499327   0.499329   0.553994   0.553994
      p         4.066907   2.869043   0.000000   0.000000   0.000000   0.000000
      d         0.051867   0.066246   0.000000   0.000000   0.000000   0.000000
    total       7.688872   6.204483   0.499327   0.499329   0.553994   0.553994
 

 Total number of electrons:   16.00000000

 item #                     3 suffix=:.drt1.state3:
 read_civout: repnuc=  -88.6565513780480     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 96% root-following 0
 MR-CISD energy:   -94.66619491    -6.00964353
 residuum:     0.00037264
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.34532008    -5.68876870
 residuum:     0.00065006
 deltae:     0.00000090
 apxde:     0.00000021
================================================================================
  Reading record                      3  of civout
 INFO:ref#  3vector#  3 method:  0 last record  1max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.30933679    -5.65278541
 residuum:     0.00067046
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821    -0.01624692
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894     0.00491632
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907     0.10910721

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873513     0.00001296     0.01234876     0.02259949     0.02297308     0.08540722    -0.24621821     0.00000000
 ref:   2    -0.00000127     0.94789898    -0.00001002     0.00034756    -0.09539000     0.05418419     0.01079894     0.00000000
 ref:   3     0.00692423    -0.00002592     0.94559557    -0.03422374    -0.04416436    -0.11950148    -0.13618907     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 computing MRCISD density
 densi: densityinfo%a4den=   1.00000000000000     
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    3
--------------------------------------------------------------------------------
================================================================================
   DYZ=     215  DYX=     335  DYW=     401
   D0Z=     147  D0Y=     776  D0X=     489  D0W=     695
  DDZI=     109 DDYI=     378 DDXI=     272 DDWI=     359
  DDZE=       0 DDYE=      92 DDXE=      71 DDWE=      98
================================================================================
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

          modens reordered block   1

               a     1        a     2        a     3        a     4        a     5        a     6        a     7        a     8
  a     1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     3    0.00000        0.00000        1.98000       1.972217E-03   3.412576E-08   5.252042E-03   4.995578E-08  -1.761125E-06
  a     4    0.00000        0.00000       1.972217E-03    1.97264       2.133288E-08  -3.001737E-03   4.302987E-08  -1.895094E-06
  a     5    0.00000        0.00000       3.412576E-08   2.133288E-08    1.97124      -2.351221E-08  -1.665566E-03   3.680608E-05
  a     6    0.00000        0.00000       5.252042E-03  -3.001737E-03  -2.351221E-08    1.96161      -3.823508E-08   4.466699E-06
  a     7    0.00000        0.00000       4.995578E-08   4.302987E-08  -1.665566E-03  -3.823508E-08    1.96986       6.204498E-05
  a     8    0.00000        0.00000      -1.761125E-06  -1.895094E-06   3.680608E-05   4.466699E-06   6.204498E-05    1.10389    
  a     9    0.00000        0.00000       9.106038E-07   1.489441E-06   4.119581E-05  -1.297894E-06  -3.265360E-06  -0.478332    
  a    10    0.00000        0.00000      -5.668207E-03   6.808700E-03  -7.931668E-07  -6.091835E-03   3.818987E-07   1.653174E-06
  a    11    0.00000        0.00000       2.165396E-07  -4.426109E-07  -1.738683E-02   5.268224E-07   1.720472E-02   4.321553E-06
  a    12    0.00000        0.00000       1.547103E-03   4.533814E-03  -7.820884E-07  -4.709294E-03   4.436680E-07   6.205163E-07
  a    13    0.00000        0.00000       2.054115E-08  -6.581790E-08   2.821152E-03  -2.249276E-08  -5.804638E-03   2.271727E-06
  a    14    0.00000        0.00000       2.879287E-03   1.422667E-02  -1.978590E-07  -1.259842E-02   2.677998E-07  -7.002143E-08
  a    15    0.00000        0.00000      -7.275483E-07   3.640076E-06  -1.180527E-05   1.208096E-06  -1.883958E-06  -4.359679E-05
  a    16    0.00000        0.00000      -2.387440E-03   1.157799E-02   1.190578E-06   4.976425E-03  -3.063617E-06  -6.268749E-08
  a    17    0.00000        0.00000       3.931254E-07  -2.529384E-06   4.934203E-03  -9.757758E-07  -1.387625E-02  -2.265898E-06
  a    18    0.00000        0.00000      -4.750448E-03  -2.050408E-03  -3.580759E-08   9.425883E-03   3.741236E-07   7.385542E-07
  a    19    0.00000        0.00000      -1.180310E-07  -1.120812E-05   7.513104E-06   4.634830E-06  -1.574058E-06  -1.908038E-02
  a    20    0.00000        0.00000      -2.611869E-04   1.114631E-02  -1.700599E-07  -4.983781E-03   2.812417E-07  -1.693269E-05
  a    21    0.00000        0.00000       5.202335E-08  -1.602499E-07  -1.401155E-02   3.592557E-07   1.753258E-02   8.846972E-07
  a    22    0.00000        0.00000       1.377665E-03  -4.553503E-03  -2.056044E-08   8.901546E-03  -4.656023E-07   3.980691E-07
  a    23    0.00000        0.00000       2.833332E-08   2.031034E-07  -1.004692E-02  -1.952900E-07  -1.247087E-02   1.426884E-06
  a    24    0.00000        0.00000      -4.338487E-05   3.037478E-03  -2.039799E-07   2.136853E-03  -2.540609E-07   4.692056E-07
  a    25    0.00000        0.00000       3.517733E-03   3.670230E-03   1.408420E-07  -2.765748E-03  -2.149938E-07  -4.759774E-07
  a    26    0.00000        0.00000      -1.838892E-07  -1.686744E-08   2.814525E-03   1.222963E-07  -3.069693E-03   1.983984E-06
  a    27    0.00000        0.00000      -3.604370E-07  -7.750851E-08   3.114608E-06   8.555623E-07   3.613133E-06  -4.336028E-03
  a    28    0.00000        0.00000      -5.421458E-07  -2.013591E-06  -5.705126E-07   5.394570E-07   6.176155E-07   5.413795E-08
  a    29    0.00000        0.00000      -5.917533E-07   2.572977E-07   7.521952E-07   7.557922E-07  -1.035142E-07   1.969489E-09
  a    30    0.00000        0.00000       1.250204E-03  -8.316019E-04   1.525036E-08   1.986079E-03  -4.764925E-08   2.079863E-07
  a    31    0.00000        0.00000      -1.016863E-07   2.217059E-07   4.456003E-03   7.066502E-10  -2.890695E-03   4.097943E-07
  a    32    0.00000        0.00000      -1.973636E-03   2.751084E-03  -2.521752E-07   7.103492E-04   1.315798E-07  -4.556600E-08
  a    33    0.00000        0.00000       1.494906E-03   4.475755E-03  -1.496727E-07  -5.802263E-04  -4.498201E-08  -1.090428E-06
  a    34    0.00000        0.00000      -4.128094E-07  -1.121585E-06   7.698385E-07  -2.113447E-07   2.974242E-06  -6.407300E-03
  a    35    0.00000        0.00000      -7.024663E-04  -1.139219E-03   4.679920E-07   1.298973E-03   7.094772E-07   5.780489E-08
  a    36    0.00000        0.00000       9.163975E-08   2.456112E-07   3.588183E-03  -1.937306E-07   5.224028E-03   3.149262E-07

               a     9        a    10        a    11        a    12        a    13        a    14        a    15        a    16
  a     3   9.106038E-07  -5.668207E-03   2.165396E-07   1.547103E-03   2.054115E-08   2.879287E-03  -7.275483E-07  -2.387440E-03
  a     4   1.489441E-06   6.808700E-03  -4.426109E-07   4.533814E-03  -6.581790E-08   1.422667E-02   3.640076E-06   1.157799E-02
  a     5   4.119581E-05  -7.931668E-07  -1.738683E-02  -7.820884E-07   2.821152E-03  -1.978590E-07  -1.180527E-05   1.190578E-06
  a     6  -1.297894E-06  -6.091835E-03   5.268224E-07  -4.709294E-03  -2.249276E-08  -1.259842E-02   1.208096E-06   4.976425E-03
  a     7  -3.265360E-06   3.818987E-07   1.720472E-02   4.436680E-07  -5.804638E-03   2.677998E-07  -1.883958E-06  -3.063617E-06
  a     8  -0.478332       1.653174E-06   4.321553E-06   6.205163E-07   2.271727E-06  -7.002143E-08  -4.359679E-05  -6.268749E-08
  a     9   0.881208       2.322967E-06   9.280176E-06  -6.418398E-08  -8.569338E-10   5.093934E-07  -4.683681E-03   1.858056E-06
  a    10   2.322967E-06   1.121025E-02  -3.434395E-09  -4.579325E-04  -6.489052E-09   7.071060E-04   1.024515E-06   3.137180E-03
  a    11   9.280176E-06  -3.434395E-09   1.025126E-02   1.405166E-08   1.316964E-04  -8.335015E-08   5.757815E-08  -6.894232E-07
  a    12  -6.418398E-08  -4.579325E-04   1.405166E-08   9.637218E-03  -1.990922E-08   1.041610E-03   3.196816E-07   1.068613E-03
  a    13  -8.569338E-10  -6.489052E-09   1.316964E-04  -1.990922E-08   9.131745E-03  -1.476070E-07  -6.736457E-08  -3.187765E-07
  a    14   5.093934E-07   7.071060E-04  -8.335015E-08   1.041610E-03  -1.476070E-07   1.172712E-02   2.054613E-06   6.970556E-03
  a    15  -4.683681E-03   1.024515E-06   5.757815E-08   3.196816E-07  -6.736457E-08   2.054613E-06   3.723609E-03   2.131105E-06
  a    16   1.858056E-06   3.137180E-03  -6.894232E-07   1.068613E-03  -3.187765E-07   6.970556E-03   2.131105E-06   1.091680E-02
  a    17   1.017614E-06  -6.538798E-07  -2.125649E-03  -3.193614E-07  -8.978743E-04  -1.610998E-06  -2.436859E-07  -4.792054E-07
  a    18   4.865679E-07   7.573000E-04   6.470979E-08  -6.036690E-04   8.150868E-08  -1.136131E-03  -7.772130E-08  -2.336751E-04
  a    19  -1.516768E-02  -3.747503E-06  -1.172709E-07  -3.668285E-06  -2.189924E-08  -5.986721E-06   1.099258E-03  -5.702919E-06
  a    20  -1.398882E-05   3.576400E-03  -2.215734E-07   3.737532E-03  -1.250452E-07   6.222761E-03   2.705101E-06   5.568161E-03
  a    21   4.029726E-06   3.064729E-07   6.896805E-03   1.971682E-07  -3.524351E-03  -4.411560E-08  -7.879062E-08  -3.239667E-07
  a    22   1.644514E-07   5.864828E-03   9.448023E-09  -4.263828E-03   2.470912E-07  -3.006999E-03  -1.865266E-07  -7.499414E-04
  a    23   1.120660E-06  -8.137441E-09   4.224853E-03   1.767423E-07   7.236386E-03   2.018260E-07   2.329766E-07  -6.853114E-08
  a    24  -6.711805E-07   2.218242E-03  -1.564706E-07   4.233283E-03   2.052704E-07  -6.650210E-03  -9.843310E-07  -3.569793E-03
  a    25   5.154443E-07   3.486773E-03  -2.757529E-07   2.084557E-03   1.381629E-08   2.361352E-03   5.096430E-07   1.622809E-03
  a    26  -2.176352E-06  -1.831536E-07  -1.096790E-03  -1.085746E-07   1.111210E-04  -6.739347E-08   3.289348E-07  -3.151697E-07
  a    27   8.776091E-04  -7.002724E-08  -2.742444E-08  -7.207144E-08  -4.272232E-07  -3.283401E-08  -1.757010E-04   4.323867E-08
  a    28   1.370559E-08   4.659803E-08   4.801475E-07   9.553767E-08   1.321562E-07  -1.230282E-07  -3.846954E-09  -7.453501E-09
  a    29  -3.665791E-08  -1.370969E-07   1.942341E-07  -1.112953E-07   9.467636E-08  -1.497377E-07  -2.203856E-09  -5.916476E-08
  a    30  -6.494830E-07   4.787308E-04  -2.995959E-08   7.061278E-04   4.902897E-10   5.964284E-04   4.980503E-08   1.995246E-04
  a    31   2.335129E-06   1.739075E-08   7.007892E-04   1.330218E-08   8.748869E-05  -2.434477E-08  -1.427535E-07   5.020315E-07
  a    32  -6.212446E-08  -6.249977E-04   7.881371E-09  -1.250234E-04   1.182273E-09  -3.079545E-04  -8.128349E-07  -2.524813E-03
  a    33  -6.186399E-07   5.813166E-04  -1.369652E-08  -2.451794E-04  -2.691880E-08   1.403440E-04   2.261521E-07   6.748578E-04
  a    34   8.301112E-05  -1.417162E-07   7.204999E-08   4.655993E-08  -1.454782E-08   3.961133E-08   7.793391E-05  -3.901377E-08
  a    35   4.346668E-07  -2.309512E-04  -1.398141E-08  -2.608358E-05   7.491262E-08   7.262157E-04   1.991654E-08   1.841504E-04
  a    36   1.010529E-07   2.346058E-08  -1.861352E-04  -1.989411E-08   6.385018E-04  -9.834577E-08  -1.686901E-08   3.756949E-08

               a    17        a    18        a    19        a    20        a    21        a    22        a    23        a    24
  a     3   3.931254E-07  -4.750448E-03  -1.180310E-07  -2.611869E-04   5.202335E-08   1.377665E-03   2.833332E-08  -4.338487E-05
  a     4  -2.529384E-06  -2.050408E-03  -1.120812E-05   1.114631E-02  -1.602499E-07  -4.553503E-03   2.031034E-07   3.037478E-03
  a     5   4.934203E-03  -3.580759E-08   7.513104E-06  -1.700599E-07  -1.401155E-02  -2.056044E-08  -1.004692E-02  -2.039799E-07
  a     6  -9.757758E-07   9.425883E-03   4.634830E-06  -4.983781E-03   3.592557E-07   8.901546E-03  -1.952900E-07   2.136853E-03
  a     7  -1.387625E-02   3.741236E-07  -1.574058E-06   2.812417E-07   1.753258E-02  -4.656023E-07  -1.247087E-02  -2.540609E-07
  a     8  -2.265898E-06   7.385542E-07  -1.908038E-02  -1.693269E-05   8.846972E-07   3.980691E-07   1.426884E-06   4.692056E-07
  a     9   1.017614E-06   4.865679E-07  -1.516768E-02  -1.398882E-05   4.029726E-06   1.644514E-07   1.120660E-06  -6.711805E-07
  a    10  -6.538798E-07   7.573000E-04  -3.747503E-06   3.576400E-03   3.064729E-07   5.864828E-03  -8.137441E-09   2.218242E-03
  a    11  -2.125649E-03   6.470979E-08  -1.172709E-07  -2.215734E-07   6.896805E-03   9.448023E-09   4.224853E-03  -1.564706E-07
  a    12  -3.193614E-07  -6.036690E-04  -3.668285E-06   3.737532E-03   1.971682E-07  -4.263828E-03   1.767423E-07   4.233283E-03
  a    13  -8.978743E-04   8.150868E-08  -2.189924E-08  -1.250452E-07  -3.524351E-03   2.470912E-07   7.236386E-03   2.052704E-07
  a    14  -1.610998E-06  -1.136131E-03  -5.986721E-06   6.222761E-03  -4.411560E-08  -3.006999E-03   2.018260E-07  -6.650210E-03
  a    15  -2.436859E-07  -7.772130E-08   1.099258E-03   2.705101E-06  -7.879062E-08  -1.865266E-07   2.329766E-07  -9.843310E-07
  a    16  -4.792054E-07  -2.336751E-04  -5.702919E-06   5.568161E-03  -3.239667E-07  -7.499414E-04  -6.853114E-08  -3.569793E-03
  a    17   8.521886E-03  -4.540058E-08   8.651947E-07  -1.260386E-06  -1.344470E-03   2.257689E-07  -7.443184E-04   7.979010E-07
  a    18  -4.540058E-08   4.301995E-03   2.404777E-07  -4.133790E-04   6.009580E-08   1.678932E-03   5.177975E-09   4.471763E-04
  a    19   8.651947E-07   2.404777E-07   9.647450E-03   1.148059E-07   2.688668E-08   9.415677E-07   2.065027E-07   1.288648E-06
  a    20  -1.260386E-06  -4.133790E-04   1.148059E-07   9.316444E-03   7.170103E-09  -1.005661E-03   7.273723E-08  -1.341658E-03
  a    21  -1.344470E-03   6.009580E-08   2.688668E-08   7.170103E-09   6.667923E-03   1.500343E-08  -3.582938E-04   1.264690E-08
  a    22   2.257689E-07   1.678932E-03   9.415677E-07  -1.005661E-03   1.500343E-08   6.755739E-03   4.722927E-08   1.385890E-03
  a    23  -7.443184E-04   5.177975E-09   2.065027E-07   7.273723E-08  -3.582938E-04   4.722927E-08   9.319607E-03  -5.126042E-08
  a    24   7.979010E-07   4.471763E-04   1.288648E-06  -1.341658E-03   1.264690E-08   1.385890E-03  -5.126042E-08   8.200262E-03
  a    25  -4.091490E-07  -2.174614E-04  -1.695031E-06   1.799894E-03  -6.488338E-08   4.721746E-04   2.534711E-08   6.247084E-05
  a    26  -1.350234E-03   4.667092E-08   5.529301E-08  -1.929656E-08  -1.030570E-03  -1.249052E-08   3.775429E-04  -2.989891E-08
  a    27  -2.705340E-07   6.481128E-08  -4.441501E-04  -5.082268E-07   2.153831E-07  -3.342065E-08  -3.996633E-07  -5.810436E-08
  a    28   7.204041E-07  -8.968613E-09  -2.473503E-09   1.138102E-07   4.289245E-07  -1.153806E-07  -7.903765E-08   2.203037E-07
  a    29   5.441531E-07  -4.620020E-08   8.447963E-10  -2.637467E-07   1.112432E-07   2.663202E-07   1.483805E-07  -3.732784E-08
  a    30   1.611694E-08   8.946002E-04  -1.286812E-06   1.459153E-03  -1.491811E-08  -4.639625E-04   3.172992E-08   2.329315E-05
  a    31   2.609848E-03  -6.815395E-08   7.727308E-08  -1.220014E-08   2.950134E-04   2.884726E-08   7.430169E-04   1.455952E-08
  a    32   4.594510E-07  -7.505950E-05   5.620609E-07  -5.349808E-04  -1.063955E-08   1.016333E-04  -5.080946E-08   3.490488E-04
  a    33  -1.681600E-07  -2.166351E-04  -1.209198E-06   1.540093E-03   1.791003E-09   3.972489E-04  -1.672475E-08  -1.905048E-04
  a    34   2.937893E-08  -7.396440E-09   8.277549E-04   4.597257E-07   4.542424E-09  -1.289817E-07  -4.176070E-09  -3.223121E-08
  a    35  -2.624912E-08  -4.041275E-04  -9.709315E-07   9.517392E-04  -6.198361E-08  -5.000676E-04   5.367354E-08  -1.138529E-03
  a    36   2.333544E-04   4.874993E-08  -8.005380E-09  -1.131315E-07  -4.191043E-04   7.913537E-08   1.152237E-04   1.522244E-07

               a    25        a    26        a    27        a    28        a    29        a    30        a    31        a    32
  a     3   3.517733E-03  -1.838892E-07  -3.604370E-07  -5.421458E-07  -5.917533E-07   1.250204E-03  -1.016863E-07  -1.973636E-03
  a     4   3.670230E-03  -1.686744E-08  -7.750851E-08  -2.013591E-06   2.572977E-07  -8.316019E-04   2.217059E-07   2.751084E-03
  a     5   1.408420E-07   2.814525E-03   3.114608E-06  -5.705126E-07   7.521952E-07   1.525036E-08   4.456003E-03  -2.521752E-07
  a     6  -2.765748E-03   1.222963E-07   8.555623E-07   5.394570E-07   7.557922E-07   1.986079E-03   7.066502E-10   7.103492E-04
  a     7  -2.149938E-07  -3.069693E-03   3.613133E-06   6.176155E-07  -1.035142E-07  -4.764925E-08  -2.890695E-03   1.315798E-07
  a     8  -4.759774E-07   1.983984E-06  -4.336028E-03   5.413795E-08   1.969489E-09   2.079863E-07   4.097943E-07  -4.556600E-08
  a     9   5.154443E-07  -2.176352E-06   8.776091E-04   1.370559E-08  -3.665791E-08  -6.494830E-07   2.335129E-06  -6.212446E-08
  a    10   3.486773E-03  -1.831536E-07  -7.002724E-08   4.659803E-08  -1.370969E-07   4.787308E-04   1.739075E-08  -6.249977E-04
  a    11  -2.757529E-07  -1.096790E-03  -2.742444E-08   4.801475E-07   1.942341E-07  -2.995959E-08   7.007892E-04   7.881371E-09
  a    12   2.084557E-03  -1.085746E-07  -7.207144E-08   9.553767E-08  -1.112953E-07   7.061278E-04   1.330218E-08  -1.250234E-04
  a    13   1.381629E-08   1.111210E-04  -4.272232E-07   1.321562E-07   9.467636E-08   4.902897E-10   8.748869E-05   1.182273E-09
  a    14   2.361352E-03  -6.739347E-08  -3.283401E-08  -1.230282E-07  -1.497377E-07   5.964284E-04  -2.434477E-08  -3.079545E-04
  a    15   5.096430E-07   3.289348E-07  -1.757010E-04  -3.846954E-09  -2.203856E-09   4.980503E-08  -1.427535E-07  -8.128349E-07
  a    16   1.622809E-03  -3.151697E-07   4.323867E-08  -7.453501E-09  -5.916476E-08   1.995246E-04   5.020315E-07  -2.524813E-03
  a    17  -4.091490E-07  -1.350234E-03  -2.705340E-07   7.204041E-07   5.441531E-07   1.611694E-08   2.609848E-03   4.594510E-07
  a    18  -2.174614E-04   4.667092E-08   6.481128E-08  -8.968613E-09  -4.620020E-08   8.946002E-04  -6.815395E-08  -7.505950E-05
  a    19  -1.695031E-06   5.529301E-08  -4.441501E-04  -2.473503E-09   8.447963E-10  -1.286812E-06   7.727308E-08   5.620609E-07
  a    20   1.799894E-03  -1.929656E-08  -5.082268E-07   1.138102E-07  -2.637467E-07   1.459153E-03  -1.220014E-08  -5.349808E-04
  a    21  -6.488338E-08  -1.030570E-03   2.153831E-07   4.289245E-07   1.112432E-07  -1.491811E-08   2.950134E-04  -1.063955E-08
  a    22   4.721746E-04  -1.249052E-08  -3.342065E-08  -1.153806E-07   2.663202E-07  -4.639625E-04   2.884726E-08   1.016333E-04
  a    23   2.534711E-08   3.775429E-04  -3.996633E-07  -7.903765E-08   1.483805E-07   3.172992E-08   7.430169E-04  -5.080946E-08
  a    24   6.247084E-05  -2.989891E-08  -5.810436E-08   2.203037E-07  -3.732784E-08   2.329315E-05   1.455952E-08   3.490488E-04
  a    25   3.752120E-03  -6.780818E-08   7.123046E-08  -6.559622E-08   4.857611E-08  -9.365955E-05  -6.189229E-08  -3.966508E-04
  a    26  -6.780818E-08   1.970901E-03  -6.516570E-07   8.831762E-07   2.360156E-07   9.923742E-09  -9.051982E-04   3.379131E-08
  a    27   7.123046E-08  -6.516570E-07   5.259271E-03  -1.444333E-08   5.726280E-09   3.175577E-08  -9.225431E-08   1.768833E-08
  a    28  -6.559622E-08   8.831762E-07  -1.444333E-08   3.763999E-03   6.716824E-04   9.934693E-08   2.930106E-07  -3.258033E-10
  a    29   4.857611E-08   2.360156E-07   5.726280E-09   6.716824E-04   3.235519E-03   3.603200E-08   3.530385E-08  -9.092877E-08
  a    30  -9.365955E-05   9.923742E-09   3.175577E-08   9.934693E-08   3.603200E-08   1.919714E-03   2.349035E-08   7.720147E-05
  a    31  -6.189229E-08  -9.051982E-04  -9.225431E-08   2.930106E-07   3.530385E-08   2.349035E-08   2.474563E-03  -1.234011E-08
  a    32  -3.966508E-04   3.379131E-08   1.768833E-08  -3.258033E-10  -9.092877E-08   7.720147E-05  -1.234011E-08   2.169921E-03
  a    33  -5.444015E-04   6.305686E-08   1.778299E-07  -1.962580E-08  -5.719133E-08   5.003896E-04  -9.085651E-09  -2.217299E-04
  a    34   1.605885E-07  -1.337343E-07   7.184697E-04  -9.437677E-09   3.920440E-09  -1.096891E-07  -4.457300E-08   2.144556E-08
  a    35   7.521162E-04  -1.142075E-07  -2.497449E-08  -1.647167E-08  -6.628432E-09  -1.058350E-04   2.682941E-08   4.288333E-05
  a    36  -1.421514E-07  -6.591901E-04  -9.229827E-08   2.675786E-07   4.782120E-08   1.893642E-08   1.692817E-04  -1.553036E-08

               a    33        a    34        a    35        a    36
  a     3   1.494906E-03  -4.128094E-07  -7.024663E-04   9.163975E-08
  a     4   4.475755E-03  -1.121585E-06  -1.139219E-03   2.456112E-07
  a     5  -1.496727E-07   7.698385E-07   4.679920E-07   3.588183E-03
  a     6  -5.802263E-04  -2.113447E-07   1.298973E-03  -1.937306E-07
  a     7  -4.498201E-08   2.974242E-06   7.094772E-07   5.224028E-03
  a     8  -1.090428E-06  -6.407300E-03   5.780489E-08   3.149262E-07
  a     9  -6.186399E-07   8.301112E-05   4.346668E-07   1.010529E-07
  a    10   5.813166E-04  -1.417162E-07  -2.309512E-04   2.346058E-08
  a    11  -1.369652E-08   7.204999E-08  -1.398141E-08  -1.861352E-04
  a    12  -2.451794E-04   4.655993E-08  -2.608358E-05  -1.989411E-08
  a    13  -2.691880E-08  -1.454782E-08   7.491262E-08   6.385018E-04
  a    14   1.403440E-04   3.961133E-08   7.262157E-04  -9.834577E-08
  a    15   2.261521E-07   7.793391E-05   1.991654E-08  -1.686901E-08
  a    16   6.748578E-04  -3.901377E-08   1.841504E-04   3.756949E-08
  a    17  -1.681600E-07   2.937893E-08  -2.624912E-08   2.333544E-04
  a    18  -2.166351E-04  -7.396440E-09  -4.041275E-04   4.874993E-08
  a    19  -1.209198E-06   8.277549E-04  -9.709315E-07  -8.005380E-09
  a    20   1.540093E-03   4.597257E-07   9.517392E-04  -1.131315E-07
  a    21   1.791003E-09   4.542424E-09  -6.198361E-08  -4.191043E-04
  a    22   3.972489E-04  -1.289817E-07  -5.000676E-04   7.913537E-08
  a    23  -1.672475E-08  -4.176070E-09   5.367354E-08   1.152237E-04
  a    24  -1.905048E-04  -3.223121E-08  -1.138529E-03   1.522244E-07
  a    25  -5.444015E-04   1.605885E-07   7.521162E-04  -1.421514E-07
  a    26   6.305686E-08  -1.337343E-07  -1.142075E-07  -6.591901E-04
  a    27   1.778299E-07   7.184697E-04  -2.497449E-08  -9.229827E-08
  a    28  -1.962580E-08  -9.437677E-09  -1.647167E-08   2.675786E-07
  a    29  -5.719133E-08   3.920440E-09  -6.628432E-09   4.782120E-08
  a    30   5.003896E-04  -1.096891E-07  -1.058350E-04   1.893642E-08
  a    31  -9.085651E-09  -4.457300E-08   2.682941E-08   1.692817E-04
  a    32  -2.217299E-04   2.144556E-08   4.288333E-05  -1.553036E-08
  a    33   1.640298E-03  -2.638960E-08   1.772314E-04  -2.068229E-08
  a    34  -2.638960E-08   1.469441E-03  -3.804400E-08   7.891041E-09
  a    35   1.772314E-04  -3.804400E-08   1.393442E-03  -2.370573E-08
  a    36  -2.068229E-08   7.891041E-09  -2.370573E-08   1.181792E-03

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98158277     1.97375001     1.97299837     1.96896568     1.95950256     1.48372426
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.50264121     0.02783423     0.01824837     0.01804959     0.01592478     0.01532832     0.00941142     0.00878780
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00553677     0.00530382     0.00492300     0.00422155     0.00409748     0.00351204     0.00277797     0.00231942
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00216144     0.00157157     0.00124386     0.00123121     0.00096384     0.00082033     0.00081141     0.00047587
              MO    33       MO    34       MO    35       MO    36
  occ(*)=     0.00041046     0.00039597     0.00027996     0.00019264


 total number of electrons =   16.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
    N1_ s       1.998472  -0.000030   1.425980   0.124798   0.000000   0.000000
    N1_ p       0.000000  -0.000292   0.002547   0.175740   0.500614   0.882170
    N1_ d       0.000000  -0.000105   0.004283   0.003111   0.000676   0.015316
    C1_ s       0.000539   1.999452   0.257071   1.018458   0.000000   0.000000
    C1_ p       0.000476   0.000003  -0.006098   0.123485   0.919043   0.319784
    C1_ d      -0.000081   0.000001   0.005608   0.005723   0.000248   0.017287
    H1_ s       0.000301   0.000003   0.118875   0.042105   0.058680   0.228762
    H2_ s       0.000301   0.000003   0.118874   0.042098   0.058677   0.228774
    H3_ s      -0.000004   0.000483   0.027220   0.219126   0.217529   0.138429
    H4_ s      -0.000004   0.000483   0.027222   0.219106   0.217533   0.138444
 
   ao class       7a         8a         9a        10a        11a        12a  
    N1_ s       0.036000   0.000000  -0.000000   0.002856   0.000000   0.003908
    N1_ p       0.962977   1.417911   0.007026   0.007068   0.010116   0.003270
    N1_ d       0.006361  -0.000740   0.005105   0.001516   0.000432   0.000616
    C1_ s       0.015844   0.000000   0.000000   0.004065   0.000000   0.001079
    C1_ p       0.686693   0.050940   0.490469   0.011762   0.000305   0.000500
    C1_ d       0.015769   0.015614   0.000042   0.000471   0.000163   0.000055
    H1_ s       0.058306   0.000000  -0.000000  -0.000011   0.003599   0.003423
    H2_ s       0.058307   0.000000  -0.000000  -0.000011   0.003601   0.003421
    H3_ s       0.059622   0.000000  -0.000000   0.000059   0.000016   0.000888
    H4_ s       0.059623   0.000000   0.000000   0.000059   0.000016   0.000888
 
   ao class      13a        14a        15a        16a        17a        18a  
    N1_ s       0.000000   0.001016   0.000000   0.000000   0.000328   0.000000
    N1_ p       0.000068   0.000307   0.002020   0.009131   0.001351   0.000031
    N1_ d       0.001069   0.000195   0.002020   0.000004   0.002245   0.003638
    C1_ s       0.000000   0.003927   0.000000   0.000000  -0.000124   0.000000
    C1_ p       0.006880   0.001485   0.003543  -0.000482   0.000805   0.000164
    C1_ d       0.000612   0.000134   0.001721   0.000136   0.000887   0.001472
    H1_ s       0.000244   0.000776   0.000060   0.000000   0.000010  -0.000000
    H2_ s       0.000244   0.000776   0.000060   0.000000   0.000010   0.000000
    H3_ s       0.003404   0.003356  -0.000006  -0.000000   0.000013   0.000000
    H4_ s       0.003404   0.003356  -0.000006  -0.000000   0.000013   0.000000
 
   ao class      19a        20a        21a        22a        23a        24a  
    N1_ s       0.001019  -0.000000   0.000052   0.000000  -0.000000   0.000000
    N1_ p       0.000187   0.000000   0.000123  -0.000128   0.000000   0.000071
    N1_ d       0.002921   0.004071   0.000312   0.000047   0.000099   0.000862
    C1_ s       0.000235   0.000000   0.001285  -0.000000   0.000000   0.000000
    C1_ p      -0.000016   0.000000   0.000339   0.003552   0.000000   0.000343
    C1_ d       0.000187   0.000151   0.001331   0.000041   0.002679   0.000343
    H1_ s       0.000111   0.000000  -0.000007   0.000000   0.000000   0.000089
    H2_ s       0.000111  -0.000000  -0.000007   0.000000   0.000000   0.000089
    H3_ s       0.000085   0.000000   0.000335   0.000000   0.000000   0.000261
    H4_ s       0.000085   0.000000   0.000335   0.000000   0.000000   0.000261
 
   ao class      25a        26a        27a        28a        29a        30a  
    N1_ s       0.000069   0.000199   0.000000   0.000000   0.000050   0.000118
    N1_ p       0.000305   0.000019   0.000141   0.000027   0.000241   0.000084
    N1_ d       0.000268   0.000289   0.000111   0.000358   0.000230   0.000134
    C1_ s       0.000216  -0.000126   0.000000   0.000000   0.000025  -0.000008
    C1_ p       0.000366   0.000041   0.000145   0.000003   0.000166  -0.000017
    C1_ d       0.000883   0.000307   0.000490   0.000843   0.000103   0.000175
    H1_ s      -0.000014   0.000285   0.000134   0.000000   0.000038   0.000057
    H2_ s      -0.000014   0.000285   0.000134   0.000000   0.000038   0.000057
    H3_ s       0.000042   0.000136   0.000045   0.000000   0.000036   0.000111
    H4_ s       0.000042   0.000136   0.000045   0.000000   0.000036   0.000111
 
   ao class      31a        32a        33a        34a        35a        36a  
    N1_ s       0.000000   0.000076   0.000000   0.000004   0.000000   0.000002
    N1_ p      -0.000012   0.000001   0.000156   0.000046   0.000004   0.000005
    N1_ d       0.000159   0.000010   0.000012   0.000018   0.000000   0.000008
    C1_ s      -0.000000   0.000076   0.000000  -0.000001   0.000000   0.000068
    C1_ p      -0.000025   0.000106   0.000060   0.000030   0.000046   0.000033
    C1_ d       0.000161   0.000086   0.000000   0.000032   0.000006   0.000048
    H1_ s       0.000115   0.000049   0.000105   0.000095   0.000025  -0.000000
    H2_ s       0.000115   0.000049   0.000105   0.000095   0.000025  -0.000000
    H3_ s       0.000149   0.000011  -0.000014   0.000039   0.000086   0.000014
    H4_ s       0.000148   0.000011  -0.000014   0.000039   0.000086   0.000014


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.594917   3.302080   0.516215   0.516217   0.671470   0.671473
      p         3.983324   2.614930   0.000000   0.000000   0.000000   0.000000
      d         0.055648   0.073725   0.000000   0.000000   0.000000   0.000000
    total       7.633889   5.990735   0.516215   0.516217   0.671470   0.671473
 

 Total number of electrons:   16.00000000

