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
  RTOLBK = 1e-4,1e-4,1e-4,
  NITER = 90
  NVCIMN = 5
  RTOLCI = 1e-4,1e-4,1e-4,
  NVCIMX = 8
  NVRFMX = 8
  NVBKMX = 8
   iden=2
  CSFPRN = 10,
 /&end
 transition
   1  2  1  3
 ------------------------------------------------------------------------
transition densities: resetting nroot to    3
lodens (list->root)=  2  3
invlodens (root->list)= -1  1  2
 bummer (warning):resetting fileloc for seriel operation0
 USING SEGMENTS OF EQUAL SIZE

****************  list of control variables  ****************
 lvlprt =    0      nroot  =    3      noldv  =   0      noldhv =   0
 nunitv =    5      nbkitr =    1      niter  =  90      davcor =  10
 csfprn =   10      ivmode =    3      istrt  =   0      vout   =   0
 iortls =    0      nvbkmx =    8      ibktv  =  -1      ibkthv =  -1
 nvcimx =    8      icitv  =   -1      icithv =  -1      frcsub =   0
 nvbkmn =    3      nvcimn =    5      maxseg = 300      nrfitr =  30
 ncorel =   12      nvrfmx =    8      nvrfmn =   5      iden   =   2
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
    1        1.000E-04    1.000E-04
    2        1.000E-04    1.000E-04
    3        1.000E-04    1.000E-04
 Computing density:                    .drt1.state2
 Computing density:                    .drt1.state3
 Computing transition density:          drt1.state2-> drt1.state3 (.trd2to3)
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
 Hermit Integral Program : SIFS version  srv-p22-13.cbls.c 12:57:23.670 24-Jun-21
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =   38.588704248                                          
 MCSCF energy =     -94.155328507                                                
 SIFS file created by program tran.      srv-p22-13.cbls.c 12:57:26.353 24-Jun-21

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
  ... transition density matrix calculation requested 
  ... setting rmu(*) to zero 

 new core energy added to the energy(*) list.
 from the integral file: h1_core= -9.226324016102E+01

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
 pro1e: eref =    0.000000000000000E+00
 total size of srtscr:                     3  records of                 524288 
 WP =              12582912 Bytes
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
 !timer: cisrt complete                  cpu_time=     0.176 walltime=     0.724
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
 energy( 2)= -9.226324016102E+01, ietype=    6,   fcore energy of type: H1(*)   

 total core repulsion energy = -5.367453591285E+01
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
 bummer (warning):transition densities: resetting ref occupation number to 00
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
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:           0    task #     2:           0    task #     3:           0    task #     4:           0
task #     5:           0    task #     6:           0    task #     7:           0    task #     8:        1047
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         524    task #    18:           0    task #    19:           0    task #    20:           0
 reference space has dimension       6
 dsyevx: computed roots 1 to    6(converged:   6)

    root           eigenvalues
    ----           ------------
       1         -94.3985323836
       2         -94.0535597859
       3         -94.0138933504
       4         -93.7714030970
       5         -93.6147885713
       6         -93.5309041036

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
 residual norm convergence criteria:               0.000100  0.000100  0.000100

          starting bk iteration   1

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2797 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        2305    task #     6:           0    task #     7:           0    task #     8:        1047
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2797 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        2305    task #     6:           0    task #     7:           0    task #     8:        1047
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2797 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        2305    task #     6:           0    task #     7:           0    task #     8:        1047
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2797 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        2305    task #     6:           0    task #     7:           0    task #     8:        1047
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :        2797 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:         265 wz:         421 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:           0    task #     4:           0
task #     5:        2305    task #     6:           0    task #     7:           0    task #     8:        1047
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1   -40.72399647
   ht   2     0.00000000   -40.37902387
   ht   3     0.00000000    -0.00000000   -40.33935744
   ht   4    -0.00000000    -0.00000000    -0.00000000   -40.09686718
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -39.94025266

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1    1.00000      -7.042027E-14   3.600792E-14   2.991463E-14  -3.469447E-18
 ref    2  -7.050132E-14   -1.00000      -2.082870E-13   9.985566E-14   8.470329E-22
 ref    3  -3.588073E-14  -2.085924E-13    1.00000       2.596568E-16   4.336809E-18

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1    1.00000        1.00000        1.00000       1.086610E-26   3.084497E-35

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     1.00000000    -0.00000000     0.00000000     0.00000000    -0.00000000
 ref:   2    -0.00000000    -1.00000000    -0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000    -0.00000000     1.00000000     0.00000000     0.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3985323836  2.1316E-14  3.1591E-01  9.7425E-01  1.0000E-04   
 mr-sdci #  1  2    -94.0535597859  2.1316E-14  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  1  3    -94.0138933504 -4.9738E-14  0.0000E+00  9.9686E-01  1.0000E-04   
 mr-sdci #  1  4    -93.7714030970 -6.3949E-14  0.0000E+00  1.0037E+00  1.0000E-04   
 mr-sdci #  1  5    -93.6147885713 -2.1316E-14  0.0000E+00  1.0213E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000167
time for cinew                         0.006707
time for eigenvalue solver             0.000000
time for vector access                 0.000001

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3985323836  2.1316E-14  3.1591E-01  9.7425E-01  1.0000E-04   
 mr-sdci #  1  2    -94.0535597859  2.1316E-14  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  1  3    -94.0138933504 -4.9738E-14  0.0000E+00  9.9686E-01  1.0000E-04   
 mr-sdci #  1  4    -93.7714030970 -6.3949E-14  0.0000E+00  1.0037E+00  1.0000E-04   
 mr-sdci #  1  5    -93.6147885713 -2.1316E-14  0.0000E+00  1.0213E+00  1.0000E-04   
 
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
 number of iterations:                                   90
 residual norm convergence criteria:               0.000100  0.000100  0.000100

          starting ci iteration   1

 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1   -40.72399647
   ht   2    -0.00000000   -40.37902387
   ht   3    -0.00000000     0.00000000   -40.33935744

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3
 ref    1   -1.00000      -8.791772E-14   4.299459E-14
 ref    2   8.800701E-14   -1.00000      -3.831200E-13
 ref    3   4.301707E-14  -3.831200E-13    1.00000    

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
 ref:   2     0.00000000    -1.00000000    -0.00000000
 ref:   3     0.00000000    -0.00000000     1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1    -94.3985323836  0.0000E+00  3.1591E-01  9.7425E-01  1.0000E-04   
 mr-sdci #  1  2    -94.0535597859  1.4211E-14  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  1  3    -94.0138933504  3.5527E-14  0.0000E+00  9.9686E-01  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000416
time for cinew                         0.005578
time for eigenvalue solver             0.000108
time for vector access                 0.000002

          starting ci iteration   2

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1   -40.72399647
   ht   2    -0.00000000   -40.37902387
   ht   3    -0.00000000     0.00000000   -40.33935744
   ht   4     0.31590748     0.00000055     0.00584293    -4.61714739

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1  -0.961890       8.303031E-07  -8.424064E-03   0.273307    
 ref    2   7.110895E-07    1.00000       1.416819E-07  -5.309764E-07
 ref    3  -7.123604E-03  -1.335572E-07   0.999958       5.750260E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1   0.925283        1.00000       0.999987       7.472989E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1    -0.96188991     0.00000083    -0.00842406     0.27330721
 ref:   2     0.00000071     1.00000000     0.00000014    -0.00000053
 ref:   3    -0.00712360    -0.00000013     0.99995809     0.00575026

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  2  1    -94.6553953055  2.5686E-01  1.1691E-02  1.8501E-01  1.0000E-04   
 mr-sdci #  2  2    -94.0535597859  4.6896E-13  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  2  3    -94.0139532738  5.9923E-05  0.0000E+00  9.9643E-01  1.0000E-04   
 mr-sdci #  2  4    -91.2172653734 -2.5541E+00  0.0000E+00  1.2100E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000092
time for cinew                         0.006326
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   3

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1   -40.72399647
   ht   2    -0.00000000   -40.37902387
   ht   3    -0.00000000     0.00000000   -40.33935744
   ht   4     0.31590748     0.00000055     0.00584293    -4.61714739
   ht   5     0.20505567    -0.00000081    -0.05832076     0.12001351    -0.18715474

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.959822      -2.448682E-07   7.957518E-04   0.103397      -0.260866    
 ref    2   2.427895E-07   -1.00000      -3.892506E-06   1.972905E-06   8.154712E-07
 ref    3  -1.168346E-03   3.828097E-06  -0.999570      -2.764275E-02  -9.706839E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.921259        1.00000       0.999141       1.145505E-02   6.814515E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.95982161    -0.00000024     0.00079575     0.10339693    -0.26086572
 ref:   2     0.00000024    -1.00000000    -0.00000389     0.00000197     0.00000082
 ref:   3    -0.00116835     0.00000383    -0.99957005    -0.02764275    -0.00970684

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  3  1    -94.6652633153  9.8680E-03  1.0204E-03  5.3718E-02  1.0000E-04   
 mr-sdci #  3  2    -94.0535597859  1.0850E-11  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  3  3    -94.0158696144  1.9163E-03  0.0000E+00  9.9335E-01  1.0000E-04   
 mr-sdci #  3  4    -91.7756517867  5.5839E-01  0.0000E+00  1.4558E+00  1.0000E-04   
 mr-sdci #  3  5    -91.1996542809 -2.4151E+00  0.0000E+00  1.1847E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.005648
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration   4

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.72399647
   ht   2    -0.00000000   -40.37902387
   ht   3    -0.00000000     0.00000000   -40.33935744
   ht   4     0.31590748     0.00000055     0.00584293    -4.61714739
   ht   5     0.20505567    -0.00000081    -0.05832076     0.12001351    -0.18715474
   ht   6    -0.07200742     0.00000123     0.01922203    -0.01848233     0.00751339    -0.01644158

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.959182       7.917079E-07  -8.487621E-03   0.121627      -0.147130      -0.208466    
 ref    2  -8.277650E-08   -1.00000      -5.361441E-05   5.319449E-06  -5.649055E-07   2.268239E-06
 ref    3   3.821149E-03   5.283264E-05  -0.995578      -8.503526E-02   1.458122E-02  -3.695081E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.920044        1.00000       0.991248       2.202409E-02   2.185995E-02   4.482354E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95918160     0.00000079    -0.00848762     0.12162688    -0.14713033    -0.20846626
 ref:   2    -0.00000008    -1.00000000    -0.00005361     0.00000532    -0.00000056     0.00000227
 ref:   3     0.00382115     0.00005283    -0.99557844    -0.08503526     0.01458122    -0.03695081

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  4  1    -94.6660894366  8.2612E-04  1.0049E-04  1.6896E-02  1.0000E-04   
 mr-sdci #  4  2    -94.0535597860  1.2049E-10  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  4  3    -94.0322127863  1.6343E-02  0.0000E+00  9.5901E-01  1.0000E-04   
 mr-sdci #  4  4    -92.1425071171  3.6686E-01  0.0000E+00  1.1994E+00  1.0000E-04   
 mr-sdci #  4  5    -91.3136379139  1.1398E-01  0.0000E+00  1.4866E+00  1.0000E-04   
 mr-sdci #  4  6    -91.0394269867  3.7365E+01  0.0000E+00  1.3388E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000107
time for cinew                         0.003967
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration   5

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.72399647
   ht   2    -0.00000000   -40.37902387
   ht   3    -0.00000000     0.00000000   -40.33935744
   ht   4     0.31590748     0.00000055     0.00584293    -4.61714739
   ht   5     0.20505567    -0.00000081    -0.05832076     0.12001351    -0.18715474
   ht   6    -0.07200742     0.00000123     0.01922203    -0.01848233     0.00751339    -0.01644158
   ht   7     0.03485246    -0.00000064    -0.01074087    -0.01171774    -0.00267344     0.00080371    -0.00184303

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958798       1.865814E-02  -1.036739E-04  -0.123512       4.113203E-02   9.502585E-02  -0.233196    
 ref    2  -2.297724E-07  -5.613395E-03  -0.999984      -1.119198E-05  -5.452310E-07   8.603206E-07   2.567371E-06
 ref    3   5.652948E-03   0.985588      -5.534528E-03   0.162730       1.468131E-02  -2.209650E-02  -3.698724E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919326       0.971764       0.999999       4.173647E-02   1.907385E-03   9.518169E-03   5.574854E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95879838     0.01865814    -0.00010367    -0.12351227     0.04113203     0.09502585    -0.23319624
 ref:   2    -0.00000023    -0.00561339    -0.99998424    -0.00001119    -0.00000055     0.00000086     0.00000257
 ref:   3     0.00565295     0.98558834    -0.00553453     0.16273042     0.01468131    -0.02209650    -0.03698724

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  5  1    -94.6661859581  9.6522E-05  7.9529E-06  4.9159E-03  1.0000E-04   
 mr-sdci #  5  2    -94.0540489604  4.8917E-04  0.0000E+00  9.0996E-01  1.0000E-04   
 mr-sdci #  5  3    -94.0535597707  2.1347E-02  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  5  4    -92.7694200014  6.2691E-01  0.0000E+00  1.1113E+00  1.0000E-04   
 mr-sdci #  5  5    -91.5595773442  2.4594E-01  0.0000E+00  1.2502E+00  1.0000E-04   
 mr-sdci #  5  6    -91.2671822979  2.2776E-01  0.0000E+00  1.7364E+00  1.0000E-04   
 mr-sdci #  5  7    -90.9413908782  3.7267E+01  0.0000E+00  1.1970E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000069
time for cinew                         0.008521
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   6

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.72399647
   ht   2    -0.00000000   -40.37902387
   ht   3    -0.00000000     0.00000000   -40.33935744
   ht   4     0.31590748     0.00000055     0.00584293    -4.61714739
   ht   5     0.20505567    -0.00000081    -0.05832076     0.12001351    -0.18715474
   ht   6    -0.07200742     0.00000123     0.01922203    -0.01848233     0.00751339    -0.01644158
   ht   7     0.03485246    -0.00000064    -0.01074087    -0.01171774    -0.00267344     0.00080371    -0.00184303
   ht   8     0.00140252     0.00000038    -0.00014959    -0.00025657    -0.00069145     0.00002424    -0.00002116    -0.00012200

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958737      -2.519972E-02   1.388296E-06   0.103741       7.923179E-02  -8.119734E-02  -1.065355E-02   0.237574    
 ref    2  -3.052588E-07   1.134221E-04    1.00000       1.875898E-05  -3.160588E-06   1.380833E-06  -4.028981E-06  -1.890853E-06
 ref    3   6.621154E-03  -0.959920       1.143600E-04  -0.268506       3.393674E-02  -2.163653E-02   6.404716E-02   2.630772E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.922082        1.00000       8.285778E-02   7.429379E-03   7.061148E-03   4.215537E-03   5.713347E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873728    -0.02519972     0.00000139     0.10374127     0.07923179    -0.08119734    -0.01065355     0.23757393
 ref:   2    -0.00000031     0.00011342     0.99999999     0.00001876    -0.00000316     0.00000138    -0.00000403    -0.00000189
 ref:   3     0.00662115    -0.95992012     0.00011436    -0.26850611     0.03393674    -0.02163653     0.06404716     0.02630772

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  6  1    -94.6661940993  8.1412E-06  7.2805E-07  1.4348E-03  1.0000E-04   
 mr-sdci #  6  2    -94.1151365478  6.1088E-02  0.0000E+00  7.7027E-01  1.0000E-04   
 mr-sdci #  6  3    -94.0535597856  1.4906E-08  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  6  4    -92.9724489433  2.0303E-01  0.0000E+00  1.1107E+00  1.0000E-04   
 mr-sdci #  6  5    -91.8558139412  2.9624E-01  0.0000E+00  1.2513E+00  1.0000E-04   
 mr-sdci #  6  6    -91.5225409689  2.5536E-01  0.0000E+00  1.2015E+00  1.0000E-04   
 mr-sdci #  6  7    -90.9772778878  3.5887E-02  0.0000E+00  1.7491E+00  1.0000E-04   
 mr-sdci #  6  8    -90.9402783842  3.7266E+01  0.0000E+00  1.1753E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001211
time for cinew                         0.008636
time for eigenvalue solver             0.000134
time for vector access                 0.000001

          starting ci iteration   7

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165819
   ht   2    -0.00000000   -40.44060063
   ht   3     0.00000000     0.00000000   -40.37902387
   ht   4     0.00000000    -0.00000000     0.00000000   -39.29791303
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -38.18127803
   ht   6    -0.00214625     0.00147262     0.00000000    -0.00263442     0.00061233    -0.00001340

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958742      -2.309527E-02   1.052478E-06  -6.377640E-02  -0.113599       3.139419E-02
 ref    2   3.276970E-07   9.936491E-05    1.00000      -2.633504E-05   1.699836E-06  -7.421394E-07
 ref    3  -6.862374E-03  -0.934070       1.023459E-04   0.349276       7.197389E-03   1.586257E-02

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
 mr-sdci #  7  1    -94.6661948534  7.5419E-07  1.0061E-07  5.4190E-04  1.0000E-04   
 mr-sdci #  7  2    -94.1564883289  4.1352E-02  0.0000E+00  6.5973E-01  1.0000E-04   
 mr-sdci #  7  3    -94.0535597856  2.2006E-11  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  7  4    -93.1491192225  1.7667E-01  0.0000E+00  1.2451E+00  1.0000E-04   
 mr-sdci #  7  5    -92.4261603738  5.7035E-01  0.0000E+00  1.0931E+00  1.0000E-04   
 mr-sdci #  7  6    -91.5155980466 -6.9429E-03  0.0000E+00  1.3287E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000271
time for cinew                         0.007762
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   8

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165819
   ht   2    -0.00000000   -40.44060063
   ht   3     0.00000000     0.00000000   -40.37902387
   ht   4     0.00000000    -0.00000000     0.00000000   -39.29791303
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -38.18127803
   ht   6    -0.00214625     0.00147262     0.00000000    -0.00263442     0.00061233    -0.00001340
   ht   7    -0.00185615     0.00004657     0.00000012     0.00064076    -0.00093865    -0.00000046    -0.00000173

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958734      -2.634541E-02   1.045385E-06   5.211192E-02   0.104184       1.621001E-02  -7.766519E-02
 ref    2   3.333912E-07   9.659903E-05    1.00000       3.141994E-05   1.326393E-07   1.354588E-06  -1.279594E-06
 ref    3  -6.914583E-03  -0.914757       1.011371E-04  -0.395176      -4.108923E-02  -2.204919E-02   1.016906E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919219       0.837474        1.00000       0.158880       1.254272E-02   7.489314E-04   6.135291E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873395    -0.02634541     0.00000105     0.05211192     0.10418445     0.01621001    -0.07766519
 ref:   2     0.00000033     0.00009660     0.99999999     0.00003142     0.00000013     0.00000135    -0.00000128
 ref:   3    -0.00691458    -0.91475688     0.00010114    -0.39517600    -0.04108923    -0.02204919     0.01016906

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  8  1    -94.6661949432  8.9756E-08  1.5035E-08  1.9097E-04  1.0000E-04   
 mr-sdci #  8  2    -94.1771293397  2.0641E-02  0.0000E+00  5.9114E-01  1.0000E-04   
 mr-sdci #  8  3    -94.0535597856  8.1002E-13  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  8  4    -93.2546208513  1.0550E-01  0.0000E+00  1.2464E+00  1.0000E-04   
 mr-sdci #  8  5    -92.6282950920  2.0213E-01  0.0000E+00  9.4369E-01  1.0000E-04   
 mr-sdci #  8  6    -91.6050797610  8.9482E-02  0.0000E+00  1.4845E+00  1.0000E-04   
 mr-sdci #  8  7    -91.2653717978  2.8809E-01  0.0000E+00  1.3964E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000889
time for cinew                         0.007635
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   9

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165819
   ht   2    -0.00000000   -40.44060063
   ht   3     0.00000000     0.00000000   -40.37902387
   ht   4     0.00000000    -0.00000000     0.00000000   -39.29791303
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -38.18127803
   ht   6    -0.00214625     0.00147262     0.00000000    -0.00263442     0.00061233    -0.00001340
   ht   7    -0.00185615     0.00004657     0.00000012     0.00064076    -0.00093865    -0.00000046    -0.00000173
   ht   8     0.00126934    -0.00040157     0.00000001     0.00031163    -0.00063535     0.00000066     0.00000007    -0.00000038

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7         v:   8

   eig(s)   7.796556E-09   4.197405E-08   3.388639E-07    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1   2.798189E-05  -4.358548E-05   5.537252E-05   0.182114       1.499898E-03   0.215929       0.625324      -0.727445    
   x:   2  -8.648083E-06   4.605922E-07  -3.309004E-05   0.925138       7.760230E-03   2.445045E-02   0.134160       0.354206    
   x:   3   4.040607E-10   2.959072E-09  -2.463766E-10   8.391521E-03  -0.999965       5.196529E-05  -5.451010E-05   7.567451E-06
   x:   4   5.388932E-06   1.819530E-05   6.125553E-05   0.325231       2.739488E-03  -0.387349      -0.637308      -0.581392    
   x:   5  -1.612046E-05  -2.229946E-05  -8.487037E-06  -7.146951E-02  -6.691051E-04  -0.895956       0.429894       8.570242E-02
   x:   6   4.013493E-02  -4.166965E-02  -0.998325       1.586348E-12   8.871945E-13  -5.550219E-06  -1.473533E-05  -8.626540E-05
   x:   7   1.883403E-02   0.998984      -4.093999E-02  -5.512021E-13  -8.094101E-12  -4.073945E-06   4.770227E-05  -2.245041E-05
   x:   8   0.999017      -1.715936E-02   4.087896E-02   2.604527E-12  -3.897118E-12  -1.790460E-05  -6.286230E-06   3.184954E-05
 bummer (warning):overlap matrix: # small eigenvalues=1

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958737       2.429035E-02   1.325294E-06   2.976687E-03   0.116137       4.349900E-02  -2.322433E-02   7.439870E-02
 ref    2  -3.349859E-07  -9.884560E-05    1.00000       3.886498E-05   1.081935E-05  -4.871071E-06  -9.575685E-07   3.445597E-06
 ref    3   6.919202E-03   0.900292       1.062529E-04  -0.365642      -0.203820       8.829989E-02   1.371619E-02  -5.408355E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919224       0.811116        1.00000       0.133703       5.503022E-02   9.689033E-03   7.275034E-04   8.460197E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873651     0.02429035     0.00000133     0.00297669     0.11613703     0.04349900    -0.02322433     0.07439870
 ref:   2    -0.00000033    -0.00009885     0.99999999     0.00003886     0.00001082    -0.00000487    -0.00000096     0.00000345
 ref:   3     0.00691920     0.90029234     0.00010625    -0.36564210    -0.20381956     0.08829989     0.01371619    -0.05408355

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1    -94.6661949572  1.3955E-08  0.0000E+00  7.5894E-05  1.0000E-04   
 mr-sdci #  9  2    -94.1826832566  5.5539E-03  1.0845E-01  5.5580E-01  1.0000E-04   
 mr-sdci #  9  3    -94.0535597857  4.2938E-11  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci #  9  4    -93.5649404341  3.1032E-01  0.0000E+00  1.0685E+00  1.0000E-04   
 mr-sdci #  9  5    -92.8427163686  2.1442E-01  0.0000E+00  9.0684E-01  1.0000E-04   
 mr-sdci #  9  6    -91.8522841914  2.4720E-01  0.0000E+00  1.2677E+00  1.0000E-04   
 mr-sdci #  9  7    -91.6016174393  3.3625E-01  0.0000E+00  1.4960E+00  1.0000E-04   
 mr-sdci #  9  8    -91.0386186094  9.8340E-02  0.0000E+00  1.5120E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000011
time for cinew                         0.011503
time for eigenvalue solver             0.000151
time for vector access                 0.000002

          starting ci iteration  10

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165904
   ht   2     0.00000000   -40.50814734
   ht   3    -0.00000000    -0.00000000   -40.37902387
   ht   4    -0.00000000     0.00000000     0.00000000   -39.89040452
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -39.16818046
   ht   6     0.03595944    -0.97867152    -0.00001350     0.20038725     0.73902141    -1.91359576

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958737       1.457480E-02   2.199249E-06  -2.134044E-02  -0.111410       3.172969E-02
 ref    2   3.349813E-07  -1.215813E-04    1.00000      -8.373033E-06   1.237942E-05   3.371603E-05
 ref    3  -6.919169E-03   0.932251       1.246879E-04   0.204014       3.131977E-02  -0.272155    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919224       0.869304        1.00000       4.207724E-02   1.339307E-02   7.507500E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873651     0.01457480     0.00000220    -0.02134044    -0.11140981     0.03172969
 ref:   2     0.00000033    -0.00012158     0.99999999    -0.00000837     0.00001238     0.00003372
 ref:   3    -0.00691917     0.93225104     0.00012469     0.20401427     0.03131977    -0.27215477

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 10  1    -94.6661949572 -4.2633E-14  0.0000E+00  7.5896E-05  1.0000E-04   
 mr-sdci # 10  2    -94.2922545635  1.0957E-01  1.4959E-02  2.0439E-01  1.0000E-04   
 mr-sdci # 10  3    -94.0535597858  1.3808E-10  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 10  4    -93.6606547702  9.5714E-02  0.0000E+00  8.8471E-01  1.0000E-04   
 mr-sdci # 10  5    -92.9725941535  1.2988E-01  0.0000E+00  7.4497E-01  1.0000E-04   
 mr-sdci # 10  6    -91.3335294986 -5.1875E-01  0.0000E+00  1.4434E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000143
time for cinew                         0.004763
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  11

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165904
   ht   2     0.00000000   -40.50814734
   ht   3    -0.00000000    -0.00000000   -40.37902387
   ht   4    -0.00000000     0.00000000     0.00000000   -39.89040452
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -39.16818046
   ht   6     0.03595944    -0.97867152    -0.00001350     0.20038725     0.73902141    -1.91359576
   ht   7    -4.70204544     1.71566224     0.00012316     0.36479757     0.76475672     0.12388410    -0.93621090

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958736       6.282649E-03   7.418913E-06  -4.562179E-02   0.114814      -6.757569E-02   3.326056E-02
 ref    2   3.332901E-07  -1.414540E-04    1.00000       6.766518E-05  -1.668965E-05  -4.848843E-06  -4.406599E-05
 ref    3  -6.913575E-03   0.943991       1.362675E-04   0.130606      -1.806976E-02   0.151091       0.234127    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919222       0.891159        1.00000       1.913939E-02   1.350873E-02   2.739504E-02   5.592195E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873555     0.00628265     0.00000742    -0.04562179     0.11481384    -0.06757569     0.03326056
 ref:   2     0.00000033    -0.00014145     0.99999999     0.00006767    -0.00001669    -0.00000485    -0.00004407
 ref:   3    -0.00691357     0.94399123     0.00013627     0.13060644    -0.01806976     0.15109124     0.23412749

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 11  1    -94.6661949574  2.6588E-10  0.0000E+00  7.1179E-05  1.0000E-04   
 mr-sdci # 11  2    -94.3044368916  1.2182E-02 -6.7573E-03  1.2116E-01  1.0000E-04   
 mr-sdci # 11  3    -94.0535597884  2.5803E-09  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 11  4    -93.8066333388  1.4598E-01  0.0000E+00  5.6055E-01  1.0000E-04   
 mr-sdci # 11  5    -92.9783901407  5.7960E-03  0.0000E+00  7.5773E-01  1.0000E-04   
 mr-sdci # 11  6    -91.4671825484  1.3365E-01  0.0000E+00  1.3749E+00  1.0000E-04   
 mr-sdci # 11  7    -91.1415865682 -4.6003E-01  0.0000E+00  1.5315E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000287
time for cinew                         0.003955
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  12

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165904
   ht   2     0.00000000   -40.50814734
   ht   3    -0.00000000    -0.00000000   -40.37902387
   ht   4    -0.00000000     0.00000000     0.00000000   -39.89040452
   ht   5     0.00000000     0.00000000     0.00000000    -0.00000000   -39.16818046
   ht   6     0.03595944    -0.97867152    -0.00001350     0.20038725     0.73902141    -1.91359576
   ht   7    -4.70204544     1.71566224     0.00012316     0.36479757     0.76475672     0.12388410    -0.93621090
   ht   8   -17.45416658     6.48495678     0.00024823    -2.95119357     0.65279559     0.17592129    -2.35972777    -9.02269964

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735       1.179314E-02  -1.009892E-05   3.568726E-02  -9.505869E-02   0.113660      -0.192264       0.119677    
 ref    2   3.331772E-07  -1.372002E-04   -1.00000      -6.153480E-05   1.356884E-05   3.342523E-06   1.412339E-06  -4.739465E-05
 ref    3  -6.914418E-03   0.946020      -1.339397E-04  -0.126382       1.630184E-02  -0.154368       2.442691E-02   0.234457    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919220       0.895093        1.00000       1.724588E-02   9.301905E-03   3.674826E-02   3.756215E-02   6.929248E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873489     0.01179314    -0.00001010     0.03568726    -0.09505869     0.11366035    -0.19226408     0.11967712
 ref:   2     0.00000033    -0.00013720    -0.99999999    -0.00006153     0.00001357     0.00000334     0.00000141    -0.00004739
 ref:   3    -0.00691442     0.94602020    -0.00013394    -0.12638154     0.01630184    -0.15436835     0.02442691     0.23445653

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 12  1    -94.6661949575  2.9992E-11  0.0000E+00  7.1715E-05  1.0000E-04   
 mr-sdci # 12  2    -94.3064607382  2.0238E-03  2.2964E-03  8.6683E-02  1.0000E-04   
 mr-sdci # 12  3    -94.0535597889  5.0731E-10  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 12  4    -93.8119807509  5.3474E-03  0.0000E+00  5.2755E-01  1.0000E-04   
 mr-sdci # 12  5    -92.9894859946  1.1096E-02  0.0000E+00  7.7358E-01  1.0000E-04   
 mr-sdci # 12  6    -91.4720452777  4.8627E-03  0.0000E+00  1.3189E+00  1.0000E-04   
 mr-sdci # 12  7    -91.3671223411  2.2554E-01  0.0000E+00  1.4918E+00  1.0000E-04   
 mr-sdci # 12  8    -91.1038920592  6.5273E-02  0.0000E+00  1.4358E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000002
time for cinew                         0.007695
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  13

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165904
   ht   2    -0.00000000   -40.63192483
   ht   3    -0.00000000     0.00000000   -40.37902388
   ht   4     0.00000000     0.00000000     0.00000000   -40.13744484
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -39.31495008
   ht   6     3.98141477     2.15018705    -0.00018075     0.16716142    -0.31494504    -0.58306816

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735       9.088823E-03   1.748933E-05   4.693497E-02   9.770378E-02  -0.117128    
 ref    2  -3.328441E-07  -1.442863E-04    1.00000      -1.074518E-04  -1.631804E-05   6.775156E-05
 ref    3   6.914368E-03   0.944103       1.444451E-04  -0.104677      -1.152158E-02  -0.167480    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919220       0.891413        1.00000       1.316014E-02   9.678776E-03   4.176863E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873460     0.00908882     0.00001749     0.04693497     0.09770378    -0.11712783
 ref:   2    -0.00000033    -0.00014429     0.99999998    -0.00010745    -0.00001632     0.00006775
 ref:   3     0.00691437     0.94410287     0.00014445    -0.10467682    -0.01152158    -0.16748043

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 13  1    -94.6661949575  1.4481E-11  0.0000E+00  7.2239E-05  1.0000E-04   
 mr-sdci # 13  2    -94.3074343240  9.7359E-04 -2.4927E-03  7.4655E-02  1.0000E-04   
 mr-sdci # 13  3    -94.0535597939  5.0077E-09  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 13  4    -93.8295618019  1.7581E-02  0.0000E+00  5.1620E-01  1.0000E-04   
 mr-sdci # 13  5    -92.9902458256  7.5983E-04  0.0000E+00  7.7055E-01  1.0000E-04   
 mr-sdci # 13  6    -91.6182616055  1.4622E-01  0.0000E+00  1.4403E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000391
time for cinew                         0.003039
time for eigenvalue solver             0.000113
time for vector access                 0.000000

          starting ci iteration  14

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165904
   ht   2    -0.00000000   -40.63192483
   ht   3    -0.00000000     0.00000000   -40.37902388
   ht   4     0.00000000     0.00000000     0.00000000   -40.13744484
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -39.31495008
   ht   6     3.98141477     2.15018705    -0.00018075     0.16716142    -0.31494504    -0.58306816
   ht   7    11.08535943     3.68222791    -0.00043767    -0.52357742    -1.05986855    -1.35653228    -3.52465816

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735       1.239689E-02   1.226600E-05   3.119626E-02  -9.428456E-02   7.341483E-02  -0.237773    
 ref    2  -3.326462E-07  -1.497939E-04    1.00000      -1.390996E-04   1.704648E-05  -6.813955E-05   9.301213E-06
 ref    3   6.914147E-03   0.944831       1.476190E-04  -9.638999E-02   1.051424E-02   0.167659      -2.005951E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919220       0.892860        1.00000       1.026426E-02   9.000127E-03   3.349939E-02   5.693833E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873486     0.01239689     0.00001227     0.03119626    -0.09428456     0.07341483    -0.23777288
 ref:   2    -0.00000033    -0.00014979     0.99999998    -0.00013910     0.00001705    -0.00006814     0.00000930
 ref:   3     0.00691415     0.94483114     0.00014762    -0.09638999     0.01051424     0.16765933    -0.02005951

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 14  1    -94.6661949575  5.5707E-12  0.0000E+00  7.2157E-05  1.0000E-04   
 mr-sdci # 14  2    -94.3083823491  9.4803E-04  8.7790E-04  5.0691E-02  1.0000E-04   
 mr-sdci # 14  3    -94.0535597967  2.8159E-09  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 14  4    -93.8452394613  1.5678E-02  0.0000E+00  4.7442E-01  1.0000E-04   
 mr-sdci # 14  5    -92.9907095713  4.6375E-04  0.0000E+00  7.7057E-01  1.0000E-04   
 mr-sdci # 14  6    -91.6373017091  1.9040E-02  0.0000E+00  1.4779E+00  1.0000E-04   
 mr-sdci # 14  7    -91.1007697330 -2.6635E-01  0.0000E+00  1.2969E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000874
time for cinew                         0.005558
time for eigenvalue solver             0.000083
time for vector access                 0.000001

          starting ci iteration  15

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165904
   ht   2    -0.00000000   -40.63192483
   ht   3    -0.00000000     0.00000000   -40.37902388
   ht   4     0.00000000     0.00000000     0.00000000   -40.13744484
   ht   5    -0.00000000     0.00000000     0.00000000     0.00000000   -39.31495008
   ht   6     3.98141477     2.15018705    -0.00018075     0.16716142    -0.31494504    -0.58306816
   ht   7    11.08535943     3.68222791    -0.00043767    -0.52357742    -1.05986855    -1.35653228    -3.52465816
   ht   8     1.96716668     0.08480478    -0.00003295    -0.13239136    -0.27504638    -0.20232330    -0.58122296    -0.12123069

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958735      -1.225101E-02  -1.724974E-05  -3.091244E-02  -9.427903E-02   4.835525E-02   2.833914E-02   0.243724    
 ref    2  -3.277684E-07   1.742386E-04   -1.00000       3.280535E-04   1.648952E-05  -1.653796E-04  -3.041223E-05   1.060037E-05
 ref    3   6.912325E-03  -0.946654      -1.484253E-04   9.280848E-02   1.053572E-02   7.218117E-02  -0.151133       8.281522E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919220       0.896304        1.00000       9.569102E-03   8.999537E-03   7.548378E-03   2.364427E-02   6.625986E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873487    -0.01225101    -0.00001725    -0.03091244    -0.09427903     0.04835525     0.02833914     0.24372423
 ref:   2    -0.00000033     0.00017424    -0.99999992     0.00032805     0.00001649    -0.00016538    -0.00003041     0.00001060
 ref:   3     0.00691232    -0.94665413    -0.00014843     0.09280848     0.01053572     0.07218117    -0.15113292     0.08281522

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 15  1    -94.6661949577  2.0832E-10  0.0000E+00  6.5703E-05  1.0000E-04   
 mr-sdci # 15  2    -94.3091893658  8.0702E-04  1.4552E-04  1.7001E-02  1.0000E-04   
 mr-sdci # 15  3    -94.0535598595  6.2823E-08  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 15  4    -93.8652964778  2.0057E-02  0.0000E+00  4.1749E-01  1.0000E-04   
 mr-sdci # 15  5    -92.9907130757  3.5044E-06  0.0000E+00  7.7001E-01  1.0000E-04   
 mr-sdci # 15  6    -92.0625884496  4.2529E-01  0.0000E+00  1.3261E+00  1.0000E-04   
 mr-sdci # 15  7    -91.2485770983  1.4781E-01  0.0000E+00  1.7314E+00  1.0000E-04   
 mr-sdci # 15  8    -91.0674814652 -3.6411E-02  0.0000E+00  1.2631E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001053
time for cinew                         0.008662
time for eigenvalue solver             0.000116
time for vector access                 0.000000

          starting ci iteration  16

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165904
   ht   2     0.00000000   -40.63465345
   ht   3     0.00000000    -0.00000000   -40.37902395
   ht   4     0.00000000     0.00000000    -0.00000000   -40.19076056
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -39.31617716
   ht   6    -0.03073926     0.29533508     0.00006911     0.07099983    -0.00500775    -0.00615177

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735       1.241830E-02  -1.482957E-05  -2.649914E-02  -9.473912E-02  -3.847651E-02
 ref    2   3.235664E-07  -1.983121E-04  -0.999999       9.549030E-04   7.625354E-06  -4.802258E-04
 ref    3  -6.912722E-03   0.945551      -2.294545E-04   6.669005E-02   1.293177E-02   0.206006    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.894221       0.999999       5.150679E-03   9.142732E-03   4.391911E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873496     0.01241830    -0.00001483    -0.02649914    -0.09473912    -0.03847651
 ref:   2     0.00000032    -0.00019831    -0.99999940     0.00095490     0.00000763    -0.00048023
 ref:   3    -0.00691272     0.94555125    -0.00022945     0.06669005     0.01293177     0.20600592

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 16  1    -94.6661949577  2.2929E-11  0.0000E+00  6.4858E-05  1.0000E-04   
 mr-sdci # 16  2    -94.3092936699  1.0430E-04  3.3493E-05  8.9416E-03  1.0000E-04   
 mr-sdci # 16  3    -94.0535603541  4.9453E-07  0.0000E+00  9.9729E-01  1.0000E-04   
 mr-sdci # 16  4    -93.8877050801  2.2409E-02  0.0000E+00  3.6128E-01  1.0000E-04   
 mr-sdci # 16  5    -92.9908112742  9.8198E-05  0.0000E+00  7.6809E-01  1.0000E-04   
 mr-sdci # 16  6    -92.1732228563  1.1063E-01  0.0000E+00  1.5805E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000080
time for cinew                         0.003001
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  17

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165904
   ht   2     0.00000000   -40.63465345
   ht   3     0.00000000    -0.00000000   -40.37902395
   ht   4     0.00000000     0.00000000    -0.00000000   -40.19076056
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -39.31617716
   ht   6    -0.03073926     0.29533508     0.00006911     0.07099983    -0.00500775    -0.00615177
   ht   7    -0.23039728     0.09867429     0.00008807     0.02928266    -0.03168382    -0.00091722    -0.00242979

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958735      -1.207729E-02   1.681776E-04  -3.236530E-02  -2.831181E-02  -9.678719E-02  -9.583179E-02
 ref    2   3.167571E-07   2.331911E-04   0.999996       2.192102E-03   1.978695E-03   1.422892E-04   2.830725E-04
 ref    3  -6.912579E-03  -0.945733       1.737368E-04   6.574631E-02  -7.414748E-02   6.348462E-03   0.228892    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919220       0.894556       0.999991       5.374895E-03   6.303323E-03   9.408083E-03   6.157517E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873480    -0.01207729     0.00016818    -0.03236530    -0.02831181    -0.09678719    -0.09583179
 ref:   2     0.00000032     0.00023319     0.99999555     0.00219210     0.00197869     0.00014229     0.00028307
 ref:   3    -0.00691258    -0.94573255     0.00017374     0.06574631    -0.07414748     0.00634846     0.22889158

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 17  1    -94.6661949577  8.6189E-12  0.0000E+00  6.4459E-05  1.0000E-04   
 mr-sdci # 17  2    -94.3093224801  2.8810E-05 -4.4624E-06  5.9592E-03  1.0000E-04   
 mr-sdci # 17  3    -94.0535642521  3.8980E-06  0.0000E+00  9.9728E-01  1.0000E-04   
 mr-sdci # 17  4    -93.8962049257  8.4998E-03  0.0000E+00  3.5034E-01  1.0000E-04   
 mr-sdci # 17  5    -93.0933973967  1.0259E-01  0.0000E+00  1.2375E+00  1.0000E-04   
 mr-sdci # 17  6    -92.9902451535  8.1702E-01  0.0000E+00  7.8118E-01  1.0000E-04   
 mr-sdci # 17  7    -91.1513763692 -9.7201E-02  0.0000E+00  1.4124E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000322
time for cinew                         0.005766
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  18

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165904
   ht   2     0.00000000   -40.63465345
   ht   3     0.00000000    -0.00000000   -40.37902395
   ht   4     0.00000000     0.00000000    -0.00000000   -40.19076056
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -39.31617716
   ht   6    -0.03073926     0.29533508     0.00006911     0.07099983    -0.00500775    -0.00615177
   ht   7    -0.23039728     0.09867429     0.00008807     0.02928266    -0.03168382    -0.00091722    -0.00242979
   ht   8     0.71919580     0.27292029     0.00001502    -0.03326016     0.06483247    -0.00166309     0.00366993    -0.01541456

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735      -1.214166E-02   3.111132E-04   3.256018E-02   4.735854E-02   9.794786E-02  -0.191879       9.658533E-02
 ref    2   3.171313E-07   2.308305E-04   0.999995      -2.204545E-03  -2.166197E-03  -1.306538E-04   5.534724E-04  -2.502170E-04
 ref    3  -6.912560E-03  -0.945782       2.805305E-04  -6.559676E-02   8.648935E-02  -6.246139E-03   9.747284E-02   0.252118    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919220       0.894650       0.999990       5.367961E-03   9.727931E-03   9.632815E-03   4.631896E-02   7.289221E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873482    -0.01214166     0.00031111     0.03256018     0.04735854     0.09794786    -0.19187938     0.09658533
 ref:   2     0.00000032     0.00023083     0.99999499    -0.00220455    -0.00216620    -0.00013065     0.00055347    -0.00025022
 ref:   3    -0.00691256    -0.94578164     0.00028053    -0.06559676     0.08648935    -0.00624614     0.09747284     0.25211787

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 18  1    -94.6661949577  3.5527E-14  0.0000E+00  6.4437E-05  1.0000E-04   
 mr-sdci # 18  2    -94.3093227933  3.1314E-07  3.6802E-06  5.8357E-03  1.0000E-04   
 mr-sdci # 18  3    -94.0535656419  1.3899E-06  0.0000E+00  9.9727E-01  1.0000E-04   
 mr-sdci # 18  4    -93.8962073417  2.4160E-06  0.0000E+00  3.5048E-01  1.0000E-04   
 mr-sdci # 18  5    -93.1094527767  1.6055E-02  0.0000E+00  1.2059E+00  1.0000E-04   
 mr-sdci # 18  6    -92.9903756113  1.3046E-04  0.0000E+00  7.7689E-01  1.0000E-04   
 mr-sdci # 18  7    -91.2044424993  5.3066E-02  0.0000E+00  1.3662E+00  1.0000E-04   
 mr-sdci # 18  8    -91.0567831709 -1.0698E-02  0.0000E+00  1.4461E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000791
time for cinew                         0.008298
time for eigenvalue solver             0.000139
time for vector access                 0.000000

          starting ci iteration  19

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165904
   ht   2     0.00000000   -40.63478688
   ht   3    -0.00000000     0.00000000   -40.37902973
   ht   4     0.00000000    -0.00000000    -0.00000000   -40.22167143
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000   -39.43491686
   ht   6    -0.52495933     0.21892818    -0.00019802     0.01084142    -0.05650427    -0.00855611

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735      -1.207950E-02  -3.716449E-04  -3.419128E-02  -6.223721E-02   0.183654    
 ref    2  -3.170231E-07   2.315585E-04  -0.999995       2.255599E-03   2.268778E-03  -6.490437E-04
 ref    3   6.912587E-03  -0.945737      -3.229312E-04   6.434989E-02  -9.701001E-02   0.125591    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919220       0.894565       0.999989       5.315040E-03   1.328956E-02   4.950223E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873478    -0.01207950    -0.00037164    -0.03419128    -0.06223721     0.18365408
 ref:   2    -0.00000032     0.00023156    -0.99999459     0.00225560     0.00226878    -0.00064904
 ref:   3     0.00691259    -0.94573737    -0.00032293     0.06434989    -0.09701001     0.12559056

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 19  1    -94.6661949577  1.3500E-13  0.0000E+00  6.4438E-05  1.0000E-04   
 mr-sdci # 19  2    -94.3093230936  3.0034E-07 -4.1749E-06  5.7490E-03  1.0000E-04   
 mr-sdci # 19  3    -94.0535658906  2.4867E-07  0.0000E+00  9.9727E-01  1.0000E-04   
 mr-sdci # 19  4    -93.8963827185  1.7538E-04  0.0000E+00  3.5032E-01  1.0000E-04   
 mr-sdci # 19  5    -93.1211653686  1.1713E-02  0.0000E+00  1.1812E+00  1.0000E-04   
 mr-sdci # 19  6    -91.3001385627 -1.6902E+00  0.0000E+00  1.4975E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000324
time for cinew                         0.006450
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  20

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165904
   ht   2     0.00000000   -40.63478688
   ht   3    -0.00000000     0.00000000   -40.37902973
   ht   4     0.00000000    -0.00000000    -0.00000000   -40.22167143
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000   -39.43491686
   ht   6    -0.52495933     0.21892818    -0.00019802     0.01084142    -0.05650427    -0.00855611
   ht   7    -0.69580048     0.26447675    -0.00029748     0.02201839    -0.07068338    -0.01109185    -0.01443826

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735      -1.237479E-02  -8.372052E-05  -2.570036E-02  -1.214104E-02   0.261178       5.778873E-02
 ref    2  -3.076638E-07   2.638068E-04  -0.999982       4.189500E-03  -4.220250E-03  -2.366610E-04   7.881689E-04
 ref    3   6.912704E-03  -0.945543      -5.577009E-04   5.496697E-02   0.114794       8.943151E-02  -9.794857E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.894205       0.999964       3.699428E-03   1.334291E-02   7.621198E-02   1.293408E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873502    -0.01237479    -0.00008372    -0.02570036    -0.01214104     0.26117795     0.05778873
 ref:   2    -0.00000031     0.00026381    -0.99998190     0.00418950    -0.00422025    -0.00023666     0.00078817
 ref:   3     0.00691270    -0.94554309    -0.00055770     0.05496697     0.11479415     0.08943151    -0.09794857

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 20  1    -94.6661949577  6.2883E-12  0.0000E+00  6.4238E-05  1.0000E-04   
 mr-sdci # 20  2    -94.3093326187  9.5251E-06  3.9097E-06  2.8157E-03  1.0000E-04   
 mr-sdci # 20  3    -94.0535765280  1.0637E-05  0.0000E+00  9.9723E-01  1.0000E-04   
 mr-sdci # 20  4    -93.9029190860  6.5364E-03  0.0000E+00  3.4266E-01  1.0000E-04   
 mr-sdci # 20  5    -93.3783986870  2.5723E-01  0.0000E+00  8.4576E-01  1.0000E-04   
 mr-sdci # 20  6    -91.4258579890  1.2572E-01  0.0000E+00  1.3621E+00  1.0000E-04   
 mr-sdci # 20  7    -91.0211045816 -1.8334E-01  0.0000E+00  1.5686E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000175
time for cinew                         0.005863
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  21

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165904
   ht   2     0.00000000   -40.63478688
   ht   3    -0.00000000     0.00000000   -40.37902973
   ht   4     0.00000000    -0.00000000    -0.00000000   -40.22167143
   ht   5    -0.00000000    -0.00000000    -0.00000000     0.00000000   -39.43491686
   ht   6    -0.52495933     0.21892818    -0.00019802     0.01084142    -0.05650427    -0.00855611
   ht   7    -0.69580048     0.26447675    -0.00029748     0.02201839    -0.07068338    -0.01109185    -0.01443826
   ht   8    -0.07037834    -0.03233324     0.00001474     0.00898103    -0.02082089    -0.00078144    -0.00107069    -0.00027552

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958735      -1.232649E-02  -3.836644E-04   2.767813E-02  -2.262897E-03   0.106510      -0.235417      -7.101344E-02
 ref    2  -3.291299E-07   3.075579E-04  -0.999900      -9.117228E-03   1.045614E-02  -2.567189E-03  -4.431988E-04  -9.922692E-04
 ref    3   6.912912E-03  -0.945680       3.658616E-05  -5.730933E-02  -5.828896E-02  -0.156830      -0.152344       7.265452E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.894463       0.999800       4.133562E-03   3.512055E-03   3.594670E-02   7.863000E-02   1.032257E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873508    -0.01232649    -0.00038366     0.02767813    -0.00226290     0.10651027    -0.23541720    -0.07101344
 ref:   2    -0.00000033     0.00030756    -0.99989980    -0.00911723     0.01045614    -0.00256719    -0.00044320    -0.00099227
 ref:   3     0.00691291    -0.94568034     0.00003659    -0.05730933    -0.05828896    -0.15683006    -0.15234351     0.07265452

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 21  1    -94.6661949577  6.3878E-12  0.0000E+00  6.4648E-05  1.0000E-04   
 mr-sdci # 21  2    -94.3093358183  3.1996E-06  9.4699E-07  1.4857E-03  1.0000E-04   
 mr-sdci # 21  3    -94.0536381446  6.1617E-05  0.0000E+00  9.9703E-01  1.0000E-04   
 mr-sdci # 21  4    -93.9091614115  6.2423E-03  0.0000E+00  3.2262E-01  1.0000E-04   
 mr-sdci # 21  5    -93.5784930135  2.0009E-01  0.0000E+00  7.5280E-01  1.0000E-04   
 mr-sdci # 21  6    -92.0501741363  6.2432E-01  0.0000E+00  1.4357E+00  1.0000E-04   
 mr-sdci # 21  7    -91.3478545796  3.2675E-01  0.0000E+00  1.3708E+00  1.0000E-04   
 mr-sdci # 21  8    -91.0111711290 -4.5612E-02  0.0000E+00  1.5401E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001322
time for cinew                         0.010361
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  22

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165904
   ht   2    -0.00000000   -40.63479991
   ht   3    -0.00000000     0.00000000   -40.37910223
   ht   4     0.00000000    -0.00000000     0.00000000   -40.23462550
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -39.90395710
   ht   6    -0.00244935     0.02171821     0.00019589    -0.00505735     0.00942649    -0.00003504

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735      -1.234745E-02  -1.522028E-05  -2.559659E-02  -1.157455E-02  -3.807254E-02
 ref    2   5.072755E-07   3.650819E-04   0.999275       1.984045E-02  -3.170954E-02  -6.519235E-03
 ref    3  -6.912018E-03  -0.945557       3.359774E-03   4.282323E-02   0.105592       0.132487    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919220       0.894230       0.998562       2.882658E-03   1.228917E-02   1.904495E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873490    -0.01234745    -0.00001522    -0.02559659    -0.01157455    -0.03807254
 ref:   2     0.00000051     0.00036508     0.99927502     0.01984045    -0.03170954    -0.00651923
 ref:   3    -0.00691202    -0.94555687     0.00335977     0.04282323     0.10559217     0.13248747

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 22  1    -94.6661949578  6.4567E-11  0.0000E+00  6.5503E-05  1.0000E-04   
 mr-sdci # 22  2    -94.3093366584  8.4018E-07  3.3709E-07  8.3389E-04  1.0000E-04   
 mr-sdci # 22  3    -94.0540558907  4.1775E-04  0.0000E+00  9.9556E-01  1.0000E-04   
 mr-sdci # 22  4    -93.9125909736  3.4296E-03  0.0000E+00  3.0668E-01  1.0000E-04   
 mr-sdci # 22  5    -93.7449168842  1.6642E-01  0.0000E+00  5.6554E-01  1.0000E-04   
 mr-sdci # 22  6    -91.6990135170 -3.5116E-01  0.0000E+00  1.5961E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000208
time for cinew                         0.004665
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  23

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165904
   ht   2    -0.00000000   -40.63479991
   ht   3    -0.00000000     0.00000000   -40.37910223
   ht   4     0.00000000    -0.00000000     0.00000000   -40.23462550
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -39.90395710
   ht   6    -0.00244935     0.02171821     0.00019589    -0.00505735     0.00942649    -0.00003504
   ht   7     0.03047704    -0.00979617    -0.00017290     0.00382732     0.00368574     0.00000627    -0.00003628

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958736       1.229912E-02   4.491638E-03  -2.675122E-02   1.846206E-02   7.989761E-02  -0.107894    
 ref    2   9.418713E-07  -4.332464E-04   0.995678       2.467434E-02  -8.543579E-02  -2.659492E-02   2.604957E-03
 ref    3  -6.912717E-03   0.945585       2.800870E-03   4.178934E-02   7.211476E-02  -8.231834E-02   0.213949    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919222       0.894283       0.991402       3.070799E-03   1.284066E-02   1.386723E-02   5.742207E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873574     0.01229912     0.00449164    -0.02675122     0.01846206     0.07989761    -0.10789446
 ref:   2     0.00000094    -0.00043325     0.99567771     0.02467434    -0.08543579    -0.02659492     0.00260496
 ref:   3    -0.00691272     0.94558536     0.00280087     0.04178934     0.07211476    -0.08231834     0.21394876

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 23  1    -94.6661949579  1.0452E-10  0.0000E+00  6.4191E-05  1.0000E-04   
 mr-sdci # 23  2    -94.3093369408  2.8238E-07 -2.2138E-07  6.6121E-04  1.0000E-04   
 mr-sdci # 23  3    -94.0560152348  1.9593E-03  0.0000E+00  9.8764E-01  1.0000E-04   
 mr-sdci # 23  4    -93.9126921884  1.0121E-04  0.0000E+00  3.0756E-01  1.0000E-04   
 mr-sdci # 23  5    -93.8314462817  8.6529E-02  0.0000E+00  5.3360E-01  1.0000E-04   
 mr-sdci # 23  6    -92.9174803183  1.2185E+00  0.0000E+00  9.5730E-01  1.0000E-04   
 mr-sdci # 23  7    -91.0185049746 -3.2935E-01  0.0000E+00  1.3906E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000035
time for cinew                         0.004990
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  24

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165904
   ht   2    -0.00000000   -40.63479991
   ht   3    -0.00000000     0.00000000   -40.37910223
   ht   4     0.00000000    -0.00000000     0.00000000   -40.23462550
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -39.90395710
   ht   6    -0.00244935     0.02171821     0.00019589    -0.00505735     0.00942649    -0.00003504
   ht   7     0.03047704    -0.00979617    -0.00017290     0.00382732     0.00368574     0.00000627    -0.00003628
   ht   8    -0.10661740    -0.03592520     0.00004472    -0.00592830    -0.00367311     0.00001427     0.00007648    -0.00032747

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958736      -1.232467E-02   4.515318E-03  -2.528120E-02   5.964496E-03  -5.338963E-02  -0.246753      -4.110719E-02
 ref    2   9.215256E-07   4.293429E-04   0.995685       2.592590E-02  -8.514306E-02   2.483069E-02   8.305740E-03   7.149206E-04
 ref    3  -6.912895E-03  -0.945605       2.817290E-03   4.221534E-02   6.224609E-02   0.100965      -9.716893E-02   0.251619    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919222       0.894320       0.991417       3.093426E-03   1.115949E-02   1.366092E-02   7.039777E-02   6.500260E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873551    -0.01232467     0.00451532    -0.02528120     0.00596450    -0.05338963    -0.24675288    -0.04110719
 ref:   2     0.00000092     0.00042934     0.99568516     0.02592590    -0.08514306     0.02483069     0.00830574     0.00071492
 ref:   3    -0.00691289    -0.94560469     0.00281729     0.04221534     0.06224609     0.10096486    -0.09716893     0.25161934

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 24  1    -94.6661949579  3.8582E-12  0.0000E+00  6.4129E-05  1.0000E-04   
 mr-sdci # 24  2    -94.3093369851  4.4310E-08  9.8716E-08  5.6017E-04  1.0000E-04   
 mr-sdci # 24  3    -94.0560152724  3.7644E-08  0.0000E+00  9.8765E-01  1.0000E-04   
 mr-sdci # 24  4    -93.9128475052  1.5532E-04  0.0000E+00  3.0792E-01  1.0000E-04   
 mr-sdci # 24  5    -93.8396669724  8.2207E-03  0.0000E+00  5.2121E-01  1.0000E-04   
 mr-sdci # 24  6    -92.9378707263  2.0390E-02  0.0000E+00  9.1347E-01  1.0000E-04   
 mr-sdci # 24  7    -91.4225610225  4.0406E-01  0.0000E+00  1.3954E+00  1.0000E-04   
 mr-sdci # 24  8    -90.9842994953 -2.6872E-02  0.0000E+00  1.3930E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000171
time for cinew                         0.009976
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  25

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2     0.00000000   -40.63480107
   ht   3     0.00000000    -0.00000000   -40.38147936
   ht   4    -0.00000000     0.00000000    -0.00000000   -40.23831159
   ht   5     0.00000000     0.00000000    -0.00000000    -0.00000000   -40.16513106
   ht   6    -0.03116644     0.01553066    -0.00013899    -0.00129487    -0.00113754    -0.00003396

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958736      -1.231103E-02   7.551560E-03   2.549333E-02  -1.464837E-02   0.125642    
 ref    2  -1.110508E-06   4.490160E-04   0.993537      -2.685310E-02   0.104955      -2.249374E-02
 ref    3   6.912409E-03  -0.945585       7.756894E-03  -4.173962E-02  -7.381565E-02   0.160390    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919222       0.894283       0.987233       3.113194E-03   1.667895E-02   4.201688E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873586    -0.01231103     0.00755156     0.02549333    -0.01464837     0.12564186
 ref:   2    -0.00000111     0.00044902     0.99353705    -0.02685310     0.10495535    -0.02249374
 ref:   3     0.00691241    -0.94558494     0.00775689    -0.04173962    -0.07381565     0.16039025

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 25  1    -94.6661949579  2.3881E-11  0.0000E+00  6.3617E-05  1.0000E-04   
 mr-sdci # 25  2    -94.3093370160  3.0847E-08 -1.0906E-07  4.6899E-04  1.0000E-04   
 mr-sdci # 25  3    -94.0573396186  1.3243E-03  0.0000E+00  9.8297E-01  1.0000E-04   
 mr-sdci # 25  4    -93.9128524317  4.9265E-06  0.0000E+00  3.0809E-01  1.0000E-04   
 mr-sdci # 25  5    -93.8517347877  1.2068E-02  0.0000E+00  4.8747E-01  1.0000E-04   
 mr-sdci # 25  6    -91.4185319734 -1.5193E+00  0.0000E+00  1.6148E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000004
time for cinew                         0.004725
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  26

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2     0.00000000   -40.63480107
   ht   3     0.00000000    -0.00000000   -40.38147936
   ht   4    -0.00000000     0.00000000    -0.00000000   -40.23831159
   ht   5     0.00000000     0.00000000    -0.00000000    -0.00000000   -40.16513106
   ht   6    -0.03116644     0.01553066    -0.00013899    -0.00129487    -0.00113754    -0.00003396
   ht   7    -0.07521865     0.02529245    -0.00063462    -0.00386919    -0.00397235    -0.00007186    -0.00016304

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735      -1.234094E-02   2.447111E-03   2.493639E-02   4.731199E-03  -0.224861       0.146295    
 ref    2  -1.351956E-06   4.736563E-04   0.990958      -2.892481E-02   0.127047       1.801089E-02   1.259729E-02
 ref    3   6.912258E-03  -0.945577       1.004923E-02  -4.096950E-02  -7.612071E-02  -0.124523      -0.106880    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.894268       0.982104       3.136968E-03   2.195778E-02   6.639276E-02   3.298426E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873514    -0.01234094     0.00244711     0.02493639     0.00473120    -0.22486080     0.14629505
 ref:   2    -0.00000135     0.00047366     0.99095769    -0.02892481     0.12704735     0.01801089     0.01259729
 ref:   3     0.00691226    -0.94557664     0.01004923    -0.04096950    -0.07612071    -0.12452303    -0.10687999

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 26  1    -94.6661949580  3.1065E-11  0.0000E+00  6.2971E-05  1.0000E-04   
 mr-sdci # 26  2    -94.3093370640  4.7992E-08  2.6317E-08  2.3752E-04  1.0000E-04   
 mr-sdci # 26  3    -94.0586400385  1.3004E-03  0.0000E+00  9.7760E-01  1.0000E-04   
 mr-sdci # 26  4    -93.9128681371  1.5705E-05  0.0000E+00  3.0884E-01  1.0000E-04   
 mr-sdci # 26  5    -93.8675997848  1.5865E-02  0.0000E+00  4.4528E-01  1.0000E-04   
 mr-sdci # 26  6    -91.4990460797  8.0514E-02  0.0000E+00  1.4776E+00  1.0000E-04   
 mr-sdci # 26  7    -91.1568458481 -2.6572E-01  0.0000E+00  1.4961E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000276
time for cinew                         0.007552
time for eigenvalue solver             0.000000
time for vector access                 0.000002

          starting ci iteration  27

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2     0.00000000   -40.63480107
   ht   3     0.00000000    -0.00000000   -40.38147936
   ht   4    -0.00000000     0.00000000    -0.00000000   -40.23831159
   ht   5     0.00000000     0.00000000    -0.00000000    -0.00000000   -40.16513106
   ht   6    -0.03116644     0.01553066    -0.00013899    -0.00129487    -0.00113754    -0.00003396
   ht   7    -0.07521865     0.02529245    -0.00063462    -0.00386919    -0.00397235    -0.00007186    -0.00016304
   ht   8    -0.00918166    -0.00197687     0.00011847    -0.00111528    -0.00158870    -0.00000665    -0.00001703    -0.00000308

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958735      -1.233827E-02  -4.937566E-03   2.492967E-02  -2.583439E-03   6.874959E-02   0.243829      -8.850974E-02
 ref    2  -2.706680E-06   5.577292E-04  -0.968589      -2.764062E-02  -0.215553      -0.118541       1.420959E-02  -2.621273E-03
 ref    3   6.912918E-03  -0.945590       1.049872E-03  -4.103076E-02   6.279041E-02  -0.102345       0.142805       0.162456    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.894294       0.938191       3.069015E-03   5.041246E-02   2.925295E-02   8.004800E-02   3.423282E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873523    -0.01233827    -0.00493757     0.02492967    -0.00258344     0.06874959     0.24382936    -0.08850974
 ref:   2    -0.00000271     0.00055773    -0.96858928    -0.02764062    -0.21555312    -0.11854078     0.01420959    -0.00262127
 ref:   3     0.00691292    -0.94559037     0.00104987    -0.04103076     0.06279041    -0.10234515     0.14280522     0.16245608

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 27  1    -94.6661949580  5.4825E-11  0.0000E+00  6.2103E-05  1.0000E-04   
 mr-sdci # 27  2    -94.3093370892  2.5224E-08  1.1112E-08  1.6226E-04  1.0000E-04   
 mr-sdci # 27  3    -94.0815294400  2.2889E-02  0.0000E+00  9.1180E-01  1.0000E-04   
 mr-sdci # 27  4    -93.9128688815  7.4438E-07  0.0000E+00  3.0864E-01  1.0000E-04   
 mr-sdci # 27  5    -93.8794013055  1.1802E-02  0.0000E+00  4.8305E-01  1.0000E-04   
 mr-sdci # 27  6    -92.5089254621  1.0099E+00  0.0000E+00  1.5310E+00  1.0000E-04   
 mr-sdci # 27  7    -91.3257409384  1.6890E-01  0.0000E+00  1.3387E+00  1.0000E-04   
 mr-sdci # 27  8    -91.1308241317  1.4652E-01  0.0000E+00  1.5501E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000023
time for cinew                         0.008229
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  28

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2     0.00000000   -40.63480118
   ht   3     0.00000000     0.00000000   -40.40699353
   ht   4    -0.00000000    -0.00000000    -0.00000000   -40.23833297
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -40.20486539
   ht   6    -0.00028087     0.00250047    -0.00003688    -0.00046523    -0.00004543    -0.00000039

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6

   eig(s)   5.912637E-09    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1   6.851843E-06   0.625489       0.503484       0.446430      -0.379674       0.108698    
   x:   2  -6.153542E-05  -3.713039E-02   8.396806E-03   0.160250      -0.141087      -0.976199    
   x:   3   1.712418E-06   9.976675E-02  -0.289060      -0.505944      -0.806085       2.716583E-02
   x:   4   1.165186E-05  -0.622494      -0.132233       0.623549      -0.414770       0.184845    
   x:   5   1.116986E-06   0.458183      -0.803363       0.360871       0.118886       1.771989E-02
   x:   6   -1.00000        0.00000        0.00000        0.00000        0.00000       6.303572E-05
 bummer (warning):overlap matrix: # small eigenvalues=1

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735       1.234068E-02   2.170656E-03  -2.562452E-02  -2.609236E-03  -3.642413E-02
 ref    2   5.902343E-06  -6.980485E-04  -0.952675      -8.053982E-03  -0.220066      -0.171401    
 ref    3  -6.912177E-03   0.945578      -2.736694E-02   4.349334E-02   6.301313E-02   0.154363    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.894270       0.908344       2.613154E-03   5.240661E-02   5.453289E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873507     0.01234068     0.00217066    -0.02562452    -0.00260924    -0.03642413
 ref:   2     0.00000590    -0.00069805    -0.95267519    -0.00805398    -0.22006623    -0.17140103
 ref:   3    -0.00691218     0.94557762    -0.02736694     0.04349334     0.06301313     0.15436273

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 28  1    -94.6661949581  5.7540E-11  0.0000E+00  6.1312E-05  1.0000E-04   
 mr-sdci # 28  2    -94.3093371006  1.1403E-08  7.6899E-09  1.3996E-04  1.0000E-04   
 mr-sdci # 28  3    -94.1536914344  7.2162E-02  0.0000E+00  7.4444E-01  1.0000E-04   
 mr-sdci # 28  4    -93.9135486965  6.7981E-04  0.0000E+00  3.0576E-01  1.0000E-04   
 mr-sdci # 28  5    -93.8794163007  1.4995E-05  0.0000E+00  4.8602E-01  1.0000E-04   
 mr-sdci # 28  6    -91.6239945829 -8.8493E-01  0.0000E+00  1.6196E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000251
time for cinew                         0.004449
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  29

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2     0.00000000   -40.63480118
   ht   3     0.00000000     0.00000000   -40.40699353
   ht   4    -0.00000000    -0.00000000    -0.00000000   -40.23833297
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -40.20486539
   ht   6    -0.00028087     0.00250047    -0.00003688    -0.00046523    -0.00004543    -0.00000039
   ht   7    -0.00345183     0.00112449     0.00081838    -0.00046587     0.00031226    -0.00000011    -0.00000052

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7

   eig(s)   4.526334E-09   6.039294E-09    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1   7.862391E-05   3.092310E-05  -3.591260E-02  -1.155585E-02   0.282790      -0.493581      -0.821574    
   x:   2  -8.685142E-06  -6.691022E-05  -0.174106       7.210691E-02  -5.005243E-03  -0.839672       0.509327    
   x:   3  -2.049779E-05  -4.406789E-06  -0.457989      -0.295420       0.796481       0.180155       0.190095    
   x:   4   7.676401E-06   1.449279E-05  -0.849498       0.345938      -0.348615       0.113393      -0.155851    
   x:   5  -8.085535E-06  -1.277061E-06   0.192346       0.887536       0.405082       7.755920E-02   7.194398E-02
   x:   6   0.289335      -0.957228      -4.947026E-13   1.297196E-13  -7.297717E-13   5.000399E-05  -3.838102E-05
   x:   7  -0.957228      -0.289335        0.00000        0.00000        0.00000      -2.141198E-05  -8.963228E-05
 bummer (warning):overlap matrix: # small eigenvalues=2

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958735       1.232818E-02   1.942822E-02   2.450916E-02   1.105408E-02  -4.688098E-02  -0.108694    
 ref    2   5.972515E-06  -1.072156E-03   0.930550       1.463986E-02   9.517309E-02   0.333375       6.829787E-02
 ref    3  -6.912186E-03   0.945585      -2.714099E-03  -3.883141E-02  -7.519759E-02  -2.256794E-02   0.208549    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.894284       0.866308       2.322903E-03   1.483479E-02   0.113846       5.997162E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873507     0.01232818     0.01942822     0.02450916     0.01105408    -0.04688098    -0.10869404
 ref:   2     0.00000597    -0.00107216     0.93055015     0.01463986     0.09517309     0.33337535     0.06829787
 ref:   3    -0.00691219     0.94558479    -0.00271410    -0.03883141    -0.07519759    -0.02256794     0.20854885

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 29  1    -94.6661949581  0.0000E+00  0.0000E+00  6.1301E-05  1.0000E-04   
 mr-sdci # 29  2    -94.3093371145  1.3894E-08 -1.0448E-09  2.1457E-04  1.0000E-04   
 mr-sdci # 29  3    -94.2550665706  1.0138E-01  0.0000E+00  4.1810E-01  1.0000E-04   
 mr-sdci # 29  4    -93.9136943729  1.4568E-04  0.0000E+00  3.0675E-01  1.0000E-04   
 mr-sdci # 29  5    -93.8925287939  1.3112E-02  0.0000E+00  3.7854E-01  1.0000E-04   
 mr-sdci # 29  6    -92.0056627210  3.8167E-01  0.0000E+00  1.6341E+00  1.0000E-04   
 mr-sdci # 29  7    -91.1338437813 -1.9190E-01  0.0000E+00  1.4106E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000002
time for cinew                         0.003284
time for eigenvalue solver             0.000120
time for vector access                 0.000001

          starting ci iteration  30

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2     0.00000000   -40.63480118
   ht   3     0.00000000     0.00000000   -40.40699353
   ht   4    -0.00000000    -0.00000000    -0.00000000   -40.23833297
   ht   5    -0.00000000    -0.00000000     0.00000000    -0.00000000   -40.20486539
   ht   6    -0.00028087     0.00250047    -0.00003688    -0.00046523    -0.00004543    -0.00000039
   ht   7    -0.00345183     0.00112449     0.00081838    -0.00046587     0.00031226    -0.00000011    -0.00000052
   ht   8    -0.02369285    -0.01017424    -0.00072521    -0.00104922     0.00020671     0.00000047    -0.00000178    -0.00001742

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7         v:   8

   eig(s)   4.393343E-09   6.039052E-09   3.018048E-08    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1  -3.680176E-05  -3.296674E-05  -5.820407E-04  -0.109787       3.010899E-02  -3.740407E-02  -0.374779      -0.919337    
   x:   2   2.641249E-05   6.618228E-05  -2.493312E-04   0.166740      -1.376427E-02   0.151949       0.892483      -0.390377    
   x:   3   2.164972E-05   4.431109E-06  -1.548983E-05  -0.167856      -0.301750      -0.922012       0.173584      -2.308793E-02
   x:   4  -5.732736E-06  -1.460016E-05  -2.652121E-05   0.961059      -0.136449      -0.161747      -0.172615      -4.228837E-02
   x:   5   7.658794E-06   1.323934E-06   6.150584E-06   9.129016E-02   0.942991      -0.315029       5.556178E-02   1.014864E-02
   x:   6  -0.291912       0.956146      -2.392662E-02   3.405236E-12  -6.870942E-12  -1.342192E-05  -5.913923E-05   1.720190E-05
   x:   7   0.953745       0.292876       6.778471E-02  -2.487104E-12   9.670424E-12   1.259970E-05  -6.232724E-05  -6.670067E-05
   x:   8  -7.181962E-02   3.032685E-03   0.997413       3.558446E-13  -1.210995E-12  -1.698668E-06   4.978809E-06  -6.306653E-04
 bummer (warning):overlap matrix: # small eigenvalues=2

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735      -1.232866E-02  -2.269957E-02   2.633374E-02  -1.487991E-02   6.689336E-02  -0.207963      -2.313903E-02
 ref    2   5.641010E-06   1.070085E-03  -0.932269       1.656758E-02  -9.445722E-02   0.283070       0.183239       1.921210E-02
 ref    3  -6.912466E-03  -0.945585       2.001115E-05  -3.904081E-02   7.027089E-02   6.326706E-02  -6.777459E-02   0.262335    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919220       0.894284       0.869641       2.492136E-03   1.408158E-02   8.860608E-02   8.141843E-02   6.972430E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873475    -0.01232866    -0.02269957     0.02633374    -0.01487991     0.06689336    -0.20796271    -0.02313903
 ref:   2     0.00000564     0.00107008    -0.93226916     0.01656758    -0.09445722     0.28307002     0.18323905     0.01921210
 ref:   3    -0.00691247    -0.94558515     0.00002001    -0.03904081     0.07027089     0.06326706    -0.06777459     0.26233524

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 30  1    -94.6661949581  8.8463E-12  0.0000E+00  6.1187E-05  1.0000E-04   
 mr-sdci # 30  2    -94.3093371145  1.4488E-11  1.0016E-09  2.1257E-04  1.0000E-04   
 mr-sdci # 30  3    -94.2558569113  7.9034E-04  0.0000E+00  4.2441E-01  1.0000E-04   
 mr-sdci # 30  4    -93.9138559461  1.6157E-04  0.0000E+00  3.0542E-01  1.0000E-04   
 mr-sdci # 30  5    -93.8938410929  1.3123E-03  0.0000E+00  3.8309E-01  1.0000E-04   
 mr-sdci # 30  6    -92.2507588100  2.4510E-01  0.0000E+00  1.5949E+00  1.0000E-04   
 mr-sdci # 30  7    -91.4323661840  2.9852E-01  0.0000E+00  1.4209E+00  1.0000E-04   
 mr-sdci # 30  8    -91.0703506720 -6.0473E-02  0.0000E+00  1.4117E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000168
time for cinew                         0.010359
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  31

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.63480120
   ht   3    -0.00000000    -0.00000000   -40.58132100
   ht   4     0.00000000     0.00000000     0.00000000   -40.23932003
   ht   5    -0.00000000     0.00000000    -0.00000000     0.00000000   -40.21930518
   ht   6    -0.02229225     0.00978901     0.00104616     0.00088810    -0.00005195    -0.00001558

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735       1.232823E-02   2.352502E-02  -2.638398E-02  -1.525783E-02   0.190507    
 ref    2  -5.640340E-06  -1.070496E-03   0.932194      -1.656076E-02  -9.436355E-02  -2.073736E-02
 ref    3   6.912477E-03   0.945585       7.123265E-04   3.900080E-02   6.992830E-02   0.170316    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919220       0.894284       0.869540       2.491435E-03   1.402725E-02   6.573047E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873474     0.01232823     0.02352502    -0.02638398    -0.01525783     0.19050666
 ref:   2    -0.00000564    -0.00107050     0.93219436    -0.01656076    -0.09436355    -0.02073736
 ref:   3     0.00691248     0.94558478     0.00071233     0.03900080     0.06992830     0.17031631

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 31  1    -94.6661949581  4.2633E-14  0.0000E+00  6.1189E-05  1.0000E-04   
 mr-sdci # 31  2    -94.3093371145  1.3358E-11 -1.0373E-09  2.1333E-04  1.0000E-04   
 mr-sdci # 31  3    -94.2559069064  4.9995E-05  0.0000E+00  4.2571E-01  1.0000E-04   
 mr-sdci # 31  4    -93.9138561034  1.5732E-07  0.0000E+00  3.0541E-01  1.0000E-04   
 mr-sdci # 31  5    -93.8938503200  9.2271E-06  0.0000E+00  3.8365E-01  1.0000E-04   
 mr-sdci # 31  6    -91.5837404687 -6.6702E-01  0.0000E+00  1.4863E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001194
time for cinew                         0.006608
time for eigenvalue solver             0.000092
time for vector access                 0.000000

          starting ci iteration  32

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.63480120
   ht   3    -0.00000000    -0.00000000   -40.58132100
   ht   4     0.00000000     0.00000000     0.00000000   -40.23932003
   ht   5    -0.00000000     0.00000000    -0.00000000     0.00000000   -40.21930518
   ht   6    -0.02229225     0.00978901     0.00104616     0.00088810    -0.00005195    -0.00001558
   ht   7     0.02359599    -0.01013769    -0.00113650    -0.00095412     0.00009591     0.00001641    -0.00001729

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7

   eig(s)   1.616837E-11   5.728338E-08    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1   1.003217E-05   7.918285E-04   9.173672E-02   5.225422E-04   0.115755       0.376704      -0.914483    
   x:   2  -6.109418E-07  -3.468069E-04   0.180548       4.421876E-02   0.132754       0.887378       0.400479    
   x:   3  -1.009899E-06  -3.801894E-05   0.547748      -0.682365       0.452977      -0.165029       4.391475E-02
   x:   4  -6.487168E-07  -3.238827E-05  -0.391109       0.285883       0.870187      -8.173431E-02   3.740823E-02
   x:   5   7.178047E-07   2.624489E-06  -0.711331      -0.671338      -8.102212E-02   0.191676      -3.039685E-03
   x:   6   0.718145      -0.695893       2.473799E-11  -8.911782E-15  -3.351008E-11  -2.608634E-06  -5.957536E-04
   x:   7   0.695894       0.718145       2.346320E-11    0.00000      -3.175738E-11  -2.472955E-06   6.284396E-04
 bummer (warning):overlap matrix: # small eigenvalues=1

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735       1.233495E-02  -1.255774E-03  -2.523897E-02  -2.598766E-03  -8.273088E-03  -0.270188    
 ref    2  -3.548394E-06  -4.595029E-03  -0.936168      -1.374976E-02  -4.664936E-02  -6.594260E-02  -1.678845E-02
 ref    3   6.912720E-03   0.945523      -1.650655E-02   4.327459E-02   5.969196E-02   0.193016      -5.214501E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.894187       0.876684       2.698751E-03   5.746046E-03   4.167193E-02   7.600236E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873510     0.01233495    -0.00125577    -0.02523897    -0.00259877    -0.00827309    -0.27018772
 ref:   2    -0.00000355    -0.00459503    -0.93616772    -0.01374976    -0.04664936    -0.06594260    -0.01678845
 ref:   3     0.00691272     0.94552297    -0.01650655     0.04327459     0.05969196     0.19301569    -0.05214501

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 32  1    -94.6661949581  1.5213E-11  0.0000E+00  6.0756E-05  1.0000E-04   
 mr-sdci # 32  2    -94.3093372173  1.0277E-07  6.1078E-07  1.1361E-03  1.0000E-04   
 mr-sdci # 32  3    -94.3029473039  4.7040E-02  0.0000E+00  2.7148E-01  1.0000E-04   
 mr-sdci # 32  4    -93.9140934316  2.3733E-04  0.0000E+00  3.0403E-01  1.0000E-04   
 mr-sdci # 32  5    -93.9014651134  7.6148E-03  0.0000E+00  3.5735E-01  1.0000E-04   
 mr-sdci # 32  6    -91.8573575560  2.7362E-01  0.0000E+00  1.6420E+00  1.0000E-04   
 mr-sdci # 32  7    -91.3371652685 -9.5201E-02  0.0000E+00  1.2970E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.005245
time for eigenvalue solver             0.000137
time for vector access                 0.000002

          starting ci iteration  33

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.63480120
   ht   3    -0.00000000    -0.00000000   -40.58132100
   ht   4     0.00000000     0.00000000     0.00000000   -40.23932003
   ht   5    -0.00000000     0.00000000    -0.00000000     0.00000000   -40.21930518
   ht   6    -0.02229225     0.00978901     0.00104616     0.00088810    -0.00005195    -0.00001558
   ht   7     0.02359599    -0.01013769    -0.00113650    -0.00095412     0.00009591     0.00001641    -0.00001729
   ht   8     0.02929349     0.00976322    -0.00806861    -0.00356188     0.00185140     0.00001461    -0.00001554    -0.00004110

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7         v:   8

   eig(s)   1.509083E-11   5.521600E-08   4.128219E-07    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1   9.473018E-06  -7.351920E-04   7.727732E-04  -0.133569       2.055130E-02  -0.131914      -3.786893E-02   0.981275    
   x:   2  -1.276546E-06   3.640707E-04   2.132015E-04  -0.226325      -4.006454E-03  -0.202656      -0.948025      -9.455207E-02
   x:   3  -7.132995E-07   2.283624E-05  -2.005436E-04  -0.455896      -0.328962      -0.761437       0.287602      -0.146428    
   x:   4  -5.278188E-07   2.554248E-05  -9.100231E-05  -0.240824       0.939247      -0.207122       0.105402      -7.622759E-02
   x:   5   6.430141E-07   9.335972E-07   4.675347E-05   0.815536       9.571452E-02  -0.564662      -7.739716E-02   3.010956E-02
   x:   6   0.717602       0.694526      -5.177578E-02  -2.259721E-10  -1.733291E-10   5.070770E-07   1.979586E-04   5.619086E-04
   x:   7   0.696452      -0.715438       5.570420E-02  -2.151459E-10  -1.656014E-10   4.848395E-07  -2.039792E-04  -5.944195E-04
   x:   8  -1.645694E-03   7.603326E-02   0.997104       1.468917E-12   1.468917E-12  -4.565051E-09   3.248243E-04  -7.157222E-04
 bummer (warning):overlap matrix: # small eigenvalues=1

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735       1.975735E-03   1.234184E-02   2.269506E-02   1.126927E-02   5.840010E-03   0.253394       9.464899E-02
 ref    2   3.428828E-07   0.939119       1.341074E-04  -3.686081E-04  -2.454583E-03   9.283582E-02  -7.214146E-03   7.889310E-02
 ref    3  -6.912330E-03  -3.098828E-03   0.945581      -7.323805E-02   3.154263E-02  -3.681548E-03   0.122315      -0.211861    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.881958       0.894276       5.879014E-03   1.127959E-03   8.666148E-03   7.922148E-02   6.006761E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873507     0.00197574     0.01234184     0.02269506     0.01126927     0.00584001     0.25339394     0.09464899
 ref:   2     0.00000034     0.93911900     0.00013411    -0.00036861    -0.00245458     0.09283582    -0.00721415     0.07889310
 ref:   3    -0.00691233    -0.00309883     0.94558108    -0.07323805     0.03154263    -0.00368155     0.12231495    -0.21186094

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 33  1    -94.6661949581  3.5122E-11  0.0000E+00  5.9835E-05  1.0000E-04   
 mr-sdci # 33  2    -94.3399392605  3.0602E-02  7.3270E-03  1.1317E-01  1.0000E-04   
 mr-sdci # 33  3    -94.3093370825  6.3898E-03  0.0000E+00  1.1341E-04  1.0000E-04   
 mr-sdci # 33  4    -93.9174341786  3.3407E-03  0.0000E+00  2.9447E-01  1.0000E-04   
 mr-sdci # 33  5    -93.9097126474  8.2475E-03  0.0000E+00  3.4424E-01  1.0000E-04   
 mr-sdci # 33  6    -93.0382316317  1.1809E+00  0.0000E+00  1.1095E+00  1.0000E-04   
 mr-sdci # 33  7    -91.3484140833  1.1249E-02  0.0000E+00  1.3386E+00  1.0000E-04   
 mr-sdci # 33  8    -91.2069185243  1.3657E-01  0.0000E+00  1.4964E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000343
time for cinew                         0.008781
time for eigenvalue solver             0.000042
time for vector access                 0.000001

          starting ci iteration  34

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2     0.00000000   -40.66540335
   ht   3    -0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24289827
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -40.23517673
   ht   6    -1.20096511    -2.54742897    -1.81600182    -0.42091195     0.10609532    -0.48638723

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735       2.300947E-03  -1.234191E-02  -1.954193E-02  -1.619109E-02   7.399981E-03
 ref    2  -5.653306E-07   0.948505      -3.927976E-05  -1.349076E-03   3.589225E-03   0.222707    
 ref    3   6.912289E-03   7.604583E-04  -0.945585       7.528346E-02  -9.975235E-03   0.112787    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.899667       0.894283       6.051306E-03   3.745393E-04   6.237394E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873508     0.00230095    -0.01234191    -0.01954193    -0.01619109     0.00739998
 ref:   2    -0.00000057     0.94850481    -0.00003928    -0.00134908     0.00358922     0.22270660
 ref:   3     0.00691229     0.00076046    -0.94558464     0.07528346    -0.00997523     0.11278720

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 34  1    -94.6661949581  3.9790E-13  0.0000E+00  5.9831E-05  1.0000E-04   
 mr-sdci # 34  2    -94.3438287863  3.8895E-03 -1.6345E-03  4.6085E-02  1.0000E-04   
 mr-sdci # 34  3    -94.3093370858  3.2631E-09  0.0000E+00  5.3832E-05  1.0000E-04   
 mr-sdci # 34  4    -93.9189730721  1.5389E-03  0.0000E+00  2.9452E-01  1.0000E-04   
 mr-sdci # 34  5    -93.9117881371  2.0755E-03  0.0000E+00  3.3001E-01  1.0000E-04   
 mr-sdci # 34  6    -91.7629700996 -1.2753E+00  0.0000E+00  1.6559E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000154
time for cinew                         0.003990
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  35

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2     0.00000000   -40.66540335
   ht   3    -0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24289827
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -40.23517673
   ht   6    -1.20096511    -2.54742897    -1.81600182    -0.42091195     0.10609532    -0.48638723
   ht   7   -17.37408823     0.35624328    -5.53324697     1.19998657    -0.06220082    -0.76626451    -8.75545511

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735       1.177152E-03  -1.234109E-02  -1.869833E-02   1.718691E-02  -1.350509E-02   0.262201    
 ref    2  -5.888518E-07   0.948725      -3.833508E-05  -1.434411E-03  -3.670202E-03  -0.221493      -5.545462E-02
 ref    3   6.912307E-03   4.457560E-04  -0.945584       7.553584E-02   9.941330E-03  -0.114418       6.878506E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.900081       0.894282       6.057348E-03   4.076905E-04   6.233310E-02   7.655604E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873501     0.00117715    -0.01234109    -0.01869833     0.01718691    -0.01350509     0.26220115
 ref:   2    -0.00000059     0.94872527    -0.00003834    -0.00143441    -0.00367020    -0.22149344    -0.05545462
 ref:   3     0.00691231     0.00044576    -0.94558441     0.07553584     0.00994133    -0.11441752     0.06878506

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 35  1    -94.6661949581  2.3448E-13  0.0000E+00  5.9844E-05  1.0000E-04   
 mr-sdci # 35  2    -94.3438847069  5.5921E-05  6.3342E-04  4.3835E-02  1.0000E-04   
 mr-sdci # 35  3    -94.3093370858  2.8813E-11  0.0000E+00  5.2962E-05  1.0000E-04   
 mr-sdci # 35  4    -93.9189962469  2.3175E-05  0.0000E+00  2.9440E-01  1.0000E-04   
 mr-sdci # 35  5    -93.9118204970  3.2360E-05  0.0000E+00  3.2924E-01  1.0000E-04   
 mr-sdci # 35  6    -91.7632247320  2.5463E-04  0.0000E+00  1.6551E+00  1.0000E-04   
 mr-sdci # 35  7    -91.2929846453 -5.5429E-02  0.0000E+00  1.3038E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.001575
time for cinew                         0.006694
time for eigenvalue solver             0.000141
time for vector access                 0.000000

          starting ci iteration  36

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2     0.00000000   -40.66540335
   ht   3    -0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24289827
   ht   5     0.00000000     0.00000000     0.00000000     0.00000000   -40.23517673
   ht   6    -1.20096511    -2.54742897    -1.81600182    -0.42091195     0.10609532    -0.48638723
   ht   7   -17.37408823     0.35624328    -5.53324697     1.19998657    -0.06220082    -0.76626451    -8.75545511
   ht   8    -6.48925971     0.45263503    -1.88513178     0.49907300    -0.03026990    -0.25970621    -3.26131460    -1.22946312

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958735       3.477708E-04  -1.234019E-02  -1.610616E-02  -1.993765E-02  -1.942023E-02   0.259925      -4.596943E-02
 ref    2  -6.544572E-07   0.947805      -5.294669E-06   4.756096E-04   6.584509E-04   3.453082E-03  -7.129795E-02  -0.241484    
 ref    3   6.912338E-03  -6.632792E-04  -0.945583       7.682816E-02  -2.936828E-03  -2.571521E-03   5.982648E-02  -0.141854    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.898334       0.894280       6.162200E-03   4.065684E-04   3.956819E-04   7.622336E-02   8.055010E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873499     0.00034777    -0.01234019    -0.01610616    -0.01993765    -0.01942023     0.25992452    -0.04596943
 ref:   2    -0.00000065     0.94780462    -0.00000529     0.00047561     0.00065845     0.00345308    -0.07129795    -0.24148382
 ref:   3     0.00691234    -0.00066328    -0.94558321     0.07682816    -0.00293683    -0.00257152     0.05982648    -0.14185370

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 36  1    -94.6661949581  9.2371E-13  0.0000E+00  5.9760E-05  1.0000E-04   
 mr-sdci # 36  2    -94.3450226653  1.1380E-03  3.3171E-04  2.6897E-02  1.0000E-04   
 mr-sdci # 36  3    -94.3093370871  1.2411E-09  0.0000E+00  3.0860E-05  1.0000E-04   
 mr-sdci # 36  4    -93.9195022560  5.0601E-04  0.0000E+00  2.9708E-01  1.0000E-04   
 mr-sdci # 36  5    -93.9128288470  1.0083E-03  0.0000E+00  3.2302E-01  1.0000E-04   
 mr-sdci # 36  6    -93.3310029512  1.5678E+00  0.0000E+00  1.0677E+00  1.0000E-04   
 mr-sdci # 36  7    -91.2934458035  4.6116E-04  0.0000E+00  1.3094E+00  1.0000E-04   
 mr-sdci # 36  8    -91.1587375736 -4.8181E-02  0.0000E+00  1.4577E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000359
time for cinew                         0.010945
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  37

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67048675
   ht   3    -0.00000000     0.00000000   -40.63480117
   ht   4    -0.00000000    -0.00000000     0.00000000   -40.24496634
   ht   5    -0.00000000     0.00000000    -0.00000000     0.00000000   -40.23829293
   ht   6    -0.00134388    -0.08503539     0.13489548    -0.02423197    -0.01934109    -0.00925697

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735       1.890454E-04  -1.234006E-02  -1.591329E-02   1.998017E-02  -1.621781E-02
 ref    2   5.974208E-07   0.948085      -8.535367E-07   6.396194E-04  -6.469694E-04   2.252505E-02
 ref    3  -6.912388E-03   1.959410E-04  -0.945584       7.599960E-02   2.726481E-03   8.719694E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.898865       0.894281       6.029582E-03   4.070596E-04   8.373701E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873500     0.00018905    -0.01234006    -0.01591329     0.01998017    -0.01621781
 ref:   2     0.00000060     0.94808474    -0.00000085     0.00063962    -0.00064697     0.02252505
 ref:   3    -0.00691239     0.00019594    -0.94558388     0.07599960     0.00272648     0.08719694

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 37  1    -94.6661949581  7.8160E-13  0.0000E+00  5.9768E-05  1.0000E-04   
 mr-sdci # 37  2    -94.3452499070  2.2724E-04  3.1884E-05  9.3200E-03  1.0000E-04   
 mr-sdci # 37  3    -94.3093370872  1.3623E-10  0.0000E+00  2.4364E-05  1.0000E-04   
 mr-sdci # 37  4    -93.9196629630  1.6071E-04  0.0000E+00  2.9757E-01  1.0000E-04   
 mr-sdci # 37  5    -93.9128298564  1.0094E-06  0.0000E+00  3.2293E-01  1.0000E-04   
 mr-sdci # 37  6    -92.1781769922 -1.1528E+00  0.0000E+00  1.7151E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000123
time for cinew                         0.006573
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  38

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67048675
   ht   3    -0.00000000     0.00000000   -40.63480117
   ht   4    -0.00000000    -0.00000000     0.00000000   -40.24496634
   ht   5    -0.00000000     0.00000000    -0.00000000     0.00000000   -40.23829293
   ht   6    -0.00134388    -0.08503539     0.13489548    -0.02423197    -0.01934109    -0.00925697
   ht   7     1.50402701     0.08922309     0.40516127    -0.11597241    -0.02876117    -0.00113129    -0.06549592

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735       4.187932E-04   1.234015E-02   1.662752E-02  -2.029075E-02  -6.501699E-02   0.230448    
 ref    2  -6.117558E-07   0.948047       8.151464E-07  -7.688955E-04   7.079700E-04  -6.926553E-03  -4.732640E-02
 ref    3   6.912411E-03   2.459338E-04   0.945584      -7.583585E-02  -2.905714E-03  -9.698237E-02   1.224188E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919220       0.898793       0.894281       6.028142E-03   4.206588E-04   1.368077E-02   5.549584E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873489     0.00041879     0.01234015     0.01662752    -0.02029075    -0.06501699     0.23044781
 ref:   2    -0.00000061     0.94804656     0.00000082    -0.00076890     0.00070797    -0.00692655    -0.04732640
 ref:   3     0.00691241     0.00024593     0.94558390    -0.07583585    -0.00290571    -0.09698237     0.01224188

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 38  1    -94.6661949581  6.3949E-13  0.0000E+00  5.9804E-05  1.0000E-04   
 mr-sdci # 38  2    -94.3452524173  2.5103E-06 -7.8947E-05  9.4605E-03  1.0000E-04   
 mr-sdci # 38  3    -94.3093370872  3.4106E-13  0.0000E+00  2.4385E-05  1.0000E-04   
 mr-sdci # 38  4    -93.9196819549  1.8992E-05  0.0000E+00  2.9768E-01  1.0000E-04   
 mr-sdci # 38  5    -93.9128342859  4.4294E-06  0.0000E+00  3.2277E-01  1.0000E-04   
 mr-sdci # 38  6    -92.2574252963  7.9248E-02  0.0000E+00  1.6750E+00  1.0000E-04   
 mr-sdci # 38  7    -91.5560266252  2.6258E-01  0.0000E+00  1.5112E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000146
time for cinew                         0.007818
time for eigenvalue solver             0.000142
time for vector access                 0.000001

          starting ci iteration  39

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67048675
   ht   3    -0.00000000     0.00000000   -40.63480117
   ht   4    -0.00000000    -0.00000000     0.00000000   -40.24496634
   ht   5    -0.00000000     0.00000000    -0.00000000     0.00000000   -40.23829293
   ht   6    -0.00134388    -0.08503539     0.13489548    -0.02423197    -0.01934109    -0.00925697
   ht   7     1.50402701     0.08922309     0.40516127    -0.11597241    -0.02876117    -0.00113129    -0.06549592
   ht   8     4.01370456     0.07503658     1.25206096    -0.27247421    -0.07328477    -0.00441731    -0.17355960    -0.46612301

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735       1.148373E-04  -1.233999E-02   1.623590E-02  -1.974670E-02   4.451592E-03   2.623548E-02   0.263757    
 ref    2   6.924719E-07   0.947574       8.581029E-08  -2.518050E-03   1.364039E-03   6.843265E-02  -5.693519E-02  -4.152841E-03
 ref    3  -6.912276E-03  -2.433507E-04  -0.945584      -7.705439E-02  -3.633470E-03   3.415626E-02  -0.110715       9.379495E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.897896       0.894281       6.207324E-03   4.049949E-04   5.869494E-03   1.618775E-02   7.838233E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873498     0.00011484    -0.01233999     0.01623590    -0.01974670     0.00445159     0.02623548     0.26375669
 ref:   2     0.00000069     0.94757387     0.00000009    -0.00251805     0.00136404     0.06843265    -0.05693519    -0.00415284
 ref:   3    -0.00691228    -0.00024335    -0.94558364    -0.07705439    -0.00363347     0.03415626    -0.11071511     0.09379495

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 39  1    -94.6661949581  4.3059E-12  0.0000E+00  5.9828E-05  1.0000E-04   
 mr-sdci # 39  2    -94.3453011324  4.8715E-05  2.2301E-05  6.4845E-03  1.0000E-04   
 mr-sdci # 39  3    -94.3093370872  1.2655E-11  0.0000E+00  2.4713E-05  1.0000E-04   
 mr-sdci # 39  4    -93.9199668689  2.8491E-04  0.0000E+00  2.9709E-01  1.0000E-04   
 mr-sdci # 39  5    -93.9128782856  4.4000E-05  0.0000E+00  3.2317E-01  1.0000E-04   
 mr-sdci # 39  6    -93.2987484298  1.0413E+00  0.0000E+00  1.1615E+00  1.0000E-04   
 mr-sdci # 39  7    -91.9120362416  3.5601E-01  0.0000E+00  1.7789E+00  1.0000E-04   
 mr-sdci # 39  8    -91.2859774264  1.2724E-01  0.0000E+00  1.2964E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000060
time for cinew                         0.010227
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  40

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67076522
   ht   3     0.00000000    -0.00000000   -40.63480117
   ht   4    -0.00000000     0.00000000     0.00000000   -40.24543096
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -40.23834237
   ht   6    -0.05097839    -0.08372075     0.12136537    -0.02712697     0.00574147    -0.00122431

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.958735       1.108472E-04  -1.233999E-02   1.623991E-02   1.971555E-02  -2.812449E-03
 ref    2  -5.350628E-07   0.947898       1.136975E-07  -1.045734E-03  -1.288270E-03   0.135380    
 ref    3   6.912418E-03   8.372787E-05  -0.945584      -7.537421E-02   3.854484E-03   0.142443    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.898511       0.894281       5.946100E-03   4.052195E-04   3.862548E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.95873498     0.00011085    -0.01233999     0.01623991     0.01971555    -0.00281245
 ref:   2    -0.00000054     0.94789796     0.00000011    -0.00104573    -0.00128827     0.13537975
 ref:   3     0.00691242     0.00008373    -0.94558402    -0.07537421     0.00385448     0.14244259

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 40  1    -94.6661949581  2.7782E-12  0.0000E+00  5.9698E-05  1.0000E-04   
 mr-sdci # 40  2    -94.3453144505  1.3318E-05 -5.4537E-06  2.6638E-03  1.0000E-04   
 mr-sdci # 40  3    -94.3093370872  1.7494E-11  0.0000E+00  2.4141E-05  1.0000E-04   
 mr-sdci # 40  4    -93.9202354095  2.6854E-04  0.0000E+00  2.9558E-01  1.0000E-04   
 mr-sdci # 40  5    -93.9128789303  6.4462E-07  0.0000E+00  3.2321E-01  1.0000E-04   
 mr-sdci # 40  6    -91.9867973928 -1.3120E+00  0.0000E+00  1.7295E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000334
time for cinew                         0.006322
time for eigenvalue solver             0.000099
time for vector access                 0.000002

          starting ci iteration  41

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67076522
   ht   3     0.00000000    -0.00000000   -40.63480117
   ht   4    -0.00000000     0.00000000     0.00000000   -40.24543096
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -40.23834237
   ht   6    -0.05097839    -0.08372075     0.12136537    -0.02712697     0.00574147    -0.00122431
   ht   7    -1.10198020     0.01970563     0.33764466     0.07833475    -0.02009242    -0.00243723    -0.03504888

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.958735       5.122280E-05   1.233994E-02   1.584076E-02   1.962300E-02  -1.672772E-03   0.265592    
 ref    2   5.314431E-07   0.947905      -1.109580E-07  -9.986922E-04  -1.277554E-03  -0.134799      -3.556175E-02
 ref    3  -6.912412E-03   7.033664E-05   0.945584      -7.546301E-02   3.848234E-03  -0.143403       5.569889E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.898524       0.894281       5.946592E-03   4.015031E-04   3.873809E-02   7.490589E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.95873501     0.00005122     0.01233994     0.01584076     0.01962300    -0.00167277     0.26559157
 ref:   2     0.00000053     0.94790532    -0.00000011    -0.00099869    -0.00127755    -0.13479949    -0.03556175
 ref:   3    -0.00691241     0.00007034     0.94558401    -0.07546301     0.00384823    -0.14340290     0.05569889

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 41  1    -94.6661949581  2.1316E-14  0.0000E+00  5.9687E-05  1.0000E-04   
 mr-sdci # 41  2    -94.3453146037  1.5323E-07  2.3170E-06  2.5457E-03  1.0000E-04   
 mr-sdci # 41  3    -94.3093370872  1.2079E-13  0.0000E+00  2.4145E-05  1.0000E-04   
 mr-sdci # 41  4    -93.9202414115  6.0020E-06  0.0000E+00  2.9545E-01  1.0000E-04   
 mr-sdci # 41  5    -93.9128792268  2.9654E-07  0.0000E+00  3.2326E-01  1.0000E-04   
 mr-sdci # 41  6    -91.9869903255  1.9293E-04  0.0000E+00  1.7291E+00  1.0000E-04   
 mr-sdci # 41  7    -91.3103651077 -6.0167E-01  0.0000E+00  1.3122E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000936
time for cinew                         0.006388
time for eigenvalue solver             0.000117
time for vector access                 0.000002

          starting ci iteration  42

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67076522
   ht   3     0.00000000    -0.00000000   -40.63480117
   ht   4    -0.00000000     0.00000000     0.00000000   -40.24543096
   ht   5    -0.00000000    -0.00000000     0.00000000     0.00000000   -40.23834237
   ht   6    -0.05097839    -0.08372075     0.12136537    -0.02712697     0.00574147    -0.00122431
   ht   7    -1.10198020     0.01970563     0.33764466     0.07833475    -0.02009242    -0.00243723    -0.03504888
   ht   8    -0.45375604     0.02218282     0.11946253     0.03749665    -0.00919736    -0.00092688    -0.01433623    -0.00592517

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958735       2.171621E-05  -1.233989E-02  -1.563772E-02  -1.975856E-02  -9.163875E-03   0.183753      -0.193106    
 ref    2  -5.473218E-07   0.947848       3.927637E-07   1.178787E-03   1.175147E-03   1.345106E-03  -0.135166      -9.178161E-02
 ref    3   6.912373E-03  -5.345086E-05  -0.945584       7.591932E-02  -3.807374E-03  -2.093572E-02  -8.782551E-02  -0.176172    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.898416       0.894281       6.009672E-03   4.062778E-04   5.240904E-04   5.974832E-02   7.675026E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873502     0.00002172    -0.01233989    -0.01563772    -0.01975856    -0.00916387     0.18375266    -0.19310558
 ref:   2    -0.00000055     0.94784826     0.00000039     0.00117879     0.00117515     0.00134511    -0.13516642    -0.09178161
 ref:   3     0.00691237    -0.00005345    -0.94558380     0.07591932    -0.00380737    -0.02093572    -0.08782551    -0.17617216

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 42  1    -94.6661949581  5.0449E-13  0.0000E+00  5.9721E-05  1.0000E-04   
 mr-sdci # 42  2    -94.3453189346  4.3309E-06  1.4869E-06  1.7867E-03  1.0000E-04   
 mr-sdci # 42  3    -94.3093370873  1.1823E-11  0.0000E+00  2.3934E-05  1.0000E-04   
 mr-sdci # 42  4    -93.9202913558  4.9944E-05  0.0000E+00  2.9549E-01  1.0000E-04   
 mr-sdci # 42  5    -93.9128961753  1.6948E-05  0.0000E+00  3.2314E-01  1.0000E-04   
 mr-sdci # 42  6    -93.5149758689  1.5280E+00  0.0000E+00  9.8794E-01  1.0000E-04   
 mr-sdci # 42  7    -91.3605331471  5.0168E-02  0.0000E+00  1.5023E+00  1.0000E-04   
 mr-sdci # 42  8    -91.2464983467 -3.9479E-02  0.0000E+00  1.4015E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000204
time for cinew                         0.009933
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  43

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67078302
   ht   3     0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24575544
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -40.23836026
   ht   6    -0.00107237     0.00001981     0.01349475    -0.00010250    -0.00065892    -0.00004006

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735       7.960927E-06  -1.233986E-02   1.550617E-02  -1.975584E-02  -2.127092E-02
 ref    2   5.434574E-07   0.947851       4.526294E-07  -1.173093E-03   1.174978E-03   3.742879E-03
 ref    3  -6.912424E-03   2.435667E-05  -0.945584      -7.516779E-02  -3.821531E-03   0.119727    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.898421       0.894281       5.892013E-03   4.062780E-04   1.480095E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873503     0.00000796    -0.01233986     0.01550617    -0.01975584    -0.02127092
 ref:   2     0.00000054     0.94785092     0.00000045    -0.00117309     0.00117498     0.00374288
 ref:   3    -0.00691242     0.00002436    -0.94558399    -0.07516779    -0.00382153     0.11972672

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 43  1    -94.6661949581  4.4764E-13  0.0000E+00  5.9716E-05  1.0000E-04   
 mr-sdci # 43  2    -94.3453199598  1.0252E-06  1.5881E-07  6.4269E-04  1.0000E-04   
 mr-sdci # 43  3    -94.3093370873  6.2030E-12  0.0000E+00  2.3682E-05  1.0000E-04   
 mr-sdci # 43  4    -93.9203665425  7.5187E-05  0.0000E+00  2.9593E-01  1.0000E-04   
 mr-sdci # 43  5    -93.9128961819  6.5497E-09  0.0000E+00  3.2315E-01  1.0000E-04   
 mr-sdci # 43  6    -92.0181227998 -1.4969E+00  0.0000E+00  1.6245E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000027
time for cinew                         0.005587
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  44

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67078302
   ht   3     0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24575544
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -40.23836026
   ht   6    -0.00107237     0.00001981     0.01349475    -0.00010250    -0.00065892    -0.00004006
   ht   7     0.10727965     0.00333000     0.02580524    -0.00911894    -0.00212796    -0.00000743    -0.00032837

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735      -2.458226E-05  -1.233984E-02   1.610067E-02  -1.987603E-02  -0.107526       0.222370    
 ref    2  -5.483007E-07  -0.947850       4.513070E-07  -1.215673E-03   1.183578E-03   5.867139E-03  -1.716440E-02
 ref    3   6.912433E-03  -2.665587E-05  -0.945584      -7.507196E-02  -3.870169E-03  -0.110978      -4.769427E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.898419       0.894281       5.896508E-03   4.114357E-04   2.391251E-02   5.201769E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873496    -0.00002458    -0.01233984     0.01610067    -0.01987603    -0.10752625     0.22236980
 ref:   2    -0.00000055    -0.94784975     0.00000045    -0.00121567     0.00118358     0.00586714    -0.01716440
 ref:   3     0.00691243    -0.00002666    -0.94558399    -0.07507196    -0.00387017    -0.11097835    -0.04769427

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 44  1    -94.6661949581  2.3448E-13  0.0000E+00  5.9745E-05  1.0000E-04   
 mr-sdci # 44  2    -94.3453199723  1.2489E-08 -3.9401E-07  6.5494E-04  1.0000E-04   
 mr-sdci # 44  3    -94.3093370873  2.8422E-14  0.0000E+00  2.3682E-05  1.0000E-04   
 mr-sdci # 44  4    -93.9203796760  1.3133E-05  0.0000E+00  2.9605E-01  1.0000E-04   
 mr-sdci # 44  5    -93.9128967863  6.0441E-07  0.0000E+00  3.2309E-01  1.0000E-04   
 mr-sdci # 44  6    -92.2377877858  2.1966E-01  0.0000E+00  1.6244E+00  1.0000E-04   
 mr-sdci # 44  7    -91.3974326844  3.6900E-02  0.0000E+00  1.4045E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000310
time for cinew                         0.007052
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration  45

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2     0.00000000   -40.67078302
   ht   3     0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24575544
   ht   5    -0.00000000     0.00000000     0.00000000    -0.00000000   -40.23836026
   ht   6    -0.00107237     0.00001981     0.01349475    -0.00010250    -0.00065892    -0.00004006
   ht   7     0.10727965     0.00333000     0.02580524    -0.00911894    -0.00212796    -0.00000743    -0.00032837
   ht   8    -0.28512018    -0.00286935    -0.08625102     0.02027972     0.00511211     0.00002500     0.00087046    -0.00234032

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.958735       8.186284E-06  -1.233984E-02   1.580712E-02  -1.972062E-02   6.405741E-03  -2.883780E-02  -0.265755    
 ref    2   5.818283E-07   0.947835       4.566412E-07  -1.791267E-03   1.285851E-03  -2.859430E-02   2.130752E-02   1.950811E-03
 ref    3  -6.912309E-03  -2.078564E-05  -0.945584      -7.658337E-02  -3.958058E-03  -5.145591E-02   0.170953      -8.544494E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.898392       0.894281       6.118086E-03   4.062226E-04   3.506378E-03   3.051040E-02   7.793046E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873501     0.00000819    -0.01233984     0.01580712    -0.01972062     0.00640574    -0.02883780    -0.26575518
 ref:   2     0.00000058     0.94783525     0.00000046    -0.00179127     0.00128585    -0.02859430     0.02130752     0.00195081
 ref:   3    -0.00691231    -0.00002079    -0.94558398    -0.07658337    -0.00395806    -0.05145591     0.17095252    -0.08544494

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 45  1    -94.6661949581  1.8687E-12  0.0000E+00  5.9782E-05  1.0000E-04   
 mr-sdci # 45  2    -94.3453202107  2.3845E-07  1.2420E-07  4.7751E-04  1.0000E-04   
 mr-sdci # 45  3    -94.3093370873 -2.1316E-14  0.0000E+00  2.3692E-05  1.0000E-04   
 mr-sdci # 45  4    -93.9205790508  1.9937E-04  0.0000E+00  2.9584E-01  1.0000E-04   
 mr-sdci # 45  5    -93.9129041159  7.3296E-06  0.0000E+00  3.2328E-01  1.0000E-04   
 mr-sdci # 45  6    -93.3769148091  1.1391E+00  0.0000E+00  1.1994E+00  1.0000E-04   
 mr-sdci # 45  7    -91.5316991316  1.3427E-01  0.0000E+00  1.5879E+00  1.0000E-04   
 mr-sdci # 45  8    -91.2860534531  3.9555E-02  0.0000E+00  1.2917E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000083
time for cinew                         0.007974
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  46

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.67078430
   ht   3     0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000    -0.00000000    -0.00000000   -40.24604314
   ht   5     0.00000000    -0.00000000    -0.00000000     0.00000000   -40.23836820
   ht   6    -0.00388054    -0.00262963     0.01047735    -0.00229234     0.00040883    -0.00000689

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735       7.257997E-06   1.233983E-02   1.575789E-02   1.971219E-02  -5.832548E-03
 ref    2   5.447085E-07   0.947847      -4.050861E-07  -1.154357E-03  -1.272640E-03   6.357876E-02
 ref    3  -6.912401E-03   8.959271E-06   0.945584      -7.492734E-02   4.025560E-03   0.163422    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.898414       0.894282       5.863750E-03   4.063950E-04   3.078300E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873501     0.00000726     0.01233983     0.01575789     0.01971219    -0.00583255
 ref:   2     0.00000054     0.94784711    -0.00000041    -0.00115436    -0.00127264     0.06357876
 ref:   3    -0.00691240     0.00000896     0.94558414    -0.07492734     0.00402556     0.16342193

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 46  1    -94.6661949581  9.0239E-13  0.0000E+00  5.9713E-05  1.0000E-04   
 mr-sdci # 46  2    -94.3453202873  7.6590E-08 -3.7860E-08  2.1723E-04  1.0000E-04   
 mr-sdci # 46  3    -94.3093370873  2.1103E-12  0.0000E+00  2.3640E-05  1.0000E-04   
 mr-sdci # 46  4    -93.9207620392  1.8299E-04  0.0000E+00  2.9475E-01  1.0000E-04   
 mr-sdci # 46  5    -93.9129041863  7.0339E-08  0.0000E+00  3.2329E-01  1.0000E-04   
 mr-sdci # 46  6    -92.1431558389 -1.2338E+00  0.0000E+00  1.7623E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000366
time for cinew                         0.007965
time for eigenvalue solver             0.000126
time for vector access                 0.000001

          starting ci iteration  47

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.67078430
   ht   3     0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000    -0.00000000    -0.00000000   -40.24604314
   ht   5     0.00000000    -0.00000000    -0.00000000     0.00000000   -40.23836820
   ht   6    -0.00388054    -0.00262963     0.01047735    -0.00229234     0.00040883    -0.00000689
   ht   7    -0.09111413     0.00069395     0.02731453     0.00666141    -0.00164940    -0.00001648    -0.00023878

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.958735      -2.172193E-06  -1.233981E-02   1.537793E-02  -1.967426E-02  -9.762975E-03   0.268942    
 ref    2  -5.437110E-07  -0.947847       4.042364E-07  -1.137540E-03   1.270930E-03  -6.267391E-02  -1.744389E-02
 ref    3   6.912398E-03  -8.038372E-06  -0.945584      -7.499632E-02  -4.023706E-03  -0.165709       3.473316E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.919221       0.898415       0.894282       5.862223E-03   4.048819E-04   3.148286E-02   7.384058E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.95873503    -0.00000217    -0.01233981     0.01537793    -0.01967426    -0.00976298     0.26894218
 ref:   2    -0.00000054    -0.94784734     0.00000040    -0.00113754     0.00127093    -0.06267391    -0.01744389
 ref:   3     0.00691240    -0.00000804    -0.94558413    -0.07499632    -0.00402371    -0.16570917     0.03473316

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 47  1    -94.6661949581 -1.4211E-14  0.0000E+00  5.9704E-05  1.0000E-04   
 mr-sdci # 47  2    -94.3453202884  1.0898E-09  1.5650E-08  2.0696E-04  1.0000E-04   
 mr-sdci # 47  3    -94.3093370873  2.8422E-14  0.0000E+00  2.3643E-05  1.0000E-04   
 mr-sdci # 47  4    -93.9207673026  5.2635E-06  0.0000E+00  2.9462E-01  1.0000E-04   
 mr-sdci # 47  5    -93.9129042355  4.9215E-08  0.0000E+00  3.2331E-01  1.0000E-04   
 mr-sdci # 47  6    -92.1460001179  2.8443E-03  0.0000E+00  1.7613E+00  1.0000E-04   
 mr-sdci # 47  7    -91.2983012130 -2.3340E-01  0.0000E+00  1.2999E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000705
time for cinew                         0.008248
time for eigenvalue solver             0.000130
time for vector access                 0.000001

          starting ci iteration  48

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7         ht   8
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.67078430
   ht   3     0.00000000     0.00000000   -40.63480117
   ht   4     0.00000000    -0.00000000    -0.00000000   -40.24604314
   ht   5     0.00000000    -0.00000000    -0.00000000     0.00000000   -40.23836820
   ht   6    -0.00388054    -0.00262963     0.01047735    -0.00229234     0.00040883    -0.00000689
   ht   7    -0.09111413     0.00069395     0.02731453     0.00666141    -0.00164940    -0.00001648    -0.00023878
   ht   8     0.03656365    -0.00077556    -0.00887304    -0.00320625     0.00075025     0.00000622     0.00009485    -0.00003812

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6         v:   7         v:   8

   eig(s)   6.735127E-09   8.856837E-08   5.263586E-07    1.00000        1.00000        1.00000        1.00000        1.00001    
 
   x:   1  -4.946537E-05   1.086585E-04  -2.393932E-03   7.708908E-02   7.169199E-02  -0.135621       0.241721      -0.955034    
   x:   2  -1.360753E-05  -6.585573E-05   1.768299E-05  -0.392055      -0.183413      -0.859079      -0.273094       7.459964E-03
   x:   3   7.185291E-05   1.950469E-04   7.230503E-04   8.240425E-02   0.109180      -0.339788       0.885097       0.287118    
   x:   4  -1.116959E-05  -7.160354E-05   1.780701E-04   0.563112       0.704909      -0.316558      -0.284075       7.142234E-02
   x:   5   1.604271E-06   1.384162E-05  -4.400291E-05  -0.718655       0.672610       0.167109       5.381104E-02  -1.762866E-02
   x:   6   4.195464E-02   0.995593       8.386624E-02   1.999247E-13  -1.303555E-11  -4.546022E-07  -2.396770E-04  -1.597210E-04
   x:   7   0.389660      -9.359805E-02   0.916187      -1.558240E-14   2.423188E-11   8.482077E-07  -3.795015E-06  -2.328470E-03
   x:   8   0.920003      -5.759084E-03  -0.391869        0.00000       5.902666E-11   2.066218E-06  -5.117496E-05   9.207239E-04
 bummer (warning):overlap matrix: # small eigenvalues=1

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1  -0.958735       1.049059E-06  -1.233980E-02   1.529614E-02  -1.971348E-02   3.609018E-03  -0.269308      -2.828624E-03
 ref    2  -5.429998E-07   0.947849       3.987874E-07  -9.806609E-04   1.321976E-03  -1.949701E-02   1.443846E-02  -6.653525E-02
 ref    3   6.912384E-03  -4.778978E-06  -0.945584      -7.548203E-02  -4.062248E-03   2.260618E-02  -4.499544E-02  -0.236616    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8
 ref    1   0.919221       0.898417       0.894281       5.932470E-03   4.068706E-04   9.041977E-04   7.475998E-02   6.042184E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7        ev    8

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7       v      8

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.95873503     0.00000105    -0.01233980     0.01529614    -0.01971348     0.00360902    -0.26930823    -0.00282862
 ref:   2    -0.00000054     0.94784877     0.00000040    -0.00098066     0.00132198    -0.01949701     0.01443846    -0.06653525
 ref:   3     0.00691238    -0.00000478    -0.94558404    -0.07548203    -0.00406225     0.02260618    -0.04499544    -0.23661551

 trial vector basis is being transformed.  new dimension:   5

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 48  1    -94.6661949581  4.2633E-14  0.0000E+00  5.9714E-05  1.0000E-04   
 mr-sdci # 48  2    -94.3453203176  2.9187E-08  1.1443E-08  1.5816E-04  1.0000E-04   
 mr-sdci # 48  3    -94.3093370873  1.3500E-12  0.0000E+00  2.3609E-05  1.0000E-04   
 mr-sdci # 48  4    -93.9208047356  3.7433E-05  0.0000E+00  2.9457E-01  1.0000E-04   
 mr-sdci # 48  5    -93.9129078609  3.6255E-06  0.0000E+00  3.2328E-01  1.0000E-04   
 mr-sdci # 48  6    -93.6413115108  1.4953E+00  0.0000E+00  9.2877E-01  1.0000E-04   
 mr-sdci # 48  7    -91.2987191147  4.1790E-04  0.0000E+00  1.2980E+00  1.0000E-04   
 mr-sdci # 48  8    -91.0423039084 -2.4375E-01  0.0000E+00  1.4316E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000106
time for cinew                         0.010283
time for eigenvalue solver             0.000164
time for vector access                 0.000000

          starting ci iteration  49

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:        6924 2x:        1665 4x:         261
All internal counts: zz :        1113 yy:        7060 xx:        7916 ww:       11868
One-external counts: yz :        2797 yx:        4713 yw:        5711
Two-external counts: yy :        1652 ww:        2588 xx:        1914 xz:         265 wz:         421 wx:        2412
Three-ext.   counts: yx :         335 yw:         401

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:        1587
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1   -40.99165905
   ht   2    -0.00000000   -40.67078440
   ht   3     0.00000000    -0.00000000   -40.63480117
   ht   4     0.00000000     0.00000000     0.00000000   -40.24626882
   ht   5    -0.00000000    -0.00000000    -0.00000000    -0.00000000   -40.23837195
   ht   6    -0.00004143     0.00044662     0.00121540    -0.00004240    -0.00004688    -0.00000031

                v:   1         v:   2         v:   3         v:   4         v:   5         v:   6

   eig(s)   7.091881E-09    1.00000        1.00000        1.00000        1.00000        1.00000    
 
   x:   1  -1.010758E-06   6.591216E-03  -6.687506E-02   0.997218      -6.146120E-03  -3.166919E-02
   x:   2   1.098159E-05  -0.910005      -0.224609       1.538767E-03  -5.522608E-02   0.344077    
   x:   3   2.991021E-05   0.343896       4.705097E-02   3.053048E-02  -1.840260E-02   0.937151    
   x:   4  -1.027141E-06   0.226940      -0.971000      -6.759405E-02   7.268103E-03  -3.218251E-02
   x:   5  -1.160155E-06   4.561580E-02   4.900658E-03  -7.279797E-03  -0.998259      -3.635057E-02
   x:   6    1.00000        0.00000        0.00000        0.00000       1.391363E-11  -3.191612E-05
 bummer (warning):overlap matrix: # small eigenvalues=1

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.958735      -1.706504E-08  -1.233979E-02  -1.513292E-02  -1.976758E-02  -1.841314E-02
 ref    2   5.444973E-07   0.947848       4.085257E-07   1.033905E-03   1.313393E-03  -8.070626E-03
 ref    3  -6.912407E-03   2.307841E-06  -0.945584       7.466998E-02  -3.755234E-03   0.122302    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.919221       0.898416       0.894282       5.805681E-03   4.065840E-04   1.536202E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.95873504    -0.00000002    -0.01233979    -0.01513292    -0.01976758    -0.01841314
 ref:   2     0.00000054     0.94784829     0.00000041     0.00103390     0.00131339    -0.00807063
 ref:   3    -0.00691241     0.00000231    -0.94558414     0.07466998    -0.00375523     0.12230225

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 49  1    -94.6661949581  1.2079E-13  0.0000E+00  5.9711E-05  1.0000E-04   
 mr-sdci # 49  2    -94.3453203257  8.0523E-09  1.4206E-09  6.0963E-05  1.0000E-04   
 mr-sdci # 49  3    -94.3093370873  1.5774E-12  0.0000E+00  2.3556E-05  1.0000E-04   
 mr-sdci # 49  4    -93.9208890913  8.4356E-05  0.0000E+00  2.9505E-01  1.0000E-04   
 mr-sdci # 49  5    -93.9129108923  3.0314E-06  0.0000E+00  3.2317E-01  1.0000E-04   
 mr-sdci # 49  6    -92.0518663385 -1.5894E+00  0.0000E+00  1.6127E+00  1.0000E-04   
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000080
time for cinew                         0.004796
time for eigenvalue solver             0.000018
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after 49 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 49  1    -94.6661949581  1.2079E-13  0.0000E+00  5.9711E-05  1.0000E-04   
 mr-sdci # 49  2    -94.3453203257  8.0523E-09  1.4206E-09  6.0963E-05  1.0000E-04   
 mr-sdci # 49  3    -94.3093370873  1.5774E-12  0.0000E+00  2.3556E-05  1.0000E-04   
 mr-sdci # 49  4    -93.9208890913  8.4356E-05  0.0000E+00  2.9505E-01  1.0000E-04   
 mr-sdci # 49  5    -93.9129108923  3.0314E-06  0.0000E+00  3.2317E-01  1.0000E-04   
 mr-sdci # 49  6    -92.0518663385 -1.5894E+00  0.0000E+00  1.6127E+00  1.0000E-04   

####################CIUDGINFO####################

   ci vector at position   1 energy=  -94.666194958137
   ci vector at position   2 energy=  -94.345320325663
   ci vector at position   3 energy=  -94.309337087262

################END OF CIUDGINFO################

 
    3 of the   7 expansion vectors are transformed.
    3 of the   6 matrix-vector products are transformed.

    3 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =       -94.6661949581

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1  0.906937                        +-   +-   +-   +-   +-   +-      
 z*  1  1       2  0.297233                        +-   +-   +-   +-   +-   +     - 
 z*  1  1       3 -0.092616                        +-   +-   +-   +-   +-        +- 
 z*  1  1       6 -0.011592                        +-   +-   +-   +-        +-   +- 
 z   1  1      10 -0.012284                        +-   +-   +-        +-   +-   +- 
 z   1  1      19  0.013628                        +-   +    +-    -   +-   +-   +- 
 z   1  1      21 -0.024936                        +-        +-   +-   +-   +-   +- 
 y   1  1      34  0.019056              6( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      65  0.039718             10( a  )   +-   +-   +-   +-   +-         - 
 y   1  1      73 -0.016598             18( a  )   +-   +-   +-   +-   +-         - 
 y   1  1     117 -0.032551              8( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     131 -0.010130             22( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     165 -0.012047              2( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     276  0.036176              5( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     278  0.032423              7( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     282  0.026657             11( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     286 -0.022951             15( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     409 -0.012179              3( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     599 -0.024568              4( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     609 -0.019636             14( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1    1082 -0.021562              1( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1084  0.011188              3( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1088  0.012774              7( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1092  0.010990             11( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1094 -0.014930             13( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1379 -0.015580              1( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1381  0.013476              3( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1383  0.014873              5( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1387 -0.020432              9( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1391 -0.012604             13( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1399 -0.016023             21( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1424 -0.012417             19( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1425  0.017207             20( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1477 -0.021895             18( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1559  0.014837             19( a  )   +-   +     -   +-   +-   +-    - 
 y   1  1    1730  0.014342              1( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1732  0.015822              3( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1734  0.012984              5( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1736  0.015629              7( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1740  0.014660             11( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    1745  0.014123             16( a  )    -   +-   +-   +-   +-   +     - 
 y   1  1    2108  0.013274              1( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2110  0.016841              3( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2114  0.012583              7( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2118  0.022992             11( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2131  0.012670             24( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2213  0.010040             25( a  )   +    +-   +-    -   +-   +-    - 
 y   1  1    2289 -0.012402             20( a  )   +    +-    -   +-   +-   +-    - 
 x   1  1    2526 -0.013880    4( a  )   6( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2573 -0.010757    6( a  )  12( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2692  0.010425    9( a  )  20( a  )   +-   +-   +-   +-    -    -      
 x   1  1    3586  0.015089    6( a  )   7( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3606  0.011037    5( a  )  10( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3666  0.012649   10( a  )  15( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3710 -0.010255    9( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    4275 -0.011951    2( a  )   5( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4286 -0.010176    4( a  )   7( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4293 -0.013988    5( a  )   8( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4295 -0.013713    7( a  )   8( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4366 -0.010617    8( a  )  15( a  )   +-   +-   +-    -    -   +-      
 x   1  1    6385 -0.011383    2( a  )   6( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6411 -0.013294    2( a  )  10( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6413 -0.011864    4( a  )  10( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6461  0.018066   10( a  )  14( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6537  0.012431   11( a  )  19( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6555  0.010352   11( a  )  20( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6582 -0.012557   19( a  )  21( a  )   +-   +-    -   +-   +-    -      
 x   1  1    7080  0.010343    2( a  )   4( a  )   +-   +-    -   +-    -   +-      
 x   1  1    9198  0.014530    2( a  )   7( a  )   +-   +-    -    -   +-   +-      
 x   1  1    9243 -0.011426    7( a  )  12( a  )   +-   +-    -    -   +-   +-      
 x   1  1   11430 -0.010804    7( a  )  18( a  )   +-    -   +-   +-   +-    -      
 x   1  1   18422 -0.010813   10( a  )  16( a  )    -   +-   +-   +-   +-    -      
 w   1  1   27454 -0.010578    6( a  )   6( a  )   +-   +-   +-   +-   +-           
 w   1  1   27484 -0.010124    6( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27488 -0.024439   10( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27604 -0.012626   18( a  )  18( a  )   +-   +-   +-   +-   +-           
 w   1  1   27623 -0.011105   19( a  )  19( a  )   +-   +-   +-   +-   +-           
 w   1  1   27830  0.011074    4( a  )   6( a  )   +-   +-   +-   +-   +     -      
 w   1  1   28570 -0.014938    2( a  )   2( a  )   +-   +-   +-   +-        +-      
 w   1  1   28573 -0.016965    3( a  )   3( a  )   +-   +-   +-   +-        +-      
 w   1  1   28575  0.015272    2( a  )   4( a  )   +-   +-   +-   +-        +-      
 w   1  1   28577 -0.016395    4( a  )   4( a  )   +-   +-   +-   +-        +-      
 w   1  1   28603 -0.016461    8( a  )   8( a  )   +-   +-   +-   +-        +-      
 w   1  1   28612 -0.014765    9( a  )   9( a  )   +-   +-   +-   +-        +-      
 w   1  1   28635 -0.017474    2( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   28637  0.015270    4( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   28645 -0.016729   12( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   28648  0.012704    3( a  )  13( a  )   +-   +-   +-   +-        +-      
 w   1  1   28654 -0.010767    9( a  )  13( a  )   +-   +-   +-   +-        +-      
 w   1  1   28662 -0.010753    4( a  )  14( a  )   +-   +-   +-   +-        +-      
 w   1  1   28670  0.010901   12( a  )  14( a  )   +-   +-   +-   +-        +-      
 w   1  1   29728  0.010097    6( a  )   7( a  )   +-   +-   +-   +    +-    -      
 w   1  1   29751 -0.010675    5( a  )  10( a  )   +-   +-   +-   +    +-    -      
 w   1  1   29766 -0.010105   10( a  )  11( a  )   +-   +-   +-   +    +-    -      
 w   1  1   29816  0.011938   10( a  )  15( a  )   +-   +-   +-   +    +-    -      
 w   1  1   30462 -0.016742    2( a  )   3( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30466  0.013549    3( a  )   4( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30501 -0.013313    8( a  )   9( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30520  0.010048    8( a  )  11( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30526 -0.013438    3( a  )  12( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30547  0.011493   12( a  )  13( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30561 -0.010665   13( a  )  14( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30564 -0.011708    2( a  )  15( a  )   +-   +-   +-   +     -   +-      
 w   1  1   30574 -0.012289   12( a  )  15( a  )   +-   +-   +-   +     -   +-      
 w   1  1   32726 -0.010274    1( a  )   1( a  )   +-   +-   +-        +-   +-      
 w   1  1   32731 -0.010873    3( a  )   3( a  )   +-   +-   +-        +-   +-      
 w   1  1   32740 -0.016088    5( a  )   5( a  )   +-   +-   +-        +-   +-      
 w   1  1   32753 -0.016450    7( a  )   7( a  )   +-   +-   +-        +-   +-      
 w   1  1   32835  0.011879    5( a  )  15( a  )   +-   +-   +-        +-   +-      
 w   1  1   32845 -0.013272   15( a  )  15( a  )   +-   +-   +-        +-   +-      
 w   1  1   32896 -0.010486   18( a  )  18( a  )   +-   +-   +-        +-   +-      
 w   1  1   35040  0.012685    2( a  )  10( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35042  0.011910    4( a  )  10( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35094  0.018415   10( a  )  14( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35757 -0.016030    2( a  )   4( a  )   +-   +-   +    +-    -   +-      
 w   1  1   35842 -0.016294    2( a  )  14( a  )   +-   +-   +    +-    -   +-      
 w   1  1   38019 -0.013270    1( a  )   2( a  )   +-   +-   +     -   +-   +-      
 w   1  1   38113  0.010187    5( a  )  14( a  )   +-   +-   +     -   +-   +-      
 w   1  1   38136 -0.014807   14( a  )  15( a  )   +-   +-   +     -   +-   +-      
 w   1  1   40286 -0.012407    1( a  )   1( a  )   +-   +-        +-   +-   +-      
 w   1  1   40288 -0.020007    2( a  )   2( a  )   +-   +-        +-   +-   +-      
 w   1  1   40293 -0.011791    2( a  )   4( a  )   +-   +-        +-   +-   +-      
 w   1  1   40351 -0.012210   11( a  )  11( a  )   +-   +-        +-   +-   +-      
 w   1  1   40353 -0.011242    2( a  )  12( a  )   +-   +-        +-   +-   +-      
 w   1  1   40378 -0.021544    2( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1   40380 -0.014456    4( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1   40390 -0.021031   14( a  )  14( a  )   +-   +-        +-   +-   +-      
 w   1  1   42571  0.012033    3( a  )   6( a  )   +-   +    +-   +-   +-    -      
 w   1  1   42599  0.011021    1( a  )  10( a  )   +-   +    +-   +-   +-    -      
 w   1  1   43314  0.010031    2( a  )   3( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43318 -0.017128    3( a  )   4( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43353  0.015037    8( a  )   9( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43378  0.014321    3( a  )  12( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43389 -0.012461    2( a  )  13( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43399 -0.015794   12( a  )  13( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43403 -0.011236    3( a  )  14( a  )   +-   +    +-   +-    -   +-      
 w   1  1   45581 -0.013835    1( a  )   3( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45583  0.010700    3( a  )   3( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45588  0.012281    1( a  )   5( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45601 -0.010822    3( a  )   7( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45658 -0.015344    3( a  )  13( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45683 -0.016329    1( a  )  15( a  )   +-   +    +-    -   +-   +-      
 w   1  1   45695 -0.019050   13( a  )  15( a  )   +-   +    +-    -   +-   +-      
 w   1  1   47847 -0.012245    1( a  )   2( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47850 -0.011472    2( a  )   3( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47937 -0.013363    1( a  )  14( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47949 -0.010636   13( a  )  14( a  )   +-   +     -   +-   +-   +-      
 w   1  1   47964 -0.012347   14( a  )  15( a  )   +-   +     -   +-   +-   +-      
 w   1  1   50119 -0.010620    3( a  )   3( a  )   +-        +-   +-   +-   +-      
 w   1  1   50204 -0.010100   13( a  )  13( a  )   +-        +-   +-   +-   +-      
 w   1  1   50233 -0.010286   15( a  )  15( a  )   +-        +-   +-   +-   +-      
 w   1  1   52427 -0.013943    1( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52429 -0.010498    3( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52431 -0.010370    5( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52446 -0.011368   10( a  )  11( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52511 -0.014717   10( a  )  16( a  )   +    +-   +-   +-   +-    -      
 w   1  1   57765  0.014435    1( a  )  14( a  )   +    +-    -   +-   +-   +-      
 w   1  1   57775  0.011989   11( a  )  14( a  )   +    +-    -   +-   +-   +-      
 w   1  1   57807  0.014397   14( a  )  16( a  )   +    +-    -   +-   +-   +-      

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             157
     0.01> rq > 0.001           3808
    0.001> rq > 0.0001          6452
   0.0001> rq > 0.00001         5167
  0.00001> rq > 0.000001        2539
 0.000001> rq                  46352
           all                 64477
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:         217
  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.906937267250    -85.613198970540
     2     2      0.297233060726    -28.053291552447
     3     3     -0.092615707866      8.758210955918
     4     4      0.000003808257     -0.000359448506
     5     5     -0.000003248050      0.000307384666
     6     6     -0.011592116111      1.098938800584

 number of reference csfs (nref) is     6.  root number (iroot) is  1.
 c0**2 =   0.91959475  c**2 (all zwalks) =   0.92070389

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.398161553085   "relaxed" cnot**2         =   0.919594745640
 eci       =    -94.666194958137   deltae = eci - eref       =  -0.268033405052
 eci+dv1   =    -94.687746252247   dv1 = (1-cnot**2)*deltae  =  -0.021551294110
 eci+dv2   =    -94.689630601079   dv2 = dv1 / cnot**2       =  -0.023435642942
 eci+dv3   =    -94.691876040103   dv3 = dv1 / (2*cnot**2-1) =  -0.025681081965
 eci+pople =    -94.686908747252   ( 12e- scaled deltae )    =  -0.288747194167


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 2) =       -94.3453203257

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       4 -0.917478                        +-   +-   +-   +-   +    +-    - 
 z*  1  1       5 -0.238069                        +-   +-   +-   +-   +     -   +- 
 z   1  1      11 -0.082409                        +-   +-   +    +-   +-   +-    - 
 z   1  1      12  0.066245                        +-   +-   +    +-   +-    -   +- 
 y   1  1      30  0.021101              2( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      88 -0.017977              6( a  )   +-   +-   +-   +-    -   +-      
 y   1  1      92  0.017985             10( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     115  0.027599              6( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     119  0.011422             10( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     127 -0.026429             18( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     146  0.028208             10( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     154 -0.013962             18( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     173  0.014109             10( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     192  0.018193              2( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     194 -0.021739              4( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     198 -0.021769              8( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     202  0.019793             12( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     204 -0.017303             14( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     207 -0.018179             17( a  )   +-   +-   +-   +-        +-    - 
 y   1  1     225  0.012452              8( a  )   +-   +-   +-   +-         -   +- 
 y   1  1     330  0.021323              5( a  )   +-   +-   +-    -   +    +-    - 
 y   1  1     336  0.010653             11( a  )   +-   +-   +-    -   +    +-    - 
 y   1  1     357 -0.023208              5( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     359 -0.024211              7( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     363 -0.017866             11( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     367  0.017263             15( a  )   +-   +-   +-    -   +     -   +- 
 y   1  1     434 -0.020455              1( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     436  0.030981              3( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     442 -0.041899              9( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     446 -0.026483             13( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     448  0.016434             15( a  )   +-   +-   +-   +     -   +-    - 
 y   1  1     495 -0.010727              8( a  )   +-   +-   +-        +-   +-    - 
 y   1  1     586 -0.012713             18( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     653 -0.013080              4( a  )   +-   +-    -   +-   +    +-    - 
 y   1  1     680  0.018023              4( a  )   +-   +-    -   +-   +     -   +- 
 y   1  1     690  0.014344             14( a  )   +-   +-    -   +-   +     -   +- 
 y   1  1     735  0.014918              5( a  )   +-   +-    -   +    +-   +-    - 
 y   1  1     821 -0.013349             10( a  )   +-   +-   +    +-   +-    -    - 
 y   1  1     840 -0.010732              2( a  )   +-   +-   +    +-    -   +-    - 
 y   1  1     846  0.025426              8( a  )   +-   +-   +    +-    -   +-    - 
 y   1  1     897 -0.014353              5( a  )   +-   +-   +     -   +-   +-    - 
 y   1  1     899 -0.017296              7( a  )   +-   +-   +     -   +-   +-    - 
 y   1  1     903 -0.018650             11( a  )   +-   +-   +     -   +-   +-    - 
 y   1  1     977  0.012881              4( a  )   +-   +-        +-   +-   +-    - 
 y   1  1     981  0.010110              8( a  )   +-   +-        +-   +-   +-    - 
 y   1  1    1136 -0.010457              1( a  )   +-    -   +-   +-   +    +-    - 
 y   1  1    1163  0.014621              1( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1165 -0.010248              3( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1169 -0.010017              7( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1175  0.012604             13( a  )   +-    -   +-   +-   +     -   +- 
 y   1  1    1298  0.013635              1( a  )   +-    -   +    +-   +-   +-    - 
 y   1  1    1304 -0.010356              7( a  )   +-    -   +    +-   +-   +-    - 
 y   1  1    1397 -0.019931             19( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1408 -0.038910              3( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1410  0.011791              5( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1414  0.041654              9( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1416 -0.022321             11( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1418  0.017822             13( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1420 -0.027433             15( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1426  0.011830             21( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1429 -0.013170             24( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1437  0.012042              5( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    1461 -0.017794              2( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1467  0.022413              8( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1471 -0.017302             12( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1531  0.013939             18( a  )   +-   +    +-    -    -   +-   +- 
 y   1  1    1551 -0.013899             11( a  )   +-   +     -   +-   +-   +-    - 
 y   1  1    1613 -0.011170             19( a  )   +-   +     -   +-    -   +-   +- 
 y   1  1    1784  0.010516              1( a  )    -   +-   +-   +-   +    +-    - 
 y   1  1    1813 -0.011095              3( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1815 -0.010029              5( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1817 -0.012118              7( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1821 -0.010058             11( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    1826 -0.010446             16( a  )    -   +-   +-   +-   +     -   +- 
 y   1  1    2127  0.010341             20( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2135  0.011898              1( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2143  0.016192              9( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2147  0.011335             13( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2155  0.010428             21( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2162  0.011393              1( a  )   +    +-   +-   +-    -    -   +- 
 y   1  1    2172  0.014944             11( a  )   +    +-   +-   +-    -    -   +- 
 x   1  1    4288 -0.013925    6( a  )   7( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4408 -0.013091    5( a  )  18( a  )   +-   +-   +-    -    -   +-      
 x   1  1    4414 -0.010904   11( a  )  18( a  )   +-   +-   +-    -    -   +-      
 x   1  1    5341 -0.012732    6( a  )   7( a  )   +-   +-   +-    -   +     -    - 
 x   1  1    5421 -0.010743   10( a  )  15( a  )   +-   +-   +-    -   +     -    - 
 x   1  1    7514  0.011946   10( a  )  14( a  )   +-   +-    -   +-    -   +     - 
 x   1  1    8166  0.011013    2( a  )  10( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8168  0.010099    4( a  )  10( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8216 -0.014932   10( a  )  14( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8292 -0.010394   11( a  )  19( a  )   +-   +-    -   +-   +     -    - 
 x   1  1    8337  0.010460   19( a  )  21( a  )   +-   +-    -   +-   +     -    - 
 x   1  1   10251 -0.014955    2( a  )   7( a  )   +-   +-    -    -   +    +-    - 
 x   1  1   10296  0.011835    7( a  )  12( a  )   +-   +-    -    -   +    +-    - 
 w   1  1   27832  0.023885    6( a  )   6( a  )   +-   +-   +-   +-   +     -      
 w   1  1   27856  0.010751    9( a  )   9( a  )   +-   +-   +-   +-   +     -      
 w   1  1   27862 -0.012056    6( a  )  10( a  )   +-   +-   +-   +-   +     -      
 w   1  1   27866 -0.013648   10( a  )  10( a  )   +-   +-   +-   +-   +     -      
 w   1  1   28020 -0.012098   19( a  )  20( a  )   +-   +-   +-   +-   +     -      
 w   1  1   28210 -0.012561    6( a  )   6( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28240 -0.011519    6( a  )  10( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28244 -0.028672   10( a  )  10( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28255 -0.011357   11( a  )  11( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28360 -0.014009   18( a  )  18( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28379 -0.011976   19( a  )  19( a  )   +-   +-   +-   +-   +          - 
 w   1  1   28584  0.014961    2( a  )   6( a  )   +-   +-   +-   +-        +-      
 w   1  1   28586 -0.017269    4( a  )   6( a  )   +-   +-   +-   +-        +-      
 w   1  1   28601  0.010446    6( a  )   8( a  )   +-   +-   +-   +-        +-      
 w   1  1   28639  0.017212    6( a  )  12( a  )   +-   +-   +-   +-        +-      
 w   1  1   30484 -0.012918    6( a  )   7( a  )   +-   +-   +-   +     -   +-      
 w   1  1   31974  0.010798    2( a  )   3( a  )   +-   +-   +-   +         +-    - 
 w   1  1   33860  0.010853    1( a  )   1( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33865  0.012147    3( a  )   3( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33874  0.016377    5( a  )   5( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33887  0.016603    7( a  )   7( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33938  0.010320    1( a  )  13( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33940 -0.010074    3( a  )  13( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33950  0.010360   13( a  )  13( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33967  0.010349    3( a  )  15( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33969 -0.011742    5( a  )  15( a  )   +-   +-   +-        +    +-    - 
 w   1  1   33979  0.013281   15( a  )  15( a  )   +-   +-   +-        +    +-    - 
 w   1  1   34030  0.010316   18( a  )  18( a  )   +-   +-   +-        +    +-    - 
 w   1  1   36228  0.012014   10( a  )  14( a  )   +-   +-   +    +-    -   +     - 
 w   1  1   36930 -0.010247    2( a  )  10( a  )   +-   +-   +    +-   +     -    - 
 w   1  1   36984 -0.014432   10( a  )  14( a  )   +-   +-   +    +-   +     -    - 
 w   1  1   37269  0.012659    2( a  )   4( a  )   +-   +-   +    +-        +-    - 
 w   1  1   37354  0.011505    2( a  )  14( a  )   +-   +-   +    +-        +-    - 
 w   1  1   39153  0.013652    1( a  )   2( a  )   +-   +-   +     -   +    +-    - 
 w   1  1   39231  0.010062    2( a  )  13( a  )   +-   +-   +     -   +    +-    - 
 w   1  1   39270  0.014056   14( a  )  15( a  )   +-   +-   +     -   +    +-    - 
 w   1  1   41420  0.012540    1( a  )   1( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41422  0.020654    2( a  )   2( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41427  0.010338    2( a  )   4( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41455  0.011329    8( a  )   8( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41485  0.012065   11( a  )  11( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41487  0.011817    2( a  )  12( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41512  0.019935    2( a  )  14( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41514  0.013790    4( a  )  14( a  )   +-   +-        +-   +    +-    - 
 w   1  1   41524  0.019684   14( a  )  14( a  )   +-   +-        +-   +    +-    - 
 w   1  1   43327 -0.014341    3( a  )   6( a  )   +-   +    +-   +-    -   +-      
 w   1  1   43393  0.012064    6( a  )  13( a  )   +-   +    +-   +-    -   +-      
 w   1  1   44461 -0.011035    3( a  )   6( a  )   +-   +    +-   +-   +     -    - 
 w   1  1   44830  0.011359    3( a  )   4( a  )   +-   +    +-   +-        +-    - 
 w   1  1   44865 -0.010910    8( a  )   9( a  )   +-   +    +-   +-        +-    - 
 w   1  1   44911  0.010536   12( a  )  13( a  )   +-   +    +-   +-        +-    - 
 w   1  1   46715  0.014916    1( a  )   3( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46717 -0.012308    3( a  )   3( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46722 -0.012701    1( a  )   5( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46735  0.010884    3( a  )   7( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46792  0.016633    3( a  )  13( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46817  0.016869    1( a  )  15( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   46829  0.019717   13( a  )  15( a  )   +-   +    +-    -   +    +-    - 
 w   1  1   48981  0.012222    1( a  )   2( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   48984  0.012811    2( a  )   3( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49042 -0.011171    8( a  )  11( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49071  0.012870    1( a  )  14( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49083  0.010422   13( a  )  14( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   49098  0.011684   14( a  )  15( a  )   +-   +     -   +-   +    +-    - 
 w   1  1   51248  0.010034    1( a  )   1( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51253  0.012383    3( a  )   3( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51326  0.010181    1( a  )  13( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51338  0.010574   13( a  )  13( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51355  0.010341    3( a  )  15( a  )   +-        +-   +-   +    +-    - 
 w   1  1   51367  0.010923   15( a  )  15( a  )   +-        +-   +-   +    +-    - 
 w   1  1   54317  0.011906    1( a  )  10( a  )   +    +-   +-   +-   +     -    - 
 w   1  1   54401  0.012992   10( a  )  16( a  )   +    +-   +-   +-   +     -    - 
 w   1  1   58899 -0.014335    1( a  )  14( a  )   +    +-    -   +-   +    +-    - 
 w   1  1   58909 -0.011873   11( a  )  14( a  )   +    +-    -   +-   +    +-    - 
 w   1  1   58941 -0.014101   14( a  )  16( a  )   +    +-    -   +-   +    +-    - 
 w   1  1   61181 -0.010171    1( a  )  15( a  )   +     -   +-   +-   +    +-    - 

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             168
     0.01> rq > 0.001           4014
    0.001> rq > 0.0001          6462
   0.0001> rq > 0.00001         3746
  0.00001> rq > 0.000001        2388
 0.000001> rq                  47697
           all                 64477
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:         217
  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.000021023972      0.001977704405
     2     2      0.000078597601     -0.007391734976
     3     3      0.000010686373     -0.001000846574
     4     4     -0.917478017220     86.291700825736
     5     5     -0.238069310028     22.392671026471
     6     6     -0.000004878201      0.000458968145

 number of reference csfs (nref) is     6.  root number (iroot) is  2.
 c0**2 =   0.89844292  c**2 (all zwalks) =   0.90962261

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.053551451474   "relaxed" cnot**2         =   0.898442915217
 eci       =    -94.345320325663   deltae = eci - eref       =  -0.291768874188
 eci+dv1   =    -94.374951521956   dv1 = (1-cnot**2)*deltae  =  -0.029631196293
 eci+dv2   =    -94.378300936663   dv2 = dv1 / cnot**2       =  -0.032980611000
 eci+dv3   =    -94.382504066622   dv3 = dv1 / (2*cnot**2-1) =  -0.037183740959
 eci+pople =    -94.374982644689   ( 12e- scaled deltae )    =  -0.321431193215


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 3) =       -94.3093370873

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7

                                          orbital     3    4    5    6    7    8    9

                                         symmetry   a    a    a    a    a    a    a  

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1 -0.304306                        +-   +-   +-   +-   +-   +-      
 z*  1  1       2  0.895563                        +-   +-   +-   +-   +-   +     - 
 z*  1  1       3 -0.063584                        +-   +-   +-   +-   +-        +- 
 z*  1  1       6  0.018174                        +-   +-   +-   +-        +-   +- 
 z   1  1      13 -0.012612                        +-   +-   +    +-    -   +-   +- 
 z   1  1      19 -0.014822                        +-   +    +-    -   +-   +-   +- 
 z   1  1      21  0.022389                        +-        +-   +-   +-   +-   +- 
 z   1  1      27  0.037088                        +     -   +-   +-   +-   +-   +- 
 y   1  1      38  0.018742             10( a  )   +-   +-   +-   +-   +-    -      
 y   1  1      65  0.026071             10( a  )   +-   +-   +-   +-   +-         - 
 y   1  1      73  0.011794             18( a  )   +-   +-   +-   +-   +-         - 
 y   1  1      86  0.010409              4( a  )   +-   +-   +-   +-    -   +-      
 y   1  1      90 -0.046369              8( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     104 -0.011774             22( a  )   +-   +-   +-   +-    -   +-      
 y   1  1     111  0.012711              2( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     117 -0.027003              8( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     121  0.012006             12( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     123 -0.010574             14( a  )   +-   +-   +-   +-    -   +     - 
 y   1  1     144 -0.018731              8( a  )   +-   +-   +-   +-    -        +- 
 y   1  1     165 -0.010045              2( a  )   +-   +-   +-   +-   +     -    - 
 y   1  1     245  0.016962              1( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     249  0.061375              5( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     251  0.051856              7( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     255  0.035186             11( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     259 -0.032324             15( a  )   +-   +-   +-    -   +-   +-      
 y   1  1     272  0.010312              1( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     276  0.032454              5( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     278  0.012554              7( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     280 -0.011166              9( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     282  0.016807             11( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     286 -0.013573             15( a  )   +-   +-   +-    -   +-   +     - 
 y   1  1     303  0.020570              5( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     305  0.020782              7( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     309  0.016488             11( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     313 -0.015645             15( a  )   +-   +-   +-    -   +-        +- 
 y   1  1     407  0.010699              1( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     409 -0.016106              3( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     415  0.012598              9( a  )   +-   +-   +-   +    +-    -    - 
 y   1  1     572 -0.044911              4( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     576 -0.013168              8( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     580  0.011166             12( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     582 -0.021045             14( a  )   +-   +-    -   +-   +-   +-      
 y   1  1     597 -0.017365              2( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     599 -0.013402              4( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     609 -0.016564             14( a  )   +-   +-    -   +-   +-   +     - 
 y   1  1     626 -0.016251              4( a  )   +-   +-    -   +-   +-        +- 
 y   1  1     636 -0.013185             14( a  )   +-   +-    -   +-   +-        +- 
 y   1  1    1055 -0.038833              1( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1057  0.020935              3( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1059  0.021235              5( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1061  0.016466              7( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1065  0.013947             11( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1067 -0.021558             13( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1069 -0.010168             15( a  )   +-    -   +-   +-   +-   +-      
 y   1  1    1082 -0.019427              1( a  )   +-    -   +-   +-   +-   +     - 
 y   1  1    1109 -0.014743              1( a  )   +-    -   +-   +-   +-        +- 
 y   1  1    1121 -0.011302             13( a  )   +-    -   +-   +-   +-        +- 
 y   1  1    1379 -0.015127              1( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1381  0.018800              3( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1387 -0.014293              9( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1391 -0.010851             13( a  )   +-   +    +-   +-   +-    -    - 
 y   1  1    1424  0.013840             19( a  )   +-   +    +-   +-    -   +-    - 
 y   1  1    1451 -0.010641             19( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    1452  0.012075             20( a  )   +-   +    +-   +-    -    -   +- 
 y   1  1    1477  0.010189             18( a  )   +-   +    +-    -   +-   +-    - 
 y   1  1    1504 -0.014432             18( a  )   +-   +    +-    -   +-    -   +- 
 y   1  1    1703  0.022927              1( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1705  0.019470              3( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1707  0.014505              5( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1709  0.018768              7( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1713  0.013511             11( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1718  0.017658             16( a  )    -   +-   +-   +-   +-   +-      
 y   1  1    1757  0.010381              1( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    1759  0.010682              3( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    1763  0.011017              7( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    1767  0.010271             11( a  )    -   +-   +-   +-   +-        +- 
 y   1  1    2110  0.015512              3( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2114  0.017255              7( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2118  0.023227             11( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2131  0.011363             24( a  )   +    +-   +-   +-   +-    -    - 
 y   1  1    2154 -0.012383             20( a  )   +    +-   +-   +-    -   +-    - 
 y   1  1    2206  0.021116             18( a  )   +    +-   +-    -   +-   +-    - 
 y   1  1    2288 -0.018696             19( a  )   +    +-    -   +-   +-   +-    - 
 y   1  1    2368  0.011180             18( a  )   +     -   +-   +-   +-   +-    - 
 x   1  1    2526 -0.010990    4( a  )   6( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2552  0.010368    4( a  )  10( a  )   +-   +-   +-   +-    -    -      
 x   1  1    2656 -0.010287    8( a  )  18( a  )   +-   +-   +-   +-    -    -      
 x   1  1    3706  0.012441    5( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3708  0.013599    7( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    3712  0.014154   11( a  )  18( a  )   +-   +-   +-    -   +-    -      
 x   1  1    4059 -0.010088    7( a  )  18( a  )   +-   +-   +-    -   +-         - 
 x   1  1    4626 -0.011809    2( a  )   5( a  )   +-   +-   +-    -    -   +     - 
 x   1  1    4644 -0.011101    5( a  )   8( a  )   +-   +-   +-    -    -   +     - 
 x   1  1    4646 -0.010958    7( a  )   8( a  )   +-   +-   +-    -    -   +     - 
 x   1  1    6411  0.014315    2( a  )  10( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6461 -0.013394   10( a  )  14( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6537 -0.010624   11( a  )  19( a  )   +-   +-    -   +-   +-    -      
 x   1  1    6762 -0.011076    2( a  )  10( a  )   +-   +-    -   +-   +-         - 
 x   1  1    6812  0.014785   10( a  )  14( a  )   +-   +-    -   +-   +-         - 
 x   1  1    6888  0.010050   11( a  )  19( a  )   +-   +-    -   +-   +-         - 
 x   1  1    9549  0.014261    2( a  )   7( a  )   +-   +-    -    -   +-   +     - 
 x   1  1    9594 -0.010839    7( a  )  12( a  )   +-   +-    -    -   +-   +     - 
 w   1  1   27434  0.010507    1( a  )   1( a  )   +-   +-   +-   +-   +-           
 w   1  1   27454 -0.011939    6( a  )   6( a  )   +-   +-   +-   +-   +-           
 w   1  1   27484  0.016138    6( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27488  0.034363   10( a  )  10( a  )   +-   +-   +-   +-   +-           
 w   1  1   27489  0.012143    1( a  )  11( a  )   +-   +-   +-   +-   +-           
 w   1  1   27499  0.011733   11( a  )  11( a  )   +-   +-   +-   +-   +-           
 w   1  1   27642  0.015840   19( a  )  20( a  )   +-   +-   +-   +-   +-           
 w   1  1   27654  0.012609   11( a  )  21( a  )   +-   +-   +-   +-   +-           
 w   1  1   27751  0.011567   18( a  )  25( a  )   +-   +-   +-   +-   +-           
 w   1  1   28948 -0.015162    2( a  )   2( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28951 -0.016461    3( a  )   3( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28953  0.014634    2( a  )   4( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28955 -0.015979    4( a  )   4( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28981 -0.016011    8( a  )   8( a  )   +-   +-   +-   +-        +     - 
 w   1  1   28990 -0.014280    9( a  )   9( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29013 -0.017202    2( a  )  12( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29015  0.014483    4( a  )  12( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29023 -0.016262   12( a  )  12( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29026  0.011667    3( a  )  13( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29032 -0.010137    9( a  )  13( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29040 -0.010592    4( a  )  14( a  )   +-   +-   +-   +-        +     - 
 w   1  1   29048  0.010288   12( a  )  14( a  )   +-   +-   +-   +-        +     - 
 w   1  1   30840 -0.016857    2( a  )   3( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30844  0.013165    3( a  )   4( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30879 -0.013530    8( a  )   9( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30904 -0.013118    3( a  )  12( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30925  0.010922   12( a  )  13( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30939 -0.010441   13( a  )  14( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30942 -0.011255    2( a  )  15( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   30952 -0.011556   12( a  )  15( a  )   +-   +-   +-   +     -   +     - 
 w   1  1   33104 -0.010320    1( a  )   1( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33109 -0.010701    3( a  )   3( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33118 -0.013486    5( a  )   5( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33131 -0.014019    7( a  )   7( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33213  0.010335    5( a  )  15( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33223 -0.012691   15( a  )  15( a  )   +-   +-   +-        +-   +     - 
 w   1  1   33274 -0.010086   18( a  )  18( a  )   +-   +-   +-        +-   +     - 
 w   1  1   35040 -0.012705    2( a  )  10( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35094 -0.013481   10( a  )  14( a  )   +-   +-   +    +-   +-    -      
 w   1  1   35418  0.011074    2( a  )  10( a  )   +-   +-   +    +-   +-         - 
 w   1  1   35472  0.014843   10( a  )  14( a  )   +-   +-   +    +-   +-         - 
 w   1  1   36135 -0.015516    2( a  )   4( a  )   +-   +-   +    +-    -   +     - 
 w   1  1   36220 -0.016227    2( a  )  14( a  )   +-   +-   +    +-    -   +     - 
 w   1  1   38397 -0.012592    1( a  )   2( a  )   +-   +-   +     -   +-   +     - 
 w   1  1   38514 -0.014781   14( a  )  15( a  )   +-   +-   +     -   +-   +     - 
 w   1  1   40664 -0.011120    1( a  )   1( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40666 -0.018528    2( a  )   2( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40671 -0.010433    2( a  )   4( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40729 -0.011417   11( a  )  11( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40731 -0.011058    2( a  )  12( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40756 -0.021054    2( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40758 -0.013530    4( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1   40768 -0.021174   14( a  )  14( a  )   +-   +-        +-   +-   +     - 
 w   1  1   42571  0.011016    3( a  )   6( a  )   +-   +    +-   +-   +-    -      
 w   1  1   42707  0.010498    1( a  )  18( a  )   +-   +    +-   +-   +-    -      
 w   1  1   43692  0.010304    2( a  )   3( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43696 -0.017092    3( a  )   4( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43731  0.015322    8( a  )   9( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43756  0.013971    3( a  )  12( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43767 -0.012382    2( a  )  13( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43777 -0.015349   12( a  )  13( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   43781 -0.011192    3( a  )  14( a  )   +-   +    +-   +-    -   +     - 
 w   1  1   45959 -0.013518    1( a  )   3( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   45961  0.010836    3( a  )   3( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   45966  0.011261    1( a  )   5( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   46036 -0.014845    3( a  )  13( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   46061 -0.015481    1( a  )  15( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   46073 -0.018063   13( a  )  15( a  )   +-   +    +-    -   +-   +     - 
 w   1  1   48225 -0.010721    1( a  )   2( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48228 -0.011004    2( a  )   3( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48315 -0.012281    1( a  )  14( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48327 -0.010248   13( a  )  14( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   48342 -0.012130   14( a  )  15( a  )   +-   +     -   +-   +-   +     - 
 w   1  1   50497 -0.010714    3( a  )   3( a  )   +-        +-   +-   +-   +     - 
 w   1  1   52427  0.011996    1( a  )  10( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52511  0.011394   10( a  )  16( a  )   +    +-   +-   +-   +-    -      
 w   1  1   52805 -0.010791    1( a  )  10( a  )   +    +-   +-   +-   +-         - 
 w   1  1   52889 -0.011070   10( a  )  16( a  )   +    +-   +-   +-   +-         - 
 w   1  1   58143  0.013436    1( a  )  14( a  )   +    +-    -   +-   +-   +     - 
 w   1  1   58153  0.011432   11( a  )  14( a  )   +    +-    -   +-   +-   +     - 
 w   1  1   58185  0.014225   14( a  )  16( a  )   +    +-    -   +-   +-   +     - 

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01             181
     0.01> rq > 0.001           4078
    0.001> rq > 0.0001          6670
   0.0001> rq > 0.00001         5242
  0.00001> rq > 0.000001        2284
 0.000001> rq                  46020
           all                 64477
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:         735 2x:           0 4x:           0
All internal counts: zz :        1113 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:         421    task #     2:         673    task #     3:        1461    task #     4:         730
task #     5:        2305    task #     6:        4013    task #     7:        4799    task #     8:        1047
task #     9:        7064    task #    10:        1039    task #    11:        1039    task #    12:        1373
task #    13:        1373    task #    14:         686    task #    15:         688    task #    16:         688
task #    17:         524    task #    18:        2022    task #    19:        1204    task #    20:         217
  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.304306297808     28.615529805360
     2     2      0.895562591233    -84.188320491329
     3     3     -0.063584163903      6.009288826433
     4     4      0.000080443520     -0.007562604159
     5     5      0.000014494617     -0.001361303339
     6     6      0.018174089160     -1.713418869980

 number of reference csfs (nref) is     6.  root number (iroot) is  3.
 c0**2 =   0.89900793  c**2 (all zwalks) =   0.90132867

 pople ci energy extrapolation is computed with 12 correlated electrons.

 eref      =    -94.011442873147   "relaxed" cnot**2         =   0.899007927798
 eci       =    -94.309337087262   deltae = eci - eref       =  -0.297894214115
 eci+dv1   =    -94.339422041243   dv1 = (1-cnot**2)*deltae  =  -0.030084953980
 eci+dv2   =    -94.342801702034   dv2 = dv1 / cnot**2       =  -0.033464614772
 eci+dv3   =    -94.347036781785   dv3 = dv1 / (2*cnot**2-1) =  -0.037699694523
 eci+pople =    -94.339420367954   ( 12e- scaled deltae )    =  -0.327977494807
maximum overlap with reference    1(overlap= 0.95874)
weight of reference states=  0.9192

 information on vector: 1 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    2(overlap= 0.94785)
weight of reference states=  0.8984

 information on vector: 2 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    3(overlap= 0.94558)
weight of reference states=  0.8943

 information on vector: 3 from unit 11 written to unit 48 filename civout              
 passed aftci ... 
 readint2: molcas,dalton2=                     0                     0
 files%faoints=aoints              
                       Size (real*8) of d2temp for two-external contributions      28539
 
                       Size (real*8) of d2temp for all-internal contributions        546
                       Size (real*8) of d2temp for one-external contributions       6804
                       Size (real*8) of d2temp for two-external contributions      28539
size_thrext:  lsym   l1    ksym   k1strt   k1       cnt3 
                1    1    1    1   28    10881
                1    2    1    1   28    10881
                1    3    1    1   28    10881
                1    4    1    1   28    10881
                1    5    1    1   28    10881
                1    6    1    1   28    10881
                1    7    1    1   28    10881
                       Size (real*8) of d2temp for three-external contributions      76167
                       Size (real*8) of d2temp for four-external contributions      82215
 enough memory for temporary d2 elements on vdisk ... 
location of d2temp files... fileloc(dd012)=       1
location of d2temp files... fileloc(d3)=       1
location of d2temp files... fileloc(d4)=       1
 files%dd012ext =  unit=  22  vdsk=   1  filestart=       1
 files%d3ext =     unit=  23  vdsk=   1  filestart=   55297
 files%d4ext =     unit=  24  vdsk=   1  filestart=  175297
            0xdiag    0ext      1ext      2ext      3ext      4ext
d2off                  6145     12289     24577         1         1
d2rec                     1         2         5         2         2
recsize                6144      6144      6144     60000     60000
d2bufferlen=          60000
maxbl3=               60000
maxbl4=               60000
  allocated                 295297  DP for d2temp 
sifcfg setup: record length 4096 DP
# d1 elements per record  3272
# d2 elements per record  2730
  The MR-CISD density will be calculated.
 item #                     1 suffix=:.drt1.state2:
 read_civout: repnuc=  -53.6745359128513     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 96% root-following 0
 MR-CISD energy:   -94.66619496   -40.99165905
 residuum:     0.00005971
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.34532033   -40.67078441
 residuum:     0.00006096
 deltae:     0.00000001
 apxde:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873504    -0.00000002    -0.01233979    -0.01513292    -0.01976758    -0.01841314    -0.26930823    -0.00282862
 ref:   2     0.00000054     0.94784829     0.00000041     0.00103390     0.00131339    -0.00807063     0.01443846    -0.06653525
 ref:   3    -0.00691241     0.00000231    -0.94558414     0.07466998    -0.00375523     0.12230225    -0.04499544    -0.23661551

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873504    -0.00000002    -0.01233979    -0.01513292    -0.01976758    -0.01841314     0.00000000     0.00000000
 ref:   2     0.00000054     0.94784829     0.00000041     0.00103390     0.00131339    -0.00807063     0.00000000     0.00000000
 ref:   3    -0.00691241     0.00000231    -0.94558414     0.07466998    -0.00375523     0.12230225     0.00000000     0.00000000

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
  DDZI=     189 DDYI=     614 DDXI=     452 DDWI=     599
  DDZE=       0 DDYE=      92 DDXE=      71 DDWE=      98
================================================================================
--------------------------------------------------------------------------------
  2e-density  for root #    2
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
 starting prpd24 .... 
================================================================================
   DYZ=     215  DYX=    4713  DYW=    5711
   D0Z=     147  D0Y=    7060  D0X=    3958  D0W=    5934
  DDZI=     735 DDYI=    2366 DDXI=    1673 DDWI=    2150
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

          modens reordered block   1

               a     1        a     2        a     3        a     4        a     5        a     6        a     7        a     8
  a     1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     3    0.00000        0.00000        1.98322       1.001876E-03   7.360936E-08   2.301955E-03   8.201598E-08   3.563811E-07
  a     4    0.00000        0.00000       1.001876E-03    1.97207       6.028850E-09   1.389971E-03  -3.347586E-07   1.072592E-07
  a     5    0.00000        0.00000       7.360936E-08   6.028850E-09    1.95895       1.615983E-07  -6.012370E-02  -5.407067E-06
  a     6    0.00000        0.00000       2.301955E-03   1.389971E-03   1.615983E-07    1.96547       1.576769E-06  -4.552090E-08
  a     7    0.00000        0.00000       8.201598E-08  -3.347586E-07  -6.012370E-02   1.576769E-06    1.00855      -8.182640E-05
  a     8    0.00000        0.00000       3.563811E-07   1.072592E-07  -5.407067E-06  -4.552090E-08  -8.182640E-05    1.90349    
  a     9    0.00000        0.00000      -2.051428E-06  -4.278993E-07  -8.670395E-06   1.591864E-08  -1.700745E-05   0.226204    
  a    10    0.00000        0.00000      -9.057358E-03  -6.687035E-03   4.915394E-07  -4.589857E-03  -7.368714E-07  -6.669239E-07
  a    11    0.00000        0.00000       4.243856E-07   5.252310E-07   7.160325E-03  -1.176652E-07  -1.681864E-02   1.197780E-05
  a    12    0.00000        0.00000      -4.225058E-03  -6.251882E-03   1.953232E-07   6.469213E-03  -5.444530E-07  -1.330042E-07
  a    13    0.00000        0.00000      -6.915123E-08   1.303801E-07   9.173621E-03  -1.918905E-07   1.868668E-02   4.549086E-06
  a    14    0.00000        0.00000      -7.016797E-03  -1.899063E-03   1.386824E-07   1.627915E-02  -3.347642E-07   5.174723E-07
  a    15    0.00000        0.00000       8.436127E-07  -2.560478E-06  -2.898127E-08  -2.306288E-06   1.801903E-06  -1.809910E-02
  a    16    0.00000        0.00000       2.894575E-03  -6.677543E-03   1.218521E-06  -5.877782E-03   3.015707E-06   6.935660E-06
  a    17    0.00000        0.00000      -6.523423E-07   1.486936E-06   5.746939E-03   1.104845E-06   1.424351E-02   3.341133E-06
  a    18    0.00000        0.00000       3.427344E-03  -3.856557E-05  -1.510948E-07  -8.977617E-03  -6.322132E-07  -1.375210E-07
  a    19    0.00000        0.00000       3.864885E-06   3.115154E-06   2.792857E-07  -5.178997E-06  -8.368844E-07  -2.463217E-03
  a    20    0.00000        0.00000      -4.487008E-03  -3.227228E-03   1.816544E-07   5.236425E-03  -2.784582E-07  -3.214683E-06
  a    21    0.00000        0.00000       3.081426E-08   7.926796E-08   1.319611E-04  -1.852413E-07  -1.926881E-02   3.841119E-07
  a    22    0.00000        0.00000      -1.503104E-03   3.274054E-03   2.892096E-08  -4.581192E-03   4.666304E-07  -3.552293E-07
  a    23    0.00000        0.00000      -4.237229E-08  -8.185056E-08   1.531451E-03   1.887098E-07   1.602062E-02   3.708900E-06
  a    24    0.00000        0.00000       7.173948E-04  -2.020106E-03   1.536986E-08   7.720475E-05   2.387290E-07  -3.042803E-07
  a    25    0.00000        0.00000      -9.433898E-05   2.766651E-03   1.708833E-07   2.018619E-03   1.041083E-06   6.715527E-07
  a    26    0.00000        0.00000       1.221448E-08  -9.662672E-08   2.659080E-03  -2.607963E-07   2.010233E-02  -1.960991E-06
  a    27    0.00000        0.00000       1.265161E-07   2.681840E-07   2.118340E-07   2.982245E-07   8.640479E-07   1.873230E-02
  a    28    0.00000        0.00000       3.443595E-07  -1.671164E-08  -1.071663E-06   1.282229E-07  -7.858905E-06  -3.027666E-07
  a    29    0.00000        0.00000      -2.028869E-10   6.685154E-08   8.418567E-07  -6.791708E-07  -1.894774E-06  -3.663844E-08
  a    30    0.00000        0.00000      -4.343969E-04  -1.506993E-03   1.175193E-07   7.242184E-03   1.367474E-10   1.596516E-07
  a    31    0.00000        0.00000      -2.724326E-07   2.667542E-07   6.015553E-03   5.580320E-07  -5.520002E-03   4.768023E-07
  a    32    0.00000        0.00000      -5.758358E-03   3.558111E-03  -1.681087E-07   1.379824E-02   2.413919E-07  -3.801645E-07
  a    33    0.00000        0.00000       7.378402E-04   4.502272E-03  -2.113336E-07  -3.752779E-03  -1.839658E-07   1.653821E-06
  a    34    0.00000        0.00000      -5.845171E-08  -9.121470E-07   2.977044E-08   5.547007E-07  -4.287184E-07   5.514457E-03
  a    35    0.00000        0.00000       3.617865E-05   1.988595E-03   2.850682E-07   5.037059E-04   6.333543E-07  -1.455443E-07
  a    36    0.00000        0.00000       5.216395E-09  -1.439573E-07   2.496475E-03  -1.137150E-07   4.458214E-03   7.320455E-08

               a     9        a    10        a    11        a    12        a    13        a    14        a    15        a    16
  a     3  -2.051428E-06  -9.057358E-03   4.243856E-07  -4.225058E-03  -6.915123E-08  -7.016797E-03   8.436127E-07   2.894575E-03
  a     4  -4.278993E-07  -6.687035E-03   5.252310E-07  -6.251882E-03   1.303801E-07  -1.899063E-03  -2.560478E-06  -6.677543E-03
  a     5  -8.670395E-06   4.915394E-07   7.160325E-03   1.953232E-07   9.173621E-03   1.386824E-07  -2.898127E-08   1.218521E-06
  a     6   1.591864E-08  -4.589857E-03  -1.176652E-07   6.469213E-03  -1.918905E-07   1.627915E-02  -2.306288E-06  -5.877782E-03
  a     7  -1.700745E-05  -7.368714E-07  -1.681864E-02  -5.444530E-07   1.868668E-02  -3.347642E-07   1.801903E-06   3.015707E-06
  a     8   0.226204      -6.669239E-07   1.197780E-05  -1.330042E-07   4.549086E-06   5.174723E-07  -1.809910E-02   6.935660E-06
  a     9    1.06026       2.639038E-07  -3.801375E-05   7.989356E-07  -3.205950E-06  -8.281276E-07  -1.438052E-02   4.235217E-06
  a    10   2.639038E-07   9.105778E-03   2.459354E-08  -3.220415E-04   2.338466E-08   2.797774E-05   7.611870E-07   2.217224E-03
  a    11  -3.801375E-05   2.459354E-08   9.345733E-03  -3.527161E-08   4.727813E-04  -4.252092E-08  -2.191746E-07  -8.949706E-07
  a    12   7.989356E-07  -3.220415E-04  -3.527161E-08   1.008296E-02   4.218184E-08   9.212211E-04   3.012715E-07   1.026344E-03
  a    13  -3.205950E-06   2.338466E-08   4.727813E-04   4.218184E-08   5.812918E-03  -8.818535E-08  -3.197502E-07   5.957603E-09
  a    14  -8.281276E-07   2.797774E-05  -4.252092E-08   9.212211E-04  -8.818535E-08   7.972288E-03   1.046640E-06   3.324901E-03
  a    15  -1.438052E-02   7.611870E-07  -2.191746E-07   3.012715E-07  -3.197502E-07   1.046640E-06   8.311532E-03  -1.251765E-10
  a    16   4.235217E-06   2.217224E-03  -8.949706E-07   1.026344E-03   5.957603E-09   3.324901E-03  -1.251765E-10   8.356981E-03
  a    17  -3.087325E-06  -5.162652E-07  -3.160805E-03  -4.189815E-07   2.808543E-04  -7.644525E-07  -1.266925E-08  -5.138119E-07
  a    18  -6.884534E-07   2.117570E-03   1.581122E-07  -3.501154E-03  -8.881395E-09  -7.419428E-04   9.898691E-09  -5.996578E-05
  a    19   2.952831E-02  -3.127931E-06  -5.018436E-07  -3.998405E-06  -1.940684E-07  -3.629251E-06   8.658820E-04  -3.658326E-06
  a    20   2.824615E-05   3.042088E-03  -2.191141E-07   4.094035E-03  -5.176338E-08   3.683625E-03   1.966970E-06   3.457017E-03
  a    21  -6.240647E-06   2.578981E-07   5.688619E-03   1.208236E-07  -1.839394E-03  -1.972110E-08   1.690897E-07  -5.954942E-07
  a    22   2.531060E-07   5.807111E-03   2.959292E-08  -4.714931E-03   1.552354E-07  -1.982310E-03  -6.841582E-08  -3.880739E-04
  a    23  -4.955306E-06  -2.269688E-08   4.104288E-03   2.102635E-07   5.240245E-03   1.337223E-07  -2.317366E-07  -5.154592E-08
  a    24   8.060159E-07   2.076684E-03  -1.934237E-07   5.101605E-03   1.613412E-07  -4.930913E-03  -5.259434E-07  -1.928187E-03
  a    25  -1.894799E-07   2.993840E-03  -2.172940E-07   1.625315E-03   1.781366E-08   2.048024E-03   2.922274E-07   9.970260E-04
  a    26   1.815754E-06  -1.407833E-07  -6.577454E-04  -6.457007E-08   1.870391E-04  -8.018475E-08   1.793265E-08  -1.116797E-07
  a    27  -4.427658E-03  -6.668804E-08   4.953612E-07  -4.598684E-08   2.765707E-07   1.696034E-09  -6.262657E-04   2.117947E-07
  a    28   2.354981E-08  -1.317271E-07   2.689256E-07  -1.020486E-07   1.103582E-07  -1.072440E-08   2.544780E-09  -3.188246E-08
  a    29   1.951310E-08  -1.374160E-07   4.275076E-08   1.321118E-07   1.066787E-07  -2.447490E-09  -1.688643E-09   2.071327E-08
  a    30   1.687189E-07   6.517179E-04   3.363524E-09  -6.005674E-04  -2.279794E-09   1.208991E-05  -8.296321E-08  -1.919438E-04
  a    31  -1.212339E-06  -1.413156E-08   1.400963E-04   1.361774E-08   3.868586E-04   1.122954E-08   2.490542E-08   2.794920E-07
  a    32   1.016270E-07  -3.080168E-04   2.047862E-08  -1.835554E-04  -1.533277E-08   3.211675E-04  -7.354573E-07  -2.240983E-03
  a    33  -6.183713E-07   2.268072E-04  -4.023982E-08   8.274772E-04  -1.260728E-08   6.413498E-05   1.775348E-07   5.490909E-04
  a    34  -4.894806E-03  -6.793316E-08   2.337267E-07  -1.688017E-07   5.554175E-08  -8.704472E-09  -1.696263E-04   1.963322E-08
  a    35  -2.951340E-07  -2.808616E-04  -1.561885E-08  -5.492032E-05   7.037668E-08   3.444557E-04  -8.219528E-08  -1.494140E-04
  a    36   8.806983E-08   2.557915E-08  -2.298680E-04  -1.031307E-08   5.659293E-04  -5.010848E-08   1.537108E-08   9.504195E-08

               a    17        a    18        a    19        a    20        a    21        a    22        a    23        a    24
  a     3  -6.523423E-07   3.427344E-03   3.864885E-06  -4.487008E-03   3.081426E-08  -1.503104E-03  -4.237229E-08   7.173948E-04
  a     4   1.486936E-06  -3.856557E-05   3.115154E-06  -3.227228E-03   7.926796E-08   3.274054E-03  -8.185056E-08  -2.020106E-03
  a     5   5.746939E-03  -1.510948E-07   2.792857E-07   1.816544E-07   1.319611E-04   2.892096E-08   1.531451E-03   1.536986E-08
  a     6   1.104845E-06  -8.977617E-03  -5.178997E-06   5.236425E-03  -1.852413E-07  -4.581192E-03   1.887098E-07   7.720475E-05
  a     7   1.424351E-02  -6.322132E-07  -8.368844E-07  -2.784582E-07  -1.926881E-02   4.666304E-07   1.602062E-02   2.387290E-07
  a     8   3.341133E-06  -1.375210E-07  -2.463217E-03  -3.214683E-06   3.841119E-07  -3.552293E-07   3.708900E-06  -3.042803E-07
  a     9  -3.087325E-06  -6.884534E-07   2.952831E-02   2.824615E-05  -6.240647E-06   2.531060E-07  -4.955306E-06   8.060159E-07
  a    10  -5.162652E-07   2.117570E-03  -3.127931E-06   3.042088E-03   2.578981E-07   5.807111E-03  -2.269688E-08   2.076684E-03
  a    11  -3.160805E-03   1.581122E-07  -5.018436E-07  -2.191141E-07   5.688619E-03   2.959292E-08   4.104288E-03  -1.934237E-07
  a    12  -4.189815E-07  -3.501154E-03  -3.998405E-06   4.094035E-03   1.208236E-07  -4.714931E-03   2.102635E-07   5.101605E-03
  a    13   2.808543E-04  -8.881395E-09  -1.940684E-07  -5.176338E-08  -1.839394E-03   1.552354E-07   5.240245E-03   1.613412E-07
  a    14  -7.644525E-07  -7.419428E-04  -3.629251E-06   3.683625E-03  -1.972110E-08  -1.982310E-03   1.337223E-07  -4.930913E-03
  a    15  -1.266925E-08   9.898691E-09   8.658820E-04   1.966970E-06   1.690897E-07  -6.841582E-08  -2.317366E-07  -5.259434E-07
  a    16  -5.138119E-07  -5.996578E-05  -3.658326E-06   3.457017E-03  -5.954942E-07  -3.880739E-04  -5.154592E-08  -1.928187E-03
  a    17   5.914550E-03   5.584198E-08  -1.198289E-07  -8.173760E-07  -2.541171E-03   1.913686E-07  -4.549470E-04   3.964115E-07
  a    18   5.584198E-08   7.502952E-03   8.916927E-07  -1.058682E-03   1.360022E-07   3.206469E-03  -6.522408E-08  -9.046575E-04
  a    19  -1.198289E-07   8.916927E-07   1.048433E-02   2.001469E-06  -2.733056E-07   6.759761E-07  -3.792769E-07  -8.549364E-08
  a    20  -8.173760E-07  -1.058682E-03   2.001469E-06   8.204491E-03  -1.333413E-08  -6.663739E-04   5.588621E-08   1.686720E-04
  a    21  -2.541171E-03   1.360022E-07  -2.733056E-07  -1.333413E-08   4.834974E-03   5.414318E-08   4.431213E-04  -1.830068E-09
  a    22   1.913686E-07   3.206469E-03   6.759761E-07  -6.663739E-04   5.414318E-08   6.801850E-03  -5.105710E-09   4.374343E-04
  a    23  -4.549470E-04  -6.522408E-08  -3.792769E-07   5.588621E-08   4.431213E-04  -5.105710E-09   7.526429E-03  -2.590731E-08
  a    24   3.964115E-07  -9.046575E-04  -8.549364E-08   1.686720E-04  -1.830068E-09   4.374343E-04  -2.590731E-08   7.711132E-03
  a    25  -2.253346E-07  -1.149224E-04  -1.296423E-06   1.316059E-03  -3.735370E-08   7.351090E-04   2.031357E-08   2.792340E-05
  a    26  -5.175630E-04   1.173355E-08   2.445278E-07  -4.176936E-09  -6.376752E-04  -2.387215E-08   5.633413E-04  -7.391846E-09
  a    27   3.491860E-08   4.189854E-08  -1.259646E-03  -1.317543E-06  -3.718985E-08  -7.639506E-08   2.922836E-07  -6.326980E-08
  a    28   2.962132E-07   6.266119E-08   7.518841E-09  -5.407145E-08   2.507350E-07   4.183691E-08  -1.760875E-07  -6.999454E-08
  a    29   2.784649E-07  -2.743171E-07  -3.028952E-10  -8.944839E-08  -2.420201E-08  -1.081887E-07   1.057068E-07   3.798260E-08
  a    30   1.118057E-07   1.940576E-03  -6.732209E-07   8.514852E-04  -7.861140E-09   4.324575E-04  -4.255420E-09  -2.367003E-04
  a    31   1.600522E-03  -7.212442E-08  -8.748007E-08   5.023967E-09  -3.051492E-04   1.486955E-08   7.010570E-04   1.007359E-08
  a    32   4.358567E-07  -2.999118E-04   2.549766E-07  -1.966394E-04   1.168977E-08   1.767975E-04  -4.052002E-08   8.603679E-05
  a    33  -1.573798E-07  -9.208896E-04  -1.495302E-06   1.808306E-03  -1.169660E-08  -2.403481E-04   1.618228E-09   3.151135E-04
  a    34   6.734189E-08   1.270895E-07   6.723581E-04   2.201192E-07  -7.205150E-09   5.433721E-09   4.141398E-08  -9.861959E-08
  a    35   6.047396E-08  -3.688619E-04  -7.673086E-07   7.188131E-04  -5.699913E-08  -4.151760E-04   4.440514E-08  -8.938542E-04
  a    36   3.079858E-04   3.663917E-08   1.065230E-08  -8.339118E-08  -3.951458E-04   6.257631E-08   1.119177E-04   1.212323E-07

               a    25        a    26        a    27        a    28        a    29        a    30        a    31        a    32
  a     3  -9.433898E-05   1.221448E-08   1.265161E-07   3.443595E-07  -2.028869E-10  -4.343969E-04  -2.724326E-07  -5.758358E-03
  a     4   2.766651E-03  -9.662672E-08   2.681840E-07  -1.671164E-08   6.685154E-08  -1.506993E-03   2.667542E-07   3.558111E-03
  a     5   1.708833E-07   2.659080E-03   2.118340E-07  -1.071663E-06   8.418567E-07   1.175193E-07   6.015553E-03  -1.681087E-07
  a     6   2.018619E-03  -2.607963E-07   2.982245E-07   1.282229E-07  -6.791708E-07   7.242184E-03   5.580320E-07   1.379824E-02
  a     7   1.041083E-06   2.010233E-02   8.640479E-07  -7.858905E-06  -1.894774E-06   1.367474E-10  -5.520002E-03   2.413919E-07
  a     8   6.715527E-07  -1.960991E-06   1.873230E-02  -3.027666E-07  -3.663844E-08   1.596516E-07   4.768023E-07  -3.801645E-07
  a     9  -1.894799E-07   1.815754E-06  -4.427658E-03   2.354981E-08   1.951310E-08   1.687189E-07  -1.212339E-06   1.016270E-07
  a    10   2.993840E-03  -1.407833E-07  -6.668804E-08  -1.317271E-07  -1.374160E-07   6.517179E-04  -1.413156E-08  -3.080168E-04
  a    11  -2.172940E-07  -6.577454E-04   4.953612E-07   2.689256E-07   4.275076E-08   3.363524E-09   1.400963E-04   2.047862E-08
  a    12   1.625315E-03  -6.457007E-08  -4.598684E-08  -1.020486E-07   1.321118E-07  -6.005674E-04   1.361774E-08  -1.835554E-04
  a    13   1.781366E-08   1.870391E-04   2.765707E-07   1.103582E-07   1.066787E-07  -2.279794E-09   3.868586E-04  -1.533277E-08
  a    14   2.048024E-03  -8.018475E-08   1.696034E-09  -1.072440E-08  -2.447490E-09   1.208991E-05   1.122954E-08   3.211675E-04
  a    15   2.922274E-07   1.793265E-08  -6.262657E-04   2.544780E-09  -1.688643E-09  -8.296321E-08   2.490542E-08  -7.354573E-07
  a    16   9.970260E-04  -1.116797E-07   2.117947E-07  -3.188246E-08   2.071327E-08  -1.919438E-04   2.794920E-07  -2.240983E-03
  a    17  -2.253346E-07  -5.175630E-04   3.491860E-08   2.962132E-07   2.784649E-07   1.118057E-07   1.600522E-03   4.358567E-07
  a    18  -1.149224E-04   1.173355E-08   4.189854E-08   6.266119E-08  -2.743171E-07   1.940576E-03  -7.212442E-08  -2.999118E-04
  a    19  -1.296423E-06   2.445278E-07  -1.259646E-03   7.518841E-09  -3.028952E-10  -6.732209E-07  -8.748007E-08   2.549766E-07
  a    20   1.316059E-03  -4.176936E-09  -1.317543E-06  -5.407145E-08  -8.944839E-08   8.514852E-04   5.023967E-09  -1.966394E-04
  a    21  -3.735370E-08  -6.376752E-04  -3.718985E-08   2.507350E-07  -2.420201E-08  -7.861140E-09  -3.051492E-04   1.168977E-08
  a    22   7.351090E-04  -2.387215E-08  -7.639506E-08   4.183691E-08  -1.081887E-07   4.324575E-04   1.486955E-08   1.767975E-04
  a    23   2.031357E-08   5.633413E-04   2.922836E-07  -1.760875E-07   1.057068E-07  -4.255420E-09   7.010570E-04  -4.052002E-08
  a    24   2.792340E-05  -7.391846E-09  -6.326980E-08  -6.999454E-08   3.798260E-08  -2.367003E-04   1.007359E-08   8.603679E-05
  a    25   3.389070E-03  -6.398046E-08   1.359708E-07  -3.878869E-08   2.062993E-08  -3.068131E-04  -3.624406E-08  -2.034814E-04
  a    26  -6.398046E-08   1.776165E-03  -6.513364E-07   1.178053E-06   3.102322E-07   1.998323E-08  -5.753551E-04   1.615825E-08
  a    27   1.359708E-07  -6.513364E-07   6.221927E-03  -1.755330E-08   6.679415E-09  -8.115943E-09  -2.449592E-08   7.811887E-11
  a    28  -3.878869E-08   1.178053E-06  -1.755330E-08   4.150967E-03   7.134999E-04   8.322264E-08   1.130029E-07  -7.722352E-09
  a    29   2.062993E-08   3.102322E-07   6.679415E-09   7.134999E-04   2.750709E-03  -1.027497E-09  -3.760842E-08   6.136179E-08
  a    30  -3.068131E-04   1.998323E-08  -8.115943E-09   8.322264E-08  -1.027497E-09   2.257841E-03   1.589002E-09  -1.520047E-04
  a    31  -3.624406E-08  -5.753551E-04  -2.449592E-08   1.130029E-07  -3.760842E-08   1.589002E-09   1.881861E-03   1.720660E-08
  a    32  -2.034814E-04   1.615825E-08   7.811887E-11  -7.722352E-09   6.136179E-08  -1.520047E-04   1.720660E-08   2.171105E-03
  a    33  -4.401466E-04   5.158637E-08   1.024903E-07  -2.811296E-08  -2.810541E-08   3.507416E-04  -4.747803E-10  -4.439798E-05
  a    34   1.293703E-07  -7.668506E-08   5.441141E-04  -5.852545E-09   1.114204E-09  -7.398849E-08   1.935767E-08  -1.362398E-09
  a    35   6.458634E-04  -9.430749E-08  -1.433336E-08   3.165511E-08   1.543402E-08  -6.062022E-05   3.725694E-08   9.884915E-05
  a    36  -1.200511E-07  -5.137202E-04  -5.436526E-08   1.928689E-07   5.287566E-08   1.331314E-08   2.250511E-04  -2.386723E-08

               a    33        a    34        a    35        a    36
  a     3   7.378402E-04  -5.845171E-08   3.617865E-05   5.216395E-09
  a     4   4.502272E-03  -9.121470E-07   1.988595E-03  -1.439573E-07
  a     5  -2.113336E-07   2.977044E-08   2.850682E-07   2.496475E-03
  a     6  -3.752779E-03   5.547007E-07   5.037059E-04  -1.137150E-07
  a     7  -1.839658E-07  -4.287184E-07   6.333543E-07   4.458214E-03
  a     8   1.653821E-06   5.514457E-03  -1.455443E-07   7.320455E-08
  a     9  -6.183713E-07  -4.894806E-03  -2.951340E-07   8.806983E-08
  a    10   2.268072E-04  -6.793316E-08  -2.808616E-04   2.557915E-08
  a    11  -4.023982E-08   2.337267E-07  -1.561885E-08  -2.298680E-04
  a    12   8.274772E-04  -1.688017E-07  -5.492032E-05  -1.031307E-08
  a    13  -1.260728E-08   5.554175E-08   7.037668E-08   5.659293E-04
  a    14   6.413498E-05  -8.704472E-09   3.444557E-04  -5.010848E-08
  a    15   1.775348E-07  -1.696263E-04  -8.219528E-08   1.537108E-08
  a    16   5.490909E-04   1.963322E-08  -1.494140E-04   9.504195E-08
  a    17  -1.573798E-07   6.734189E-08   6.047396E-08   3.079858E-04
  a    18  -9.208896E-04   1.270895E-07  -3.688619E-04   3.663917E-08
  a    19  -1.495302E-06   6.723581E-04  -7.673086E-07   1.065230E-08
  a    20   1.808306E-03   2.201192E-07   7.188131E-04  -8.339118E-08
  a    21  -1.169660E-08  -7.205150E-09  -5.699913E-08  -3.951458E-04
  a    22  -2.403481E-04   5.433721E-09  -4.151760E-04   6.257631E-08
  a    23   1.618228E-09   4.141398E-08   4.440514E-08   1.119177E-04
  a    24   3.151135E-04  -9.861959E-08  -8.938542E-04   1.212323E-07
  a    25  -4.401466E-04   1.293703E-07   6.458634E-04  -1.200511E-07
  a    26   5.158637E-08  -7.668506E-08  -9.430749E-08  -5.137202E-04
  a    27   1.024903E-07   5.441141E-04  -1.433336E-08  -5.436526E-08
  a    28  -2.811296E-08  -5.852545E-09   3.165511E-08   1.928689E-07
  a    29  -2.810541E-08   1.114204E-09   1.543402E-08   5.287566E-08
  a    30   3.507416E-04  -7.398849E-08  -6.062022E-05   1.331314E-08
  a    31  -4.747803E-10   1.935767E-08   3.725694E-08   2.250511E-04
  a    32  -4.439798E-05  -1.362398E-09   9.884915E-05  -2.386723E-08
  a    33   1.529850E-03   1.082293E-08   9.576031E-05  -1.107972E-08
  a    34   1.082293E-08   1.585710E-03  -2.386643E-08   7.004229E-09
  a    35   9.576031E-05  -2.386643E-08   1.310794E-03  -3.690165E-08
  a    36  -1.107972E-08   7.004229E-09  -3.690165E-08   1.011568E-03

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98372233     1.97227684     1.96535163     1.96284386     1.96073568     1.00669747
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     1.00447814     0.02025159     0.01786924     0.01652254     0.01574186     0.01057814     0.01048959     0.00740995
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00602144     0.00577735     0.00569166     0.00501611     0.00453974     0.00445047     0.00245121     0.00199973
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00185050     0.00143410     0.00137402     0.00105104     0.00088078     0.00066026     0.00053604     0.00033777
              MO    33       MO    34       MO    35       MO    36
  occ(*)=     0.00032560     0.00025571     0.00019385     0.00018377


 total number of electrons =   16.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
    N1_ s       1.998472  -0.000030   1.441541   0.083474   0.037885   0.000000
    N1_ p       0.000000  -0.000292  -0.002595   0.720367   0.434739   1.209708
    N1_ d       0.000000  -0.000105   0.004476   0.005717   0.001588   0.002347
    C1_ s       0.000539   1.999452   0.271769   0.782629   0.204490   0.000000
    C1_ p       0.000476   0.000003   0.040867   0.012153   0.793641   0.212500
    C1_ d      -0.000081   0.000001   0.006621   0.002986   0.012798   0.002838
    H1_ s       0.000301   0.000003   0.096056   0.101878   0.027096   0.231350
    H2_ s       0.000301   0.000003   0.096056   0.101871   0.027093   0.231362
    H3_ s      -0.000004   0.000483   0.014465   0.080601   0.213014   0.036368
    H4_ s      -0.000004   0.000483   0.014465   0.080601   0.213008   0.036370
 
   ao class       7a         8a         9a        10a        11a        12a  
    N1_ s       0.000000   0.000000   0.000000   0.001656   0.003226   0.000000
    N1_ p       1.265493   0.068605   0.339661   0.002310   0.000433   0.006629
    N1_ d       0.004347   0.006436   0.006058   0.000402   0.000872   0.001469
    C1_ s       0.000000   0.000000   0.000000   0.006924   0.000020   0.000000
    C1_ p       0.681684   0.450258   0.652941   0.000536   0.006857   0.000135
    C1_ d       0.009211   0.009822   0.005817   0.001959   0.001163   0.000686
    H1_ s       0.000000   0.034488   0.000000   0.000094   0.000963   0.003360
    H2_ s       0.000000   0.034489   0.000000   0.000094   0.000963   0.003361
    H3_ s       0.000000   0.201298   0.000001   0.003138   0.001686   0.000442
    H4_ s       0.000000   0.201302   0.000001   0.003138   0.001686   0.000442
 
   ao class      13a        14a        15a        16a        17a        18a  
    N1_ s       0.001845   0.000000   0.000000   0.000000   0.000516   0.000000
    N1_ p       0.005466   0.010131   0.002733   0.000080   0.001142   0.000114
    N1_ d       0.000714   0.000005   0.001108   0.000004   0.002332   0.003685
    C1_ s       0.000022   0.000000   0.000000   0.000000   0.000345   0.000000
    C1_ p       0.001357   0.000149   0.002782   0.007229   0.000755   0.000013
    C1_ d       0.000266   0.000293   0.000948   0.000096   0.000940   0.001966
    H1_ s       0.002768   0.000000   0.000065   0.000000  -0.000001  -0.000000
    H2_ s       0.002768   0.000000   0.000065   0.000000  -0.000001   0.000000
    H3_ s       0.000268   0.000000   0.001394   0.000000  -0.000003  -0.000000
    H4_ s       0.000269   0.000000   0.001394   0.000000  -0.000003   0.000000
 
   ao class      19a        20a        21a        22a        23a        24a  
    N1_ s       0.000225  -0.000000   0.000743   0.000000   0.000000   0.000244
    N1_ p       0.000219   0.000635   0.000418   0.000000   0.000000   0.000259
    N1_ d       0.001180   0.000761   0.001924   0.003795   0.000361   0.000237
    C1_ s       0.002195   0.000000   0.000601   0.000000  -0.000000  -0.000013
    C1_ p       0.000686   0.002543   0.000570   0.000000   0.000000   0.000205
    C1_ d       0.001125   0.000540   0.000293   0.000656   0.002090   0.001052
    H1_ s       0.000051   0.000018  -0.000016  -0.000000   0.000000   0.000005
    H2_ s       0.000051   0.000018  -0.000016  -0.000000   0.000000   0.000005
    H3_ s      -0.000020   0.000251   0.000011  -0.000000  -0.000000   0.000003
    H4_ s      -0.000020   0.000251   0.000011   0.000000   0.000000   0.000003
 
   ao class      25a        26a        27a        28a        29a        30a  
    N1_ s      -0.000000   0.000288   0.000000   0.000000  -0.000002   0.000070
    N1_ p       0.000062   0.000016   0.000039   0.000198   0.000266  -0.000006
    N1_ d       0.000713   0.000435   0.000485   0.000141   0.000241   0.000062
    C1_ s      -0.000000   0.000107   0.000000  -0.000000   0.000035   0.000032
    C1_ p       0.000286   0.000107   0.000014   0.000098   0.000057  -0.000013
    C1_ d       0.000403   0.000131   0.000836   0.000279   0.000052   0.000252
    H1_ s       0.000046   0.000129   0.000000   0.000143   0.000099   0.000019
    H2_ s       0.000046   0.000129   0.000000   0.000143   0.000099   0.000019
    H3_ s       0.000147   0.000045  -0.000000   0.000025   0.000017   0.000113
    H4_ s       0.000147   0.000045   0.000000   0.000025   0.000017   0.000113
 
   ao class      31a        32a        33a        34a        35a        36a  
    N1_ s      -0.000000   0.000000   0.000016   0.000000   0.000012   0.000001
    N1_ p       0.000019   0.000035   0.000024   0.000043   0.000006   0.000012
    N1_ d       0.000064  -0.000000   0.000006   0.000002   0.000004   0.000007
    C1_ s       0.000000   0.000000   0.000016   0.000000   0.000038   0.000007
    C1_ p      -0.000001   0.000005   0.000005   0.000103   0.000007   0.000079
    C1_ d       0.000128   0.000006   0.000015   0.000002   0.000016   0.000038
    H1_ s       0.000034   0.000126   0.000107   0.000009   0.000012   0.000004
    H2_ s       0.000034   0.000126   0.000107   0.000009   0.000012   0.000004
    H3_ s       0.000129   0.000019   0.000015   0.000044   0.000043   0.000016
    H4_ s       0.000129   0.000019   0.000015   0.000044   0.000043   0.000016


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.570180   3.269210   0.499208   0.499210   0.554008   0.554008
      p         4.066972   2.869088   0.000000   0.000000   0.000000   0.000000
      d         0.051870   0.066246   0.000000   0.000000   0.000000   0.000000
    total       7.689021   6.204544   0.499208   0.499210   0.554008   0.554008
 

 Total number of electrons:   16.00000000

 item #                     2 suffix=:.drt1.state3:
 read_civout: repnuc=  -53.6745359128513     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 96% root-following 0
 MR-CISD energy:   -94.66619496   -40.99165905
 residuum:     0.00005971
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.34532033   -40.67078441
 residuum:     0.00006096
 deltae:     0.00000001
 apxde:     0.00000000
================================================================================
  Reading record                      3  of civout
 INFO:ref#  3vector#  3 method:  0 last record  1max overlap with ref# 95% root-following 0
 MR-CISD energy:   -94.30933709   -40.63480117
 residuum:     0.00002356
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873504    -0.00000002    -0.01233979    -0.01513292    -0.01976758    -0.01841314    -0.26930823    -0.00282862
 ref:   2     0.00000054     0.94784829     0.00000041     0.00103390     0.00131339    -0.00807063     0.01443846    -0.06653525
 ref:   3    -0.00691241     0.00000231    -0.94558414     0.07466998    -0.00375523     0.12230225    -0.04499544    -0.23661551

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.95873504    -0.00000002    -0.01233979    -0.01513292    -0.01976758    -0.01841314     0.00000000     0.00000000
 ref:   2     0.00000054     0.94784829     0.00000041     0.00103390     0.00131339    -0.00807063     0.00000000     0.00000000
 ref:   3    -0.00691241     0.00000231    -0.94558414     0.07466998    -0.00375523     0.12230225     0.00000000     0.00000000

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
  DDZI=     189 DDYI=     614 DDXI=     452 DDWI=     599
  DDZE=       0 DDYE=      92 DDXE=      71 DDWE=      98
================================================================================
--------------------------------------------------------------------------------
  2e-density  for root #    3
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
 starting prpd24 .... 
================================================================================
   DYZ=     215  DYX=    4713  DYW=    5711
   D0Z=     147  D0Y=    7060  D0X=    3958  D0W=    5934
  DDZI=     735 DDYI=    2366 DDXI=    1673 DDWI=    2150
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
Trace of MO density:    12.000000
   12  correlated and     4  frozen core electrons

          modens reordered block   1

               a     1        a     2        a     3        a     4        a     5        a     6        a     7        a     8
  a     1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  a     3    0.00000        0.00000        1.98000       1.975906E-03   3.403934E-08   5.252655E-03   4.979930E-08  -1.755362E-06
  a     4    0.00000        0.00000       1.975906E-03    1.97264       2.132504E-08  -2.996247E-03   4.317438E-08  -1.882925E-06
  a     5    0.00000        0.00000       3.403934E-08   2.132504E-08    1.97124      -2.341826E-08  -1.663264E-03   8.574344E-06
  a     6    0.00000        0.00000       5.252655E-03  -2.996247E-03  -2.341826E-08    1.96161      -3.812053E-08   4.506332E-06
  a     7    0.00000        0.00000       4.979930E-08   4.317438E-08  -1.663264E-03  -3.812053E-08    1.96988       7.546046E-05
  a     8    0.00000        0.00000      -1.755362E-06  -1.882925E-06   8.574344E-06   4.506332E-06   7.546046E-05    1.10389    
  a     9    0.00000        0.00000       9.120473E-07   1.508928E-06   1.029676E-05  -1.296819E-06   2.059630E-05  -0.478186    
  a    10    0.00000        0.00000      -5.640247E-03   6.815638E-03  -7.953412E-07  -6.119494E-03   3.829854E-07   1.648989E-06
  a    11    0.00000        0.00000       2.160978E-07  -4.440128E-07  -1.745113E-02   5.291075E-07   1.729070E-02   1.235026E-06
  a    12    0.00000        0.00000       1.562301E-03   4.535229E-03  -7.845595E-07  -4.727652E-03   4.465796E-07   6.179693E-07
  a    13    0.00000        0.00000       2.041570E-08  -6.614140E-08   2.815377E-03  -2.306496E-08  -5.761102E-03   5.261298E-07
  a    14    0.00000        0.00000       2.894210E-03   1.421899E-02  -1.982925E-07  -1.260011E-02   2.680655E-07  -7.495976E-08
  a    15    0.00000        0.00000      -7.279585E-07   3.641829E-06  -5.884597E-07   1.208895E-06  -2.561794E-06  -5.697554E-05
  a    16    0.00000        0.00000      -2.385303E-03   1.157387E-02   1.193198E-06   4.987646E-03  -3.069569E-06  -6.647012E-08
  a    17    0.00000        0.00000       3.930772E-07  -2.529592E-06   4.956794E-03  -9.784082E-07  -1.388070E-02  -1.728482E-06
  a    18    0.00000        0.00000      -4.751477E-03  -2.061632E-03  -3.575282E-08   9.448834E-03   3.736964E-07   7.363552E-07
  a    19    0.00000        0.00000      -1.376785E-07  -1.121599E-05   5.208950E-07   4.642353E-06  -1.714163E-07  -1.906014E-02
  a    20    0.00000        0.00000      -2.414058E-04   1.115197E-02  -1.770150E-07  -4.987385E-03   2.846261E-07  -1.692000E-05
  a    21    0.00000        0.00000       5.240491E-08  -1.608559E-07  -1.403161E-02   3.600186E-07   1.751895E-02   1.082708E-06
  a    22    0.00000        0.00000       1.383223E-03  -4.546034E-03  -2.077605E-08   8.909777E-03  -4.653866E-07   3.983092E-07
  a    23    0.00000        0.00000       2.860067E-08   2.023005E-07  -1.007396E-02  -1.956371E-07  -1.245258E-02  -2.835219E-07
  a    24    0.00000        0.00000      -4.639720E-05   3.053211E-03  -2.042666E-07   2.136889E-03  -2.550976E-07   4.705791E-07
  a    25    0.00000        0.00000       3.528719E-03   3.666105E-03   1.408136E-07  -2.769624E-03  -2.157042E-07  -4.788251E-07
  a    26    0.00000        0.00000      -1.843263E-07  -1.603277E-08   2.811997E-03   1.221747E-07  -3.079456E-03   7.346391E-07
  a    27    0.00000        0.00000      -3.615193E-07  -7.761118E-08   3.381289E-07   8.593768E-07   1.324532E-06  -4.350580E-03
  a    28    0.00000        0.00000      -1.080503E-07  -1.468037E-08  -5.700422E-07   5.733249E-08   6.192186E-07   5.446415E-08
  a    29    0.00000        0.00000      -1.226259E-07   1.121910E-07   7.488657E-07  -1.481404E-07  -1.041773E-07   1.681084E-09
  a    30    0.00000        0.00000       1.243839E-03  -8.314305E-04   1.505350E-08   1.979771E-03  -4.781091E-08   2.076626E-07
  a    31    0.00000        0.00000      -1.018073E-07   2.215458E-07   4.435024E-03   9.898722E-10  -2.894703E-03   1.129713E-08
  a    32    0.00000        0.00000      -1.974032E-03   2.751802E-03  -2.512876E-07   7.102997E-04   1.318997E-07  -4.386571E-08
  a    33    0.00000        0.00000       1.493906E-03   4.469467E-03  -1.492860E-07  -5.835565E-04  -4.524771E-08  -1.089612E-06
  a    34    0.00000        0.00000      -4.127698E-07  -1.120108E-06   4.567027E-07  -2.096555E-07   5.150717E-07  -6.404741E-03
  a    35    0.00000        0.00000      -6.985103E-04  -1.141256E-03   4.670115E-07   1.301790E-03   7.089996E-07   5.780877E-08
  a    36    0.00000        0.00000       9.108913E-08   2.456921E-07   3.579221E-03  -1.939498E-07   5.216734E-03   1.831761E-07

               a     9        a    10        a    11        a    12        a    13        a    14        a    15        a    16
  a     3   9.120473E-07  -5.640247E-03   2.160978E-07   1.562301E-03   2.041570E-08   2.894210E-03  -7.279585E-07  -2.385303E-03
  a     4   1.508928E-06   6.815638E-03  -4.440128E-07   4.535229E-03  -6.614140E-08   1.421899E-02   3.641829E-06   1.157387E-02
  a     5   1.029676E-05  -7.953412E-07  -1.745113E-02  -7.845595E-07   2.815377E-03  -1.982925E-07  -5.884597E-07   1.193198E-06
  a     6  -1.296819E-06  -6.119494E-03   5.291075E-07  -4.727652E-03  -2.306496E-08  -1.260011E-02   1.208895E-06   4.987646E-03
  a     7   2.059630E-05   3.829854E-07   1.729070E-02   4.465796E-07  -5.761102E-03   2.680655E-07  -2.561794E-06  -3.069569E-06
  a     8  -0.478186       1.648989E-06   1.235026E-06   6.179693E-07   5.261298E-07  -7.495976E-08  -5.697554E-05  -6.647012E-08
  a     9   0.881198       2.332921E-06   4.530439E-06  -6.552188E-08  -1.587943E-06   5.089790E-07  -4.654622E-03   1.849916E-06
  a    10   2.332921E-06   1.120764E-02  -3.664779E-09  -4.617319E-04  -6.694452E-09   7.153029E-04   1.026449E-06   3.143115E-03
  a    11   4.530439E-06  -3.664779E-09   1.025076E-02   1.422748E-08   1.328688E-04  -8.327064E-08  -6.833091E-09  -6.928899E-07
  a    12  -6.552188E-08  -4.617319E-04   1.422748E-08   9.635328E-03  -1.984675E-08   1.041633E-03   3.189248E-07   1.067252E-03
  a    13  -1.587943E-06  -6.694452E-09   1.328688E-04  -1.984675E-08   9.129343E-03  -1.474342E-07   7.698667E-08  -3.188931E-07
  a    14   5.089790E-07   7.153029E-04  -8.327064E-08   1.041633E-03  -1.474342E-07   1.172030E-02   2.052783E-06   6.964863E-03
  a    15  -4.654622E-03   1.026449E-06  -6.833091E-09   3.189248E-07   7.698667E-08   2.052783E-06   3.725603E-03   2.129519E-06
  a    16   1.849916E-06   3.143115E-03  -6.928899E-07   1.067252E-03  -3.188931E-07   6.964863E-03   2.129519E-06   1.091472E-02
  a    17   1.487076E-06  -6.558267E-07  -2.140241E-03  -3.196060E-07  -8.988619E-04  -1.609938E-06   3.338010E-08  -4.790292E-07
  a    18   4.884523E-07   7.548863E-04   6.507708E-08  -6.031493E-04   8.154731E-08  -1.135498E-03  -7.778097E-08  -2.333149E-04
  a    19  -1.519949E-02  -3.750608E-06  -1.069102E-07  -3.666435E-06   5.051921E-08  -5.983189E-06   1.097906E-03  -5.699504E-06
  a    20  -1.401073E-05   3.578878E-03  -2.213807E-07   3.735262E-03  -1.250016E-07   6.219659E-03   2.702728E-06   5.565221E-03
  a    21   2.022611E-06   3.064002E-07   6.896532E-03   1.971068E-07  -3.523175E-03  -4.407009E-08  -1.071987E-07  -3.240605E-07
  a    22   1.646455E-07   5.862028E-03   9.432890E-09  -4.264884E-03   2.469477E-07  -3.004907E-03  -1.861325E-07  -7.486950E-04
  a    23   3.834261E-08  -8.085415E-09   4.222923E-03   1.767090E-07   7.235286E-03   2.017306E-07   1.322840E-07  -6.827566E-08
  a    24  -6.728128E-07   2.212298E-03  -1.564047E-07   4.231929E-03   2.052406E-07  -6.647784E-03  -9.841179E-07  -3.568792E-03
  a    25   5.177664E-07   3.486557E-03  -2.756648E-07   2.083881E-03   1.374298E-08   2.361198E-03   5.093642E-07   1.622064E-03
  a    26  -1.047494E-06  -1.830288E-07  -1.095359E-03  -1.084497E-07   1.096179E-04  -6.738048E-08   6.885455E-08  -3.144761E-07
  a    27   8.794627E-04  -7.026205E-08  -1.095127E-07  -7.217691E-08  -3.062607E-08  -3.295707E-08  -1.752927E-04   4.301543E-08
  a    28   1.298688E-08   9.157656E-09   4.792422E-07   8.277717E-08   1.328182E-07  -7.362724E-09  -3.723963E-09   7.655659E-10
  a    29  -3.701608E-08  -1.902981E-08   1.932510E-07  -1.202730E-07   9.497651E-08  -6.774405E-08  -2.130520E-09  -8.007858E-09
  a    30  -6.494383E-07   4.795349E-04  -3.000145E-08   7.048095E-04   4.926783E-10   5.965725E-04   4.982580E-08   1.996471E-04
  a    31   9.969449E-07   1.723350E-08   6.967029E-04   1.311305E-08   8.824674E-05  -2.441236E-08   2.296371E-08   5.020079E-07
  a    32  -6.287428E-08  -6.250760E-04   8.042500E-09  -1.249875E-04   1.154231E-09  -3.077671E-04  -8.131018E-07  -2.525554E-03
  a    33  -6.212229E-07   5.808190E-04  -1.368236E-08  -2.444797E-04  -2.691803E-08   1.398466E-04   2.260775E-07   6.745926E-04
  a    34   7.605299E-05  -1.415132E-07   6.833702E-09   4.650788E-08  -1.318718E-08   3.978478E-08   7.813892E-05  -3.892851E-08
  a    35   4.360384E-07  -2.302414E-04  -1.397214E-08  -2.602250E-05   7.491729E-08   7.255403E-04   1.974637E-08   1.835920E-04
  a    36   9.642075E-08   2.336092E-08  -1.861470E-04  -1.989228E-08   6.385424E-04  -9.825121E-08  -4.522488E-09   3.781420E-08

               a    17        a    18        a    19        a    20        a    21        a    22        a    23        a    24
  a     3   3.930772E-07  -4.751477E-03  -1.376785E-07  -2.414058E-04   5.240491E-08   1.383223E-03   2.860067E-08  -4.639720E-05
  a     4  -2.529592E-06  -2.061632E-03  -1.121599E-05   1.115197E-02  -1.608559E-07  -4.546034E-03   2.023005E-07   3.053211E-03
  a     5   4.956794E-03  -3.575282E-08   5.208950E-07  -1.770150E-07  -1.403161E-02  -2.077605E-08  -1.007396E-02  -2.042666E-07
  a     6  -9.784082E-07   9.448834E-03   4.642353E-06  -4.987385E-03   3.600186E-07   8.909777E-03  -1.956371E-07   2.136889E-03
  a     7  -1.388070E-02   3.736964E-07  -1.714163E-07   2.846261E-07   1.751895E-02  -4.653866E-07  -1.245258E-02  -2.550976E-07
  a     8  -1.728482E-06   7.363552E-07  -1.906014E-02  -1.692000E-05   1.082708E-06   3.983092E-07  -2.835219E-07   4.705791E-07
  a     9   1.487076E-06   4.884523E-07  -1.519949E-02  -1.401073E-05   2.022611E-06   1.646455E-07   3.834261E-08  -6.728128E-07
  a    10  -6.558267E-07   7.548863E-04  -3.750608E-06   3.578878E-03   3.064002E-07   5.862028E-03  -8.085415E-09   2.212298E-03
  a    11  -2.140241E-03   6.507708E-08  -1.069102E-07  -2.213807E-07   6.896532E-03   9.432890E-09   4.222923E-03  -1.564047E-07
  a    12  -3.196060E-07  -6.031493E-04  -3.666435E-06   3.735262E-03   1.971068E-07  -4.264884E-03   1.767090E-07   4.231929E-03
  a    13  -8.988619E-04   8.154731E-08   5.051921E-08  -1.250016E-07  -3.523175E-03   2.469477E-07   7.235286E-03   2.052406E-07
  a    14  -1.609938E-06  -1.135498E-03  -5.983189E-06   6.219659E-03  -4.407009E-08  -3.004907E-03   2.017306E-07  -6.647784E-03
  a    15   3.338010E-08  -7.778097E-08   1.097906E-03   2.702728E-06  -1.071987E-07  -1.861325E-07   1.322840E-07  -9.841179E-07
  a    16  -4.790292E-07  -2.333149E-04  -5.699504E-06   5.565221E-03  -3.240605E-07  -7.486950E-04  -6.827566E-08  -3.568792E-03
  a    17   8.520336E-03  -4.529389E-08   2.675518E-07  -1.260339E-06  -1.344537E-03   2.255407E-07  -7.433572E-04   7.976686E-07
  a    18  -4.529389E-08   4.305237E-03   2.409349E-07  -4.140519E-04   6.009444E-08   1.678267E-03   5.176588E-09   4.462378E-04
  a    19   2.675518E-07   2.409349E-07   9.649637E-03   1.185072E-07  -5.796120E-08   9.417195E-07   1.084670E-07   1.288869E-06
  a    20  -1.260339E-06  -4.140519E-04   1.185072E-07   9.315017E-03   7.141644E-09  -1.005966E-03   7.259184E-08  -1.342302E-03
  a    21  -1.344537E-03   6.009444E-08  -5.796120E-08   7.141644E-09   6.667608E-03   1.501744E-08  -3.582530E-04   1.256815E-08
  a    22   2.255407E-07   1.678267E-03   9.417195E-07  -1.005966E-03   1.501744E-08   6.755039E-03   4.718814E-08   1.383922E-03
  a    23  -7.433572E-04   5.176588E-09   1.084670E-07   7.259184E-08  -3.582530E-04   4.718814E-08   9.319390E-03  -5.120966E-08
  a    24   7.976686E-07   4.462378E-04   1.288869E-06  -1.342302E-03   1.256815E-08   1.383922E-03  -5.120966E-08   8.199085E-03
  a    25  -4.088638E-07  -2.173984E-04  -1.694174E-06   1.799065E-03  -6.489799E-08   4.718061E-04   2.531656E-08   6.172985E-05
  a    26  -1.347662E-03   4.657717E-08   1.246098E-08  -1.933628E-08  -1.030473E-03  -1.248285E-08   3.767177E-04  -2.993447E-08
  a    27  -2.963325E-07   6.493427E-08  -4.432370E-04  -5.074559E-07  -8.680932E-08  -3.340242E-08  -1.719010E-08  -5.810109E-08
  a    28   7.196040E-07  -4.373113E-08  -2.399485E-09   3.931276E-08   4.288603E-07  -2.075581E-08  -7.866546E-08   4.341909E-08
  a    29   5.441066E-07  -3.884221E-08   7.211674E-10  -1.546871E-07   1.110940E-07   1.127572E-07   1.485230E-07  -2.045398E-08
  a    30   1.610277E-08   8.952438E-04  -1.286934E-06   1.459340E-03  -1.492746E-08  -4.635497E-04   3.173226E-08   2.300436E-05
  a    31   2.610188E-03  -6.816555E-08   3.328623E-08  -1.226179E-08   2.942684E-04   2.891483E-08   7.431450E-04   1.458490E-08
  a    32   4.595772E-07  -7.517924E-05   5.619715E-07  -5.349328E-04  -1.060651E-08   1.020096E-04  -5.083533E-08   3.495518E-04
  a    33  -1.681331E-07  -2.168379E-04  -1.209207E-06   1.540246E-03   1.808639E-09   3.971989E-04  -1.672584E-08  -1.901946E-04
  a    34   1.660919E-08  -7.371667E-09   8.278218E-04   4.598946E-07  -6.861568E-09  -1.289664E-07  -1.037847E-08  -3.228604E-08
  a    35  -2.603216E-08  -4.041121E-04  -9.710260E-07   9.518501E-04  -6.199805E-08  -4.997584E-04   5.366710E-08  -1.138418E-03
  a    36   2.340865E-04   4.873257E-08  -1.265131E-08  -1.131394E-07  -4.193349E-04   7.909460E-08   1.152495E-04   1.522076E-07

               a    25        a    26        a    27        a    28        a    29        a    30        a    31        a    32
  a     3   3.528719E-03  -1.843263E-07  -3.615193E-07  -1.080503E-07  -1.226259E-07   1.243839E-03  -1.018073E-07  -1.974032E-03
  a     4   3.666105E-03  -1.603277E-08  -7.761118E-08  -1.468037E-08   1.121910E-07  -8.314305E-04   2.215458E-07   2.751802E-03
  a     5   1.408136E-07   2.811997E-03   3.381289E-07  -5.700422E-07   7.488657E-07   1.505350E-08   4.435024E-03  -2.512876E-07
  a     6  -2.769624E-03   1.221747E-07   8.593768E-07   5.733249E-08  -1.481404E-07   1.979771E-03   9.898722E-10   7.102997E-04
  a     7  -2.157042E-07  -3.079456E-03   1.324532E-06   6.192186E-07  -1.041773E-07  -4.781091E-08  -2.894703E-03   1.318997E-07
  a     8  -4.788251E-07   7.346391E-07  -4.350580E-03   5.446415E-08   1.681084E-09   2.076626E-07   1.129713E-08  -4.386571E-08
  a     9   5.177664E-07  -1.047494E-06   8.794627E-04   1.298688E-08  -3.701608E-08  -6.494383E-07   9.969449E-07  -6.287428E-08
  a    10   3.486557E-03  -1.830288E-07  -7.026205E-08   9.157656E-09  -1.902981E-08   4.795349E-04   1.723350E-08  -6.250760E-04
  a    11  -2.756648E-07  -1.095359E-03  -1.095127E-07   4.792422E-07   1.932510E-07  -3.000145E-08   6.967029E-04   8.042500E-09
  a    12   2.083881E-03  -1.084497E-07  -7.217691E-08   8.277717E-08  -1.202730E-07   7.048095E-04   1.311305E-08  -1.249875E-04
  a    13   1.374298E-08   1.096179E-04  -3.062607E-08   1.328182E-07   9.497651E-08   4.926783E-10   8.824674E-05   1.154231E-09
  a    14   2.361198E-03  -6.738048E-08  -3.295707E-08  -7.362724E-09  -6.774405E-08   5.965725E-04  -2.441236E-08  -3.077671E-04
  a    15   5.093642E-07   6.885455E-08  -1.752927E-04  -3.723963E-09  -2.130520E-09   4.982580E-08   2.296371E-08  -8.131018E-07
  a    16   1.622064E-03  -3.144761E-07   4.301543E-08   7.655659E-10  -8.007858E-09   1.996471E-04   5.020079E-07  -2.525554E-03
  a    17  -4.088638E-07  -1.347662E-03  -2.963325E-07   7.196040E-07   5.441066E-07   1.610277E-08   2.610188E-03   4.595772E-07
  a    18  -2.173984E-04   4.657717E-08   6.493427E-08  -4.373113E-08  -3.884221E-08   8.952438E-04  -6.816555E-08  -7.517924E-05
  a    19  -1.694174E-06   1.246098E-08  -4.432370E-04  -2.399485E-09   7.211674E-10  -1.286934E-06   3.328623E-08   5.619715E-07
  a    20   1.799065E-03  -1.933628E-08  -5.074559E-07   3.931276E-08  -1.546871E-07   1.459340E-03  -1.226179E-08  -5.349328E-04
  a    21  -6.489799E-08  -1.030473E-03  -8.680932E-08   4.288603E-07   1.110940E-07  -1.492746E-08   2.942684E-04  -1.060651E-08
  a    22   4.718061E-04  -1.248285E-08  -3.340242E-08  -2.075581E-08   1.127572E-07  -4.635497E-04   2.891483E-08   1.020096E-04
  a    23   2.531656E-08   3.767177E-04  -1.719010E-08  -7.866546E-08   1.485230E-07   3.173226E-08   7.431450E-04  -5.083533E-08
  a    24   6.172985E-05  -2.993447E-08  -5.810109E-08   4.341909E-08  -2.045398E-08   2.300436E-05   1.458490E-08   3.495518E-04
  a    25   3.753055E-03  -6.782152E-08   7.123361E-08  -2.184769E-08   1.090219E-08  -9.368485E-05  -6.189548E-08  -3.967748E-04
  a    26  -6.782152E-08   1.971171E-03  -4.910880E-07   8.836066E-07   2.362821E-07   9.918127E-09  -9.051683E-04   3.379822E-08
  a    27   7.123361E-08  -4.910880E-07   5.261209E-03  -1.451715E-08   5.719525E-09   3.179523E-08  -1.498655E-07   1.773152E-08
  a    28  -2.184769E-08   8.836066E-07  -1.451715E-08   3.765201E-03   6.721664E-04   6.084752E-08   2.929927E-07   6.953313E-09
  a    29   1.090219E-08   2.362821E-07   5.719525E-09   6.721664E-04   3.236389E-03   1.219818E-07   3.531533E-08  -2.643865E-08
  a    30  -9.368485E-05   9.918127E-09   3.179523E-08   6.084752E-08   1.219818E-07   1.920253E-03   2.348995E-08   7.729295E-05
  a    31  -6.189548E-08  -9.051683E-04  -1.498655E-07   2.929927E-07   3.531533E-08   2.348995E-08   2.475218E-03  -1.234415E-08
  a    32  -3.967748E-04   3.379822E-08   1.773152E-08   6.953313E-09  -2.643865E-08   7.729295E-05  -1.234415E-08   2.170822E-03
  a    33  -5.445742E-04   6.307547E-08   1.778486E-07  -1.631621E-08  -5.824444E-08   5.004242E-04  -9.069666E-09  -2.216628E-04
  a    34   1.606614E-07  -1.077544E-07   7.185735E-04  -9.442756E-09   3.928424E-09  -1.096746E-07  -5.799971E-09   2.142523E-08
  a    35   7.524728E-04  -1.142772E-07  -2.497032E-08  -1.123265E-09   7.023934E-09  -1.057780E-04   2.685212E-08   4.284382E-05
  a    36  -1.422205E-07  -6.596222E-04  -8.121189E-08   2.677502E-07   4.788876E-08   1.893207E-08   1.694814E-04  -1.553109E-08

               a    33        a    34        a    35        a    36
  a     3   1.493906E-03  -4.127698E-07  -6.985103E-04   9.108913E-08
  a     4   4.469467E-03  -1.120108E-06  -1.141256E-03   2.456921E-07
  a     5  -1.492860E-07   4.567027E-07   4.670115E-07   3.579221E-03
  a     6  -5.835565E-04  -2.096555E-07   1.301790E-03  -1.939498E-07
  a     7  -4.524771E-08   5.150717E-07   7.089996E-07   5.216734E-03
  a     8  -1.089612E-06  -6.404741E-03   5.780877E-08   1.831761E-07
  a     9  -6.212229E-07   7.605299E-05   4.360384E-07   9.642075E-08
  a    10   5.808190E-04  -1.415132E-07  -2.302414E-04   2.336092E-08
  a    11  -1.368236E-08   6.833702E-09  -1.397214E-08  -1.861470E-04
  a    12  -2.444797E-04   4.650788E-08  -2.602250E-05  -1.989228E-08
  a    13  -2.691803E-08  -1.318718E-08   7.491729E-08   6.385424E-04
  a    14   1.398466E-04   3.978478E-08   7.255403E-04  -9.825121E-08
  a    15   2.260775E-07   7.813892E-05   1.974637E-08  -4.522488E-09
  a    16   6.745926E-04  -3.892851E-08   1.835920E-04   3.781420E-08
  a    17  -1.681331E-07   1.660919E-08  -2.603216E-08   2.340865E-04
  a    18  -2.168379E-04  -7.371667E-09  -4.041121E-04   4.873257E-08
  a    19  -1.209207E-06   8.278218E-04  -9.710260E-07  -1.265131E-08
  a    20   1.540246E-03   4.598946E-07   9.518501E-04  -1.131394E-07
  a    21   1.808639E-09  -6.861568E-09  -6.199805E-08  -4.193349E-04
  a    22   3.971989E-04  -1.289664E-07  -4.997584E-04   7.909460E-08
  a    23  -1.672584E-08  -1.037847E-08   5.366710E-08   1.152495E-04
  a    24  -1.901946E-04  -3.228604E-08  -1.138418E-03   1.522076E-07
  a    25  -5.445742E-04   1.606614E-07   7.524728E-04  -1.422205E-07
  a    26   6.307547E-08  -1.077544E-07  -1.142772E-07  -6.596222E-04
  a    27   1.778486E-07   7.185735E-04  -2.497032E-08  -8.121189E-08
  a    28  -1.631621E-08  -9.442756E-09  -1.123265E-09   2.677502E-07
  a    29  -5.824444E-08   3.928424E-09   7.023934E-09   4.788876E-08
  a    30   5.004242E-04  -1.096746E-07  -1.057780E-04   1.893207E-08
  a    31  -9.069666E-09  -5.799971E-09   2.685212E-08   1.694814E-04
  a    32  -2.216628E-04   2.142523E-08   4.284382E-05  -1.553109E-08
  a    33   1.640552E-03  -2.634381E-08   1.773528E-04  -2.069896E-08
  a    34  -2.634381E-08   1.469707E-03  -3.807296E-08  -1.023511E-09
  a    35   1.773528E-04  -3.807296E-08   1.393845E-03  -2.370450E-08
  a    36  -2.069896E-08  -1.023511E-09  -2.370450E-08   1.182206E-03

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98158353     1.97374362     1.97300751     1.96898227     1.95950785     1.48357336
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.50277583     0.02782597     0.01824793     0.01804083     0.01592310     0.01532576     0.00940582     0.00878891
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00553880     0.00530571     0.00492462     0.00422310     0.00410074     0.00351467     0.00277849     0.00232101
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00216266     0.00157147     0.00124342     0.00123155     0.00096427     0.00082115     0.00081168     0.00047573
              MO    33       MO    34       MO    35       MO    36
  occ(*)=     0.00041042     0.00039524     0.00028015     0.00019282


 total number of electrons =   16.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a   partial gross atomic populations
   ao class       1a         2a         3a         4a         5a         6a  
    N1_ s       1.998472  -0.000030   1.426462   0.124309   0.000000   0.000000
    N1_ p       0.000000  -0.000292   0.002615   0.175779   0.497542   0.885317
    N1_ d       0.000000  -0.000105   0.004284   0.003111   0.000691   0.015304
    C1_ s       0.000539   1.999452   0.256486   1.019043   0.000000   0.000000
    C1_ p       0.000476   0.000003  -0.006104   0.123387   0.921758   0.317141
    C1_ d      -0.000081   0.000001   0.005607   0.005721   0.000258   0.017276
    H1_ s       0.000301   0.000003   0.118959   0.042038   0.058083   0.229328
    H2_ s       0.000301   0.000003   0.118959   0.042031   0.058080   0.229339
    H3_ s      -0.000004   0.000483   0.027157   0.219172   0.218295   0.137632
    H4_ s      -0.000004   0.000483   0.027159   0.219152   0.218300   0.137647
 
   ao class       7a         8a         9a        10a        11a        12a  
    N1_ s       0.036080   0.000000  -0.000000   0.002860   0.000000   0.003892
    N1_ p       0.962866   1.417799   0.007015   0.007058   0.010109   0.003272
    N1_ d       0.006360  -0.000741   0.005108   0.001516   0.000434   0.000614
    C1_ s       0.015838   0.000000   0.000000   0.004059   0.000000   0.001092
    C1_ p       0.686786   0.050898   0.490611   0.011765   0.000305   0.000501
    C1_ d       0.015773   0.015617   0.000042   0.000470   0.000165   0.000056
    H1_ s       0.058276   0.000000   0.000000  -0.000010   0.003600   0.003412
    H2_ s       0.058276   0.000000   0.000000  -0.000010   0.003602   0.003410
    H3_ s       0.059626   0.000000   0.000000   0.000059   0.000016   0.000896
    H4_ s       0.059627   0.000000   0.000000   0.000059   0.000016   0.000896
 
   ao class      13a        14a        15a        16a        17a        18a  
    N1_ s       0.000000   0.001025   0.000000   0.000000   0.000328   0.000000
    N1_ p       0.000070   0.000311   0.002025   0.009132   0.001351   0.000031
    N1_ d       0.001073   0.000196   0.002016   0.000004   0.002246   0.003638
    C1_ s       0.000000   0.003919   0.000000   0.000000  -0.000124   0.000000
    C1_ p       0.006869   0.001480   0.003552  -0.000482   0.000806   0.000164
    C1_ d       0.000616   0.000134   0.001718   0.000136   0.000887   0.001473
    H1_ s       0.000244   0.000782   0.000056   0.000000   0.000010   0.000000
    H2_ s       0.000244   0.000783   0.000056   0.000000   0.000010   0.000000
    H3_ s       0.003404   0.003348  -0.000008   0.000000   0.000013  -0.000000
    H4_ s       0.003404   0.003348  -0.000008  -0.000000   0.000013   0.000000
 
   ao class      19a        20a        21a        22a        23a        24a  
    N1_ s       0.001018   0.000000   0.000053   0.000000   0.000000   0.000000
    N1_ p       0.000187   0.000000   0.000123  -0.000129   0.000000   0.000071
    N1_ d       0.002920   0.004072   0.000313   0.000047   0.000099   0.000863
    C1_ s       0.000236   0.000000   0.001286  -0.000000   0.000000   0.000000
    C1_ p      -0.000016   0.000000   0.000339   0.003554   0.000000   0.000344
    C1_ d       0.000188   0.000151   0.001331   0.000041   0.002679   0.000343
    H1_ s       0.000111  -0.000000  -0.000007  -0.000000   0.000000   0.000090
    H2_ s       0.000111  -0.000000  -0.000007  -0.000000   0.000000   0.000090
    H3_ s       0.000085   0.000000   0.000336   0.000000   0.000000   0.000261
    H4_ s       0.000085   0.000000   0.000336   0.000000   0.000000   0.000260
 
   ao class      25a        26a        27a        28a        29a        30a  
    N1_ s       0.000069   0.000199   0.000000   0.000000   0.000050   0.000118
    N1_ p       0.000305   0.000020   0.000142   0.000027   0.000241   0.000084
    N1_ d       0.000269   0.000288   0.000110   0.000358   0.000230   0.000134
    C1_ s       0.000215  -0.000126   0.000000   0.000000   0.000025  -0.000008
    C1_ p       0.000367   0.000041   0.000145   0.000003   0.000166  -0.000018
    C1_ d       0.000883   0.000308   0.000491   0.000843   0.000103   0.000174
    H1_ s      -0.000014   0.000284   0.000133   0.000000   0.000038   0.000057
    H2_ s      -0.000014   0.000284   0.000133   0.000000   0.000038   0.000057
    H3_ s       0.000042   0.000137   0.000045   0.000000   0.000037   0.000111
    H4_ s       0.000042   0.000137   0.000045   0.000000   0.000037   0.000111
 
   ao class      31a        32a        33a        34a        35a        36a  
    N1_ s       0.000000   0.000076   0.000000   0.000004   0.000000   0.000002
    N1_ p      -0.000012   0.000001   0.000155   0.000046   0.000005   0.000005
    N1_ d       0.000159   0.000010   0.000011   0.000018   0.000000   0.000008
    C1_ s      -0.000000   0.000076   0.000000  -0.000001   0.000000   0.000068
    C1_ p      -0.000024   0.000107   0.000061   0.000029   0.000046   0.000033
    C1_ d       0.000160   0.000086   0.000000   0.000031   0.000006   0.000048
    H1_ s       0.000115   0.000048   0.000105   0.000096   0.000025  -0.000000
    H2_ s       0.000115   0.000048   0.000106   0.000096   0.000025  -0.000000
    H3_ s       0.000149   0.000011  -0.000014   0.000038   0.000087   0.000014
    H4_ s       0.000148   0.000011  -0.000014   0.000038   0.000087   0.000014


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.594987   3.302073   0.516163   0.516165   0.671425   0.671427
      p         3.983269   2.615093   0.000000   0.000000   0.000000   0.000000
      d         0.055661   0.073735   0.000000   0.000000   0.000000   0.000000
    total       7.633918   5.990902   0.516163   0.516165   0.671425   0.671427
 

 Total number of electrons:   16.00000000

 accstate=                     4
 accpdens=                     3
logrecs(*)=   1   2logvrecs(*)=   2   3
 item #                     3 suffix=:.trd2to3:
 computing final density
 computing MRCISD density
 densi: densityinfo%a4den=   1.00000000000000     
 =========== Executing IN-CORE method ==========
1e-transition density (   2,   3)
--------------------------------------------------------------------------------
  1e-transition density (root #    2 -> root #   3)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
  1e-transition density (root #    2 -> root #   3)
--------------------------------------------------------------------------------
================================================================================
   DYZ=     215  DYX=     335  DYW=     401
   D0Z=     147  D0Y=     776  D0X=     489  D0W=     695
  DDZI=     189 DDYI=     614 DDXI=     452 DDWI=     599
  DDZE=       0 DDYE=      92 DDXE=      71 DDWE=      98
================================================================================
2e-transition density (   2,   3)
--------------------------------------------------------------------------------
  2e-transition density (root #    2 -> root #   3)
--------------------------------------------------------------------------------
 call to prpd23 iwx=1
 call to prpd23 iwx=2
 starting prpd24 .... 
 starting prpd24 .... 
================================================================================
   DYZ=     215  DYX=    4713  DYW=    5711
   D0Z=     147  D0Y=    7060  D0X=    3958  D0W=    5934
  DDZI=     735 DDYI=    2366 DDXI=    1673 DDWI=    2150
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
  entering dump_fourext .... 
 maximum diagonal element=  8.246287230909838E-005
