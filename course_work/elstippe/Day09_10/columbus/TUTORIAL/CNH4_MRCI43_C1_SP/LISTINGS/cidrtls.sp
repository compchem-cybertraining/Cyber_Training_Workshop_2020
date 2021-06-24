 
 program cidrt 7.0  

 distinct row table construction, reference csf selection, and internal
 walk selection for multireference single- and double-excitation
configuration interaction.

 references:  r. shepard, i. shavitt, r. m. pitzer, d. c. comeau, m. pepper
                  h. lischka, p. g. szalay, r. ahlrichs, f. b. brown, and
                  j.-g. zhao, int. j. quantum chem. symp. 22, 149 (1988).
              h. lischka, r. shepard, f. b. brown, and i. shavitt,
                  int. j. quantum chem. symp. 15, 91 (1981).

 based on the initial version by  Ron Shepard

 extended for spin-orbit CI calculations ( Russ Pitzer, OSU)

 and large active spaces (Thomas Müller, FZ(21 Juelich)

 This Version of Program CIDRT is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de

*********************** File revision status: ***********************
* cidrt1.F9 Revision: 1.1.6.2           Date: 2013/04/11 14:37:29   * 
* cidrt2.F9 Revision: 1.1.6.6           Date: 2015/02/26 17:04:32   * 
* cidrt3.F9 Revision: 1.1.6.2           Date: 2013/04/11 14:37:29   * 
* cidrt4.F9 Revision: 1.1.6.2           Date: 2013/04/11 14:37:29   * 
********************************************************************

 workspace allocation parameters: lencor= 222822400 mem1=         0 ifirst=         1
 expanded "keystrokes" are being written to file:
 /panfs/panfs.cbls.ccr.buffalo.edu/scratch/grp-cyberwksp21/elstippe/columbus/TUTO
 Spin-Orbit CI Calculation?(y,[n]) Spin-Free Calculation
 
 input the spin multiplicity [  0]: spin multiplicity, smult            :   1    singlet 
 input the total number of electrons [  0]: total number of electrons, nelt     :    16
 input the number of irreps (1:8) [  0]: point group dimension, nsym         :     1
 enter symmetry labels:(y,[n]) enter 1 labels (a4):
 enter symmetry label, default=   1
 symmetry labels: (symmetry, slabel)
 ( 1,  a  ) 
 input nmpsy(*):
 nmpsy(*)=        36
 
   symmetry block summary
 block(*)=         1
 slabel(*)=      a  
 nmpsy(*)=        36
 
 total molecular orbitals            :    36
 input the molecular spatial symmetry (irrep 1:nsym) [  0]: state spatial symmetry label        :  a  
 
 input the frozen core orbitals (sym(i),rmo(i),i=1,nfct):
 total frozen core orbitals, nfct    :     2
 
 fcorb(*)=         1   2
 slabel(*)=      a   a  
 
 number of frozen core orbitals      :     2
 number of frozen core electrons     :     4
 number of internal electrons        :    12
 
 input the frozen virtual orbitals (sym(i),rmo(i),i=1,nfvt):
 total frozen virtual orbitals, nfvt :     0

 no frozen virtual orbitals entered
 
 input the internal orbitals (sym(i),rmo(i),i=1,niot):
 niot                                :     7
 
 modrt(*)=         3   4   5   6   7   8   9
 slabel(*)=      a   a   a   a   a   a   a  
 
 total number of orbitals            :    36
 number of frozen core orbitals      :     2
 number of frozen virtual orbitals   :     0
 number of internal orbitals         :     7
 number of external orbitals         :    27
 
 orbital-to-level mapping vector
 map(*)=          -1  -1  28  29  30  31  32  33  34   1   2   3   4   5   6
                   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
                  22  23  24  25  26  27
 
 input the number of ref-csf doubly-occupied orbitals [  0]: (ref) doubly-occupied orbitals      :     4
 
 no. of internal orbitals            :     7
 no. of doubly-occ. (ref) orbitals   :     4
 no. active (ref) orbitals           :     3
 no. of active electrons             :     4
 
 input the active-orbital, active-electron occmnr(*):
   7  8  9
 input the active-orbital, active-electron occmxr(*):
   7  8  9
 
 actmo(*) =        7   8   9
 occmnr(*)=        0   0   4
 occmxr(*)=        4   4   4
 reference csf cumulative electron occupations:
 modrt(*)=         3   4   5   6   7   8   9
 occmnr(*)=        2   4   6   8   8   8  12
 occmxr(*)=        2   4   6   8  12  12  12
 
 input the active-orbital bminr(*):
   7  8  9
 input the active-orbital bmaxr(*):
   7  8  9
 reference csf b-value constraints:
 modrt(*)=         3   4   5   6   7   8   9
 bminr(*)=         0   0   0   0   0   0   0
 bmaxr(*)=         0   0   0   0   4   4   4
 input the active orbital smaskr(*):
   7  8  9
 modrt:smaskr=
   3:1000   4:1000   5:1000   6:1000   7:1111   8:1111   9:1111
 
 input the maximum excitation level from the reference csfs [  2]: maximum excitation from ref. csfs:  :     2
 number of internal electrons:       :    12
 
 input the internal-orbital mrsdci occmin(*):
   3  4  5  6  7  8  9
 input the internal-orbital mrsdci occmax(*):
   3  4  5  6  7  8  9
 mrsdci csf cumulative electron occupations:
 modrt(*)=         3   4   5   6   7   8   9
 occmin(*)=        0   0   0   0   0   0  10
 occmax(*)=       12  12  12  12  12  12  12
 
 input the internal-orbital mrsdci bmin(*):
   3  4  5  6  7  8  9
 input the internal-orbital mrsdci bmax(*):
   3  4  5  6  7  8  9
 mrsdci b-value constraints:
 modrt(*)=         3   4   5   6   7   8   9
 bmin(*)=          0   0   0   0   0   0   0
 bmax(*)=         12  12  12  12  12  12  12
 
 input the internal-orbital smask(*):
   3  4  5  6  7  8  9
 modrt:smask=
   3:1111   4:1111   5:1111   6:1111   7:1111   8:1111   9:1111
 
 internal orbital summary:
 block(*)=         1   1   1   1   1   1   1
 slabel(*)=      a   a   a   a   a   a   a  
 rmo(*)=           3   4   5   6   7   8   9
 modrt(*)=         3   4   5   6   7   8   9
 
 reference csf info:
 occmnr(*)=        2   4   6   8   8   8  12
 occmxr(*)=        2   4   6   8  12  12  12
 
 bminr(*)=         0   0   0   0   0   0   0
 bmaxr(*)=         0   0   0   0   4   4   4
 
 
 mrsdci csf info:
 occmin(*)=        0   0   0   0   0   0  10
 occmax(*)=       12  12  12  12  12  12  12
 
 bmin(*)=          0   0   0   0   0   0   0
 bmax(*)=         12  12  12  12  12  12  12
 

 a priori removal of distinct rows:

 input the level, a, and b values for the vertices 
 to be removed (-1/ to end).

 input level, a, and b (-1/ to end):
 no vertices marked for removal
 
 impose generalized interacting space restrictions?(y,[n]) generalized interacting space restrictions will be imposed.
 multp(*)=
  hmult                     0
 lxyzir   0   0   0
 symmetry of spin functions (spnir)
       --------------------------Ms ----------------------------
   S     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
   1     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   2     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   3     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   4     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   5     1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   6     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   7     0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   8     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   9     1  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0
  10     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  11     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
  12     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  13     1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  0  0  0
  14     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  15     0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0
  16     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  17     1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0  0
  18     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  19     0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0

 number of rows in the drt :  42
     4 arcs removed due to generalized interacting space restrictions.

 manual arc removal step:


 input the level, a, b, and step values 
 for the arcs to be removed (-1/ to end).

 input the level, a, b, and step (-1/ to end):
 remarc:   0 arcs removed out of   0 specified.

 xbarz=          28
 xbary=         112
 xbarx=         141
 xbarw=         196
        --------
 nwalk=         477
 input the range of drt levels to print (l1,l2):
 levprt(*)        -1   0

 reference-csf selection step 1:
 total number of z-walks in the drt, nzwalk=      28

 input the list of allowed reference symmetries:
 allowed reference symmetries:             1
 allowed reference symmetry labels:      a  
 keep all of the z-walks as references?(y,[n]) all z-walks are initially deleted.
 
 generate walks while applying reference drt restrictions?([y],n) reference drt restrictions will be imposed on the z-walks.
 
 impose additional orbital-group occupation restrictions?(y,[n]) 
 apply primary reference occupation restrictions?(y,[n]) 
 manually select individual walks?(y,[n])
 step 1 reference csf selection complete.
        6 csfs initially selected from      28 total walks.

 beginning step-vector based selection.
 enter [internal_orbital_step_vector/disposition] pairs:

 enter internal orbital step vector, (-1/ to end):
   3  4  5  6  7  8  9

 step 2 reference csf selection complete.
        6 csfs currently selected from      28 total walks.

 beginning numerical walk based selection.
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end:

 input reference walk number (0 to end) [  0]:
 numerical walk-number based selection complete.
        6 reference csfs selected from      28 total z-walks.
 
 input the reference occupations, mu(*):
 reference occupations:
 mu(*)=            2   2   2   2   0   0   0
 
 interacting space determination:
 checking diagonal loops...
 checking 2-internal loops...
 checking 3-internal loops...
 checking 4-internal loops...
 limint: nvalw=                    28                    92
                    71                    98
 post-limint icd(*)=                     1                     1
                    31                   510                     0
                     0
 
 this is an obsolete prompt.(y,[n])
 final mrsdci walk selection step:

 nvalw(*)=      28      92      71      98 nvalwt=     289

 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input mrsdci walk number (0 to end) [  0]:
 end of manual mrsdci walk selection.
 number added=   0 number removed=   0

 nvalw(*)=      28      92      71      98 nvalwt=     289

 lprune input numv1,nwalk=                   289                   477
 lprune input xbar(1,1),nref=                    28                     6

 lprune: l(*,*,*) pruned with nwalk=     477 nvalwt=     289=  28  92  71  98
 lprune:  z-drt, nprune=    57
 lprune:  y-drt, nprune=    39
 lprune: wx-drt, nprune=    39

 xbarz=          28
 xbary=          92
 xbarx=          71
 xbarw=          98
        --------
 nwalk=         289
 levprt(*)        -1   0

 beginning the reference csf index recomputation...

     iref   iwalk  step-vector
   ------  ------  ------------
        1       1  3333330
        2       2  3333312
        3       3  3333303
        4       4  3333132
        5       5  3333123
        6       6  3333033
 indx01:     6 elements set in vec01(*)

 beginning the valid upper walk index recomputation...
 indx01:   289 elements set in vec01(*)

 beginning the final csym(*) computation...

  number of valid internal walks of each symmetry:

       a  
      ----
 z              28
 y              92
 x              71
 w              98

 csfs grouped by internal walk symmetry:

       a  
      ----
 z              28
 y            2484
 x           24921
 w           37044

 total csf counts:
 z-vertex:              28
 y-vertex:            2484
 x-vertex:           24921
 w-vertex:           37044
           --------
 total:           64477
 
 input a title card, default=cidrt_title
 title card:
  cidrt_title                                                                   
  
 
 input a drt file name, default=cidrtfl
 drt and indexing arrays will be written to file:
 /panfs/panfs.cbls.ccr.buffalo.edu/scratch/grp-cyberwksp21/elstippe/columbus/TUTO
 
 write the drt file?([y],n) drt file is being written...
 wrtstr:  a  
nwalk=     289 cpos=      36 maxval=    9 cmprfactor=   87.54 %.
nwalk=     289 cpos=       3 maxval=   99 cmprfactor=   97.92 %.
nwalk=     289 cpos=       1 maxval=  999 cmprfactor=   98.96 %.
nwalk=     289 cpos=       1 maxval= 9999 cmprfactor=   98.62 %.
 compressed with: nwalk=     289 cpos=       1 maxval=  999 cmprfactor=   98.96 %.
initial index vector length:       289
compressed index vector length:         1reduction:  99.65%
nwalk=      28 cpos=       4 maxval=    9 cmprfactor=   85.71 %.
nwalk=      28 cpos=       2 maxval=   99 cmprfactor=   85.71 %.
 compressed with: nwalk=      28 cpos=       4 maxval=    9 cmprfactor=   85.71 %.
initial ref vector length:        28
compressed ref vector length:         4reduction:  85.71%
