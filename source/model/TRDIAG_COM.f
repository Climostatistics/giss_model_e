#include "rundeck_opts.h"

      MODULE TRACER_DIAG_COM
!@sum Tracer diagnostic arrays
!@+    Mostly tracer independent, but this may deped on applications
!@auth Jean Lerner
!ver   1.0
      USE MODEL_COM, only: im,jm,lm
      USE DAGCOM, only: npts !npts are conservation quantities
      USE TRACER_COM, only: ntm
#ifdef TRACERS_ON
     *     , ntsurfsrcmax, nt3Dsrcmax
#endif
      IMPLICIT NONE
      SAVE

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>

!@dbparam to_volume_MixRat: For printout of tracer concentration
!@+   to_volume_MixRat=1: printout is in Volume Mixing Ratio
!@+   to_volume_MixRat=0: printout is in Mass Mixing Ratio
#ifdef TRACERS_SPECIAL_Shindell
      INTEGER, DIMENSION(NTM) :: to_volume_MixRat=1
#else
      INTEGER, DIMENSION(NTM) :: to_volume_MixRat=0
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@dbparam to_per_mil For printout of tracer concentration in permil
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
#endif
#ifdef TRACERS_ON
!@var MMR_to_VMR: converts tracer mass mixing ratio to volume mr
      REAL*8, DIMENSION(NTM) :: MMR_to_VMR

C**** TAIJLN
!@var TAIJLN 3D tracer diagnostics (all tracers)
      real*8, dimension(im,jm,lm,ntm) :: taijln
!@var SNAME_IJT, UNITS_IJT: Names and units of lat-sigma tracer IJ diags
      character(len=30), dimension(lm,ntm) :: sname_ijt,units_ijt
!@var LNAME_IJT: descriptions of tracer IJ diagnostics
      character(len=80), dimension(lm,ntm) :: lname_ijt = 'unused'
!@var SCALE_IJT: printout scaling factor for tracer IJ diagnostics
      REAL*8, dimension(lm,ntm) :: scale_ijt
!@var IR_IJT: range index of IJ diagnosts
      integer, dimension(ntm) :: ir_ijt
!@var IA_IJT: idacc-number for tracer IJ diags (same for all tracers)
      integer ia_ijt
!@var IJTC_POWER: power of 10 used for tracer IJ concentration diags
      integer, dimension(ntm) :: ijtc_power
!@var IJTM_POWER: power of 10 used for tracer IJ mass unit diags
      integer, dimension(ntm) :: ijtm_power
!@var IJT_XX Names for TAIJLN diagnostics (Not currently used)

C**** TAIJN
!@param KTAIJ number of 3D diagnostics for each tracer
#ifdef TRACERS_WATER
      integer, parameter :: ktaij=11
#else
      integer, parameter :: ktaij=3
#endif
!@var IJT_XX names for taijn diagnostics
      integer tij_conc,tij_surf,tij_mass
#ifdef TRACERS_WATER
!@var IJT_XX names for water-based taijn diagnostics
      integer tij_rvr,tij_seaice,tij_prec,tij_evap,tij_grnd,tij_lk1
     *     ,tij_lk2,tij_soil
#endif
!@var TAIJN lat/lon tracer diagnostics (all tracers)
      real*8, dimension(im,jm,ktaij,ntm) :: taijn
!@var SCALE_TIJ: printout scaling factor for tracer IJK diagnostics
      REAL*8, dimension(ktaij,ntm) :: scale_tij
!@var SNAME_TIJ,UNITS_TIJ: Names and units of lat-sigma tracer diags
      character(len=30), dimension(ktaij,ntm) :: sname_tij,units_tij
!@var LNAME_TIJ: descriptions of tracer IJK diags
      character(len=80), dimension(ktaij,ntm) :: lname_tij = 'unused'

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
!@parm KTAIJS number of special lat/lon tracer diagnostics
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ktaijs=26
#else
      integer, parameter :: ktaijs=25
#endif
!@var TAIJS  lat/lon special tracer diagnostics; sources, sinks, etc.
      REAL*8, DIMENSION(IM,JM,ktaijs) :: TAIJS
!@var ijts_source tracer independent array for TAIJS surface src. diags
      INTEGER ijts_source(ntsurfsrcmax,ntm)
!@var SNAME_IJTS, UNITS_IJTS: Names & units of lat-sigma tracer diags
      character(len=30), dimension(ktaijs) :: sname_ijts,units_ijts
!@var LNAME_IJTS: descriptions of tracer IJTS diags
      character(len=80), dimension(ktaijs) :: lname_ijts = 'unused'
!@var SCALE_IJTS: printout scaling factor for tracer IJTS diagnostics
      REAL*8, dimension(ktaijs) :: scale_ijts
!@var IA_IJTS: idacc-number for tracer source/sink IJ diags
      integer ia_ijts(ktaijs)
!@var ijts_power: power of 10 used for tracer IJ source/sink diags
      INTEGER, DIMENSION(ktaijs) :: ijts_power
!@var ijts_index: tracer index associated with a TAIJS diagnostic
      INTEGER, DIMENSION(ktaijs) :: ijts_index

C**** TAJLN
!@parm ktajl,ktajlx number of TAJL tracer diagnostics;
!@+          ktajlx includes composits
      INTEGER, PARAMETER :: ktajl=8, ktajlx=ktajl+2
!@var TAJLN  vertical tracer diagnostics (all tracers)
      REAL*8, DIMENSION(JM,LM,ktajlx,NTM) :: TAJLN
!@var jlnt_xx Names for TAJLN diagnostics
      INTEGER jlnt_conc,jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,
     &  jlnt_vt_mm,jlnt_mc,jlnt_turb,jlnt_lscond
!@var SNAME_JLN: Names of lat-sigma tracer JL diagnostics
      character(len=30), dimension(ktajlx,ntm) :: sname_jln
!@var LNAME_JLN,UNITS_JLN: descriptions/units of tracer JL diagnostics
      character(len=50), dimension(ktajlx,ntm) :: lname_jln = 'unused'
     *     ,units_jln
!@var SCALE_JLQ: printout scaling factor for tracer JL diagnostics
      REAL*8, dimension(ktajlx) :: scale_jlq
!@var IA_JLQ,JGRID_JLQ: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajlx) :: ia_jlq,jgrid_jlq
!@var JLQ_POWER: power of 10 used for tracer JL diagnostics
!@+        It is associated with a specific physical process
      integer jlq_power(ktajlx)
      integer ktajlp !@var helps keeps track of
!@var scale_jln: Scale for jl maps
      REAL*8, DIMENSION(ntm) :: scale_jln

C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>
!@parm ktajls number of source/sink TAJLS tracer diagnostics;
#ifdef TRACERS_SPECIAL_Shindell
      INTEGER, PARAMETER :: ktajls=58
#else
      INTEGER, PARAMETER :: ktajls=34
#endif
!@var TAJLS  JL special tracer diagnostics for sources, sinks, etc
      REAL*8, DIMENSION(JM,LM,ktajls) :: TAJLS
!@var jls_source tracer independent array for TAJLS surface src. diags
      INTEGER jls_source(ntsurfsrcmax,ntm)
!@var jls_3Dsource tracer independent array for TAJLS 3D source diags
      INTEGER jls_3Dsource(nt3Dsrcmax,ntm)
!@var jls_decay tracer independent array for radioactive sinks
      INTEGER, DIMENSION(NTM) :: jls_decay
!@var jls_grav tracer independent array for grav. settling sink
      INTEGER, DIMENSION(NTM) :: jls_grav
!@var jls_index: Number of source/sink tracer JL diagnostics/tracer
      integer, dimension(ktajls) :: jls_index
!@var SNAME_JLS: Names of lat-sigma tracer JL sources/sinks
      character(len=30), dimension(ktajls) :: sname_jls
!@var LNAME_JLS,UNITS_JLS: descriptions/units of tracer JLS diags
      character(len=50), dimension(ktajls) :: lname_jls = 'unused'
     *     ,units_jls
!@var SCALE_JLS: printout scaling factor for tracer JLS diagnostics
      REAL*8, dimension(ktajls) :: scale_jls
!@var IA_JLS,JGRID_JLS: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajls) :: ia_jls,jgrid_jls
!@var JLS_POWER: power of 10 used for tracer JLS diagnostics
      integer jls_power(ktajls)
!@var jls_ltop: Top layer for this diagnostic
      integer jls_ltop(ktajls)

C**** TCONSRV
!@parm NTCONS Maximum Number of special tracer conservation points
      INTEGER, PARAMETER :: ntcons=20
!@parm KTCON total number of conservation diagnostics for tracers
      INTEGER, PARAMETER :: KTCON=npts+ntcons+2
!@var TCONSRV conservation diagnostics for tracers
      REAL*8, DIMENSION(JM,ktcon,ntm) :: TCONSRV
!@var SCALE_TCON scales for tracer conservation diagnostics
      REAL*8, DIMENSION(ktcon,ntm) :: SCALE_TCON
!@var TITLE_TCON titles for tracer conservation diagnostics
      CHARACTER*38, DIMENSION(ktcon,ntm) :: TITLE_TCON
!@var IA_TCON IDACC numbers for tracer conservation diagnostics
      INTEGER, DIMENSION(ktcon,ntm) ::  IA_TCON
!@var NSUM_TCON indices for summation of conservation diagnostics
      INTEGER, DIMENSION(ktcon,ntm) :: NSUM_TCON
!@var NOFMT indices for TCONSRV array
      INTEGER, DIMENSION(ktcon,ntm) :: NOFMT
!@var CONPTS names of special processes for tracer conservation diags
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS
!@var kt_power_inst,kt_power_change: Exponents for tracer conservation
      INTEGER, DIMENSION(ntm) :: kt_power_inst,kt_power_change
!@var name_tconsrv,lname_tconsrv,units_tconsrv: for tracer conservation
      character(len=20), dimension(ktcon,ntm) ::
     &   name_tconsrv,units_tconsrv
      character(len=80), dimension(ktcon,ntm) :: lname_tconsrv
!@var SCALE_INST,SCALE_CHANGE: Scale factors for tracer conservation
      REAL*8, dimension(ntm) :: SCALE_INST,SCALE_CHANGE
!@var itcon_surf Index array for surface source/sink conservation diags
      INTEGER, DIMENSION(ntsurfsrcmax,NTM) :: itcon_surf
!@var itcon_3D Index array for 3D source/sink conservation diags
      INTEGER, DIMENSION(nt3Dsrcmax,NTM) :: itcon_3Dsrc
!@var itcon_decay Index array for decay conservation diags
      INTEGER, DIMENSION(NTM) :: itcon_decay
!@var itcon_grav Index array for gravitational settling conserv. diags
      INTEGER, DIMENSION(NTM) :: itcon_grav
!@var itcon_mc Index array for moist convection conserv. diags
      INTEGER, DIMENSION(NTM) :: itcon_mc
!@var itcon_grav Index array for large-scale condensation conserv. diags
      INTEGER, DIMENSION(NTM) :: itcon_ss

C----------------------------------------------------
!@param KTACC total number of tracer diagnostic words
!@var TACC: Contains all tracer diagnostic accumulations
      INTEGER, PARAMETER ::
     *  ktacc=IM*JM*LM*NTM + IM*JM*ktaij*NTM + IM*JM*ktaijs +
     *        JM*LM*ktajlx*NTM + JM*LM*ktajls + JM*NTM*ktcon
      COMMON /TACCUM/ TAIJLN,TAIJN,TAIJS,TAJLN,TAJLS,TCONSRV
      REAL*8, DIMENSION(KTACC) :: TACC
      EQUIVALENCE (TACC,TAIJLN)
C----------------------------------------------------
#endif
      END MODULE TRACER_DIAG_COM

#ifdef TRACERS_ON
      SUBROUTINE SET_TCON(QCON,NAME_CON,QSUM,INST_UNIT,SUM_UNIT
     *     ,INST_SC,CHNG_SC, itr,CONPTs)
!@sum  SET_TCON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc,nfiltr
      USE DAGCOM, only: npts,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_src,conpt0
      USE TRACER_DIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
     *     ,nofmt,ia_tcon,name_tconsrv,lname_tconsrv,units_tconsrv
     *     ,ntcons
      IMPLICIT NONE
!@var QCON denotes at which points conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QCON
!@var QSUM sets whether each diag is included in final sum
!@+   should be zero for diags set using DIAGTCB (i.e. using difference)
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QSUM
      LOGICAL, DIMENSION(ktcon-1) :: QSUM_CON   ! local version
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*20, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*20, INTENT(IN) :: SUM_UNIT
!@var ITR index for the tracer
      INTEGER, INTENT(IN) :: ITR
!@var CONPTS conservation diag points for special tracer diags
      CHARACTER*16, INTENT(IN), DIMENSION(ntcons) :: CONPTS
      CHARACTER*11 CHGSTR
      INTEGER NI,NM,NS,N,k

C**** remove spaces in NAME_CON for netcdf names
      sname=TRIM(NAME_CON)
      do k=1,len_trim(NAME_CON)
        if (sname(k:k).eq." ") sname(k:k)="_"
      end do
C****
      NI=1
      NOFMT(1,itr) = NI
      TITLE_TCON(NI,itr) =
     *     "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_TCON(NI,itr) = INST_SC
      name_tconsrv(NI,itr) ="inst_"//sname
      lname_tconsrv(NI,itr) = "INSTANT "//TRIM(NAME_CON)
      units_tconsrv(NI,itr) = INST_UNIT
      NSUM_TCON(NI,itr) = -1
      IA_TCON(NI,itr) = 12
      NM=NI
      DO N=2,ktcon-1
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFMT(N,itr) = NM
          QSUM_CON(NM)=.FALSE.
          IF (QSUM(N)) QSUM_CON(NM)=.TRUE.
          CHGSTR=" CHANGE OF "
          if (n.le.npts+1) then
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *         CONPT0(N-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPT0(N-1)(1:3))
          else
            IF (.not. QSUM(N)) CHGSTR="     DELTA "
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *           CONPTs(N-npts-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPTs(N-npts-1)(1:3))
          end if
          lname_tconsrv(NM,itr) = TITLE_TCON(NM,itr)
          units_tconsrv(NM,itr) = SUM_UNIT
          SELECT CASE (N)
          CASE (2)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5d
          CASE (3,4,5,6,7,9,11,12)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5s
          CASE (8)
            SCALE_TCON(NM,itr) = CHNG_SC/(NFILTR*DTSRC)
            IA_TCON(NM,itr) = ia_filt
          CASE (10)
            SCALE_TCON(NM,itr) = CHNG_SC*2./SDAY
            IA_TCON(NM,itr) = ia_12hr
          CASE (13:)   ! special tracer sources
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_src
          END SELECT
        ELSE
          NOFMT(N,itr) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KTCON) THEN
        WRITE(6,*) "KTCON not large enough for extra conserv diags",
     *       KTCON,NI,NM,NS,NAME_CON
        call stop_model(
     &       "Change KTCON in tracer diagnostic common block",255)
      END IF
      DO NM=NI+1,NS-1
        NSUM_TCON(NM,itr) = -1
        IF (QSUM_CON(NM)) NSUM_TCON(NM,itr) = NS
      END DO
      TITLE_TCON(NS,itr) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_Tconsrv(NS,itr) ="sum_chg_"//trim(sname)
      lname_Tconsrv(NS,itr) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_Tconsrv(NS,itr) = SUM_UNIT
      SCALE_TCON(NS,itr) = 1.
      IA_TCON(NS,itr) = 12
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcon

#endif
