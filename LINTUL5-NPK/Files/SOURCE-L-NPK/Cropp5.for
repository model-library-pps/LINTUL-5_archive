*------------------------------------------------------------------------- *
*  SUBROUTINE CROPP                                                        *
*  Author: Joost Wolf                                                      *
*  Date: Crop model developed on the basis of LINTUL3.fst in August 2011   * 
*  Purpose: This subroutine simulates the dry matter increase of a         *
*           crop as function of intercepted radiation, temperature,        *
*           radiation use efficiency, water and N, P and K availability.   *
*                                                                          *
*                                                                          *
*  FORMAL PARAMETERS:(I= input, O= output, C= control, IN= init., T-time)  *
*  name     meaning                                  units       class     *
*  ----     -------                                  -----       -----     *
*  ICROP    number of crop data set                   -           I        *
*  INITI    indicates initialization of run           -           I,C      *
*  IOPT     indicates optimal (=1), water limited (=2),                    *
*           water and N limited (=3) and				    				 *
*           water and N, P and K limited run (=4)     -           I,C      *
*  IDAY     julian day number                         -           I        *
*  IDEM     date of emergence                         -           I        *
*  IDEMERG  date of emergence                         -           O        *
*  IDPL     date of planting/sowing                   -           I        *
*  IDFLOW   date of flowering                                     O        *
*  IDHALT   end date of crop growth                               O        *
*  PL       indicates planting as start of simulation -           I,C      *
*  TERMIN   indicates terminal section                -           I,C      *
*  EMERG    indicates crop emergence                  -           I,C      *
*  TMIN     minimum air temperature                   C           I        *
*  TMAX     maximum air temperature                   C           I        *
*  AVRAD    daily total irradiation                   J m-2 d-1   I        *
*  CO       atmospheric CO2 concentration             ppmv        I        *
*  TRANRF   reduction factor due to drought/wetness   -           I        *  
*  RDMSO    soil related maxiumum rooting depth       cm          I        *
*  DAYLP    photoperiodically active daylength        h           I        *
*  TAGB     total above-ground biomass                kg DM ha-1  O        *
*  WLVG      weight of living leaves                  kg DM ha-1  O        *
*  WLVD     weight of dead leaves                     kg DM ha-1  O        *
*  WST      weight of stems                           kg DM ha-1  O        *
*  WRT      weight of roots                           kg DM ha-1  O        *
*  WSO      weight of storage organs                  kg DM ha-1  O        *
*  RD       actual rooting depth                      cm          O        *
*  RDMCR    crop specific maximum rooting depth       cm          O        *
*  RR       root growth rate                          cm d-1      O        *
*  RDM      soil/crop related maximal rooting depth   cm          I        *
*  LAI      leaf area index                           m2 m-2      O        *
*  DEPNR   crop group number for soil water depletion -           O        *
*  CFET    crop specific correction for transpiration -           O        * 
*  IAIRDU  air ducts in roots present (=1) or not(=0)  -          O,C      *
*  TSUM     temperature sum from emergence            C d         O        *
*  DVS      development stage                         -           O        *
*  DVSEND   development stage at end of growth period -           O        *
*  TSULP    temperature sum from sowing/planting      C d         O        *
*  FINT     fractional light interception for PAR     -           O        *
*  FINTT    fractional light interception for total radiation -   O        *
*  TPARINT  total intercepted radiation (PAR)         MJ m-2      O        *
*  TPAR     total photosynthetically active radiation MJ m-2      O        *
*  TSUML    temperature sum from emergence incl. dayl.effect C.d  O        *
*  NNI      nitrogen nutrition index                  -           O        *
*  NPKI     NPK nutrition index (=minimum of N/P/K-index) -       O        *
*  NMINT    total mineral N from soil and fertiliser  kg N ha-1   O        *
*  NMIN     mineral N available from soil for crop    kg N ha-1   O        *
*  NUPTT    total N uptake by crop from soil          kg N ha-1   O        *
*  NFIXTT   total N uptake by crop from biol.fixation kg N ha-1   O        *
*  NLIVT    Amount of N in living crop organs         kg N ha-1   O        *
*  NLOSST   Amount of N in dead crop organs           kg N ha-1   O        *
*  PMINT    total mineral P from soil and fertiliser  kg P ha-1   O        *
*  PMIN     mineral P available from soil for crop    kg P ha-1   O        *
*  PUPTT    total P uptake by crop from soil          kg P ha-1   O        *
*  PLIVT    Amount of P in living crop organs         kg P ha-1   O        *
*  PLOSST   Amount of P in dead crop organs           kg P ha-1   O        *
*  KMINT    total mineral K from soil and fertiliser  kg K ha-1   O        *
*  KMIN     mineral K available from soil for crop    kg K ha-1   O        *
*  KUPTT    total K uptake by crop from soil          kg K ha-1   O        *
*  KLIVT    Amount of K in living crop organs         kg K ha-1   O        *
*  KLOSST   Amount of K in dead crop organs           kg K ha-1   O        *
*  YCH      indicates year change                     -           I        * 
*                                                                          *
*------------------------------------------------------------------------
       SUBROUTINE CROPP(ICROP,INITI,IOPT, IDAY,IDEM,IDEMERG,IDPL,IDFLOW,
     $    IDHALT,PL,TERMIN,EMERG,TMIN,TMAX,AVRAD,CO,TRANRF,RDMSO,DAYLP,
     $	    TAGB,WLVG, WLVD, WST,WRT,WSO,RD,RDMCR,RR,RDM,LAI,
     $      CFET,DEPNR,IAIRDU,TSUM,DVS,DVSEND,TSULP,FINT,FINTT,TPARINT,
     $        TPAR,TSUML,NNI,NPKI,NMINT,NMIN,NUPTT,NFIXTT,NLIVT,NLOSST,
     $        PMINT,PMIN,PUPTT,PLIVT,PLOSST,
     $        KMINT,KMIN,KUPTT,KLIVT,KLOSST,YCH)

       IMPLICIT REAL (A-Z)
       INTEGER ICROP,IDAY,IDEM,IDEMERG,IDPL,IDFLOW,IDHALT,IDSL
       INTEGER ILCO,IOPT,IAIRDU,ILDTSM,ILSLA,ILSSA,ILKDIF,ILRUE
	 INTEGER ILTMPF,ILTMNF,ILFR,ILFL,ILFS,ILFO,ILRDRL
	 INTEGER ILRDRR,ILRDRS,ILNMXL,ILPHOT,ILFERN,ILNRFT 
	 INTEGER ILPMXL,ILKMXL, ILFERP, ILFERK,ILPRFT,ILKRFT

       LOGICAL TERMIN,INITI,PL,FLOW,EMERG,YCH
       REAL DTSMTB (30), SLATB (30),SSATB (30),KDIFTB (30)  
       REAL TMPFTB (30), TMNFTB (30), RUETB (30), FLTB (30), RDRLTB (30) 
       REAL COTB (30), RDRRTB (30), RDRSTB (30), NMXLV (30), PHOTTB (30)
	 REAL PMXLV (30), KMXLV (30)
	 REAL FSTB (30), FOTB (30), FERNTAB(30), NRFTAB(30),PRFTAB(30)
	 REAL FRTB (30), FERPTAB(30), FERKTAB(30),KRFTAB(30) 

       COMMON /CROPOUT/ TAGB1,WSO1,TPARINT1, HI1, RUEC
	 COMMON /CROPOUT/ NUPTT1,NFIXTT1,NLIV1,NLOSS1,NROOT1,TRANRF1,NNI1
       COMMON /CROPOUT/ PUPTT1,PLIV1,PLOSS1,PROOT1,NPKI1
       COMMON /CROPOUT/ KUPTT1,KLIV1,KLOSS1,KROOT1
*----------------------------------------------------
*      INITIALIZATION
*----------------------------------------------------

       IF (INITI) THEN
         FLOW= .FALSE.
         EMERG= .FALSE.
       ELSE
         GOTO 200
       ENDIF

*---- Reading crop data from file CROPP.DAT
       
       IF (ICROP .EQ. 1) CALL RDINIT (14,0, 'CROPP1.DAT')
       IF (ICROP .EQ. 2) CALL RDINIT (14,0, 'CROPP2.DAT')
       IF (ICROP .EQ. 3) CALL RDINIT (14,0, 'CROPP3.DAT')
       IF (ICROP .EQ. 4) CALL RDINIT (14,0, 'CROPP4.DAT')
       IF (ICROP .EQ. 5) CALL RDINIT (14,0, 'CROPP5.DAT')
       IF (ICROP .EQ. 6) CALL RDINIT (14,0, 'CROPP6.DAT')

      
!    Read initial states and parameter values
	call RDSINT ('IDSL', IDSL)
      CALL RDSREA ('TSUM1',TSUM1)
      CALL RDSREA ('TSUM2',TSUM2)

      CALL RDSREA ('DVSI', DVSI)
      CALL RDSREA ('DVSEND',DVSEND)
      CALL RDSREA ('TDWI',TDWI)
      CALL RDSREA ('RGRLAI',RGRLAI)
      
      CALL RDSREA ('SPA',SPA)
      CALL RDAREA ('DTSMTB',DTSMTB,30,ILDTSM)      
      CALL RDAREA ('SLATB',SLATB,30,ILSLA)
      CALL RDAREA ('SSATB',SSATB,30,ILSSA)
      CALL RDSREA ('TBASE',TBASE)
	CALL RDAREA ('KDIFTB',KDIFTB,30,ILKDIF)
      CALL RDAREA ('RUETB',RUETB,30,ILRUE)
      CALL RDAREA ('TMPFTB',TMPFTB,30,ILTMPF)
      CALL RDAREA ('TMNFTB',TMNFTB,30,ILTMNF)
      CALL RDAREA ('COTB',COTB,30,ILCO)
      CALL RDAREA ('FRTB',FRTB,30,ILFR)
      CALL RDAREA ('FLTB',FLTB,30,ILFL)
      CALL RDAREA ('FSTB',FSTB,30,ILFS)
      CALL RDAREA ('FOTB',FOTB,30,ILFO)
      CALL RDSREA ('RDRL',RDRL)
      CALL RDAREA ('RDRLTB',RDRLTB,30,ILRDRL)
      call RDSREA ('RDRSHM', RDRSHM)
      call RDSREA ('RDRNS', RDRNS)
      CALL RDAREA ('RDRRTB',RDRRTB,30,ILRDRR)
      CALL RDAREA ('RDRSTB',RDRSTB,30,ILRDRS)
      CALL RDSREA ('CFET',CFET)
      CALL RDSREA ('DEPNR',DEPNR)
      CALL RDSINT ('IAIRDU',IAIRDU)
      CALL RDSREA ('RDI',RDI)
      CALL RDSREA ('RRI',RRI)
      CALL RDSREA ('RDMCR',RDMCR)
      CALL RDSREA ('DVSDR', DVSDR)
      CALL RDSREA ('DVSDLT', DVSDLT)
      CALL RDSREA ('DVSNLT', DVSNLT)
      CALL RDSREA ('DVSNT', DVSNT)
      CALL RDSREA ('TBASEM',TBASEM)
      CALL RDSREA ('TEFFMX',TEFFMX)
      CALL RDSREA ('TSUMEM',TSUMEM)
      call RDSREA ('FNTRT', FNTRT)
      call RDSREA ('FRNX', FRNX)
      call RDSREA ('FRPX', FRPX)
      call RDSREA ('FRKX', FRKX)
      call RDSREA ('LAICR', LAICR)
      call RDSREA ('LRNR', LRNR)
      call RDSREA ('LSNR', LSNR)
      call RDSREA ('LRPR', LRPR)
      call RDSREA ('LSPR', LSPR)
      call RDSREA ('LRKR', LRKR)
      call RDSREA ('LSKR', LSKR)
      call RDSREA ('NLAI', NLAI)
      call RDSREA ('NLUE', NLUE)
      call RDSREA ('NMAXSO', NMAXSO)
      call RDSREA ('PMAXSO', PMAXSO)
      call RDSREA ('KMAXSO', KMAXSO)
      call RDSREA ('NPART', NPART)
      call RDSREA ('NFIXF', NFIXF)
      call RDSREA ('NSLA', NSLA)
      call RDSREA ('RNFLV', RNFLV)
      call RDSREA ('RNFRT', RNFRT)
      call RDSREA ('RNFST', RNFST)
      call RDSREA ('TCNT', TCNT) 
      call RDAREA ('NMXLV', NMXLV, 30, ILNMXL)
      call RDSREA ('RPFLV', RPFLV)
      call RDSREA ('RPFRT', RPFRT)
      call RDSREA ('RPFST', RPFST)
      call RDSREA ('TCPT', TCPT) 
      call RDAREA ('PMXLV', PMXLV, 30, ILPMXL)
      call RDSREA ('RKFLV', RKFLV)
      call RDSREA ('RKFRT', RKFRT)
      call RDSREA ('RKFST', RKFST)
      call RDSREA ('TCKT', TCKT) 
      call RDAREA ('KMXLV', KMXLV, 30, ILKMXL)
      call RDAREA ('PHOTTB', PHOTTB, 30, ILPHOT)

      
      CLOSE (14, STATUS= 'DELETE')

*     Read management data from file MANAGE.DAT
      Call RDINIT (15,0, 'MANAGE.DAT')
      call RDAREA ('FERNTAB', FERNTAB, 30, ILFERN) 
      call RDAREA ('NRFTAB', NRFTAB, 30, ILNRFT)        
	call RDSREA ('NMINS', NMINS)
	call RDSREA  ('RTNMINS', RTNMINS)
      call RDAREA ('FERPTAB', FERPTAB, 30, ILFERP) 
      call RDAREA ('PRFTAB', PRFTAB, 30, ILPRFT)        
	call RDSREA ('PMINS', PMINS)
	call RDSREA  ('RTPMINS', RTPMINS)
      call RDAREA ('FERKTAB', FERKTAB, 30, ILFERK) 
      call RDAREA ('KRFTAB', KRFTAB, 30, ILKRFT)        
	call RDSREA ('KMINS', KMINS)
	call RDSREA  ('RTKMINS', RTKMINS)
      CLOSE (15, STATUS= 'DELETE')


*----- initialization of state variables

       TDW= TDWI
	 DVS= DVSI
	 DVR= 0.
	 NMINI=NMINS
	 NMIN=NMINI
	 PMINI=PMINS
	 PMIN=PMINI
	 KMINI=KMINS
	 KMIN=KMINI
	 WRTI=  LINT (FRTB, ILFR, DVSI) * TDW
	 WRT= WRTI
	 TAGB= TDW- WRT 
	 WLVGI= LINT (FLTB, ILFL, DVSI) * TAGB
	 WLVG= WLVGI
	 LAII= WLVGI * LINT(SLATB,ILSLA,DVSI) 
	 LAI= LAII
	 RLAI= 0.
	 WLVD= 0.
	 WSTI=  LINT (FSTB, ILFS, DVSI) * TAGB
	 WST= WSTI 
       WSOI=  LINT (FOTB, ILFO, DVSI) * TAGB
	 WSO= WSOI
	 WSTD= 0.
	 WRTD= 0.
	 RWLVG= 0.
	 RWST=0.
	 RWRT=0.
	 RWSO= 0.
	 DLV=0.
	 DRST=0.
	 DRRT= 0.
       GRT= 0.
       TSULP= 0. 
       TSUM= 0.
	 TSUML= 0.
       DTSULP= 0. 
       DTSUM= 0.
	 DTSUML= 0.
	 DVRED= 1.
       TPAR= 0.
       TPARINT= 0.
       PAR= 0.
       PARINT= 0. 
	 ATN= 0. 
	 ATP= 0.
	 ATK= 0.
	 GTSUM= 0.  
	 RD= RDI 
	 DELT= 1. 
	 NMINT= 0.
	 PMINT= 0.
	 KMINT= 0.
	 NMAXLVI= LINT (NMXLV, ILNMXL, DVSI)
       NMAXSTI= LSNR * NMAXLVI
	 NMAXRTI= LRNR * NMAXLVI
       ANLVI= NMAXLVI * WLVGI
	 ANSTI= NMAXSTI * WSTI
	 ANRTI= NMAXRTI * WRTI
	 ANLV= ANLVI
	 ANST= ANSTI
	 ANRT= ANRTI
	 ANSOI= 0.
	 ANSO= ANSOI
	 NLOSSL= 0.
	 NLOSSR= 0.
	 NLOSSS= 0.
	 NUPTT= 0.
	 NFIXTT= 0.
	 RNMINS= 0.
	 RNMINT= 0.
	 PMAXLVI= LINT (PMXLV, ILPMXL, DVSI)
       PMAXSTI= LSPR * PMAXLVI
	 PMAXRTI= LRPR * PMAXLVI
       APLVI= PMAXLVI * WLVGI
	 APSTI= PMAXSTI * WSTI
	 APRTI= PMAXRTI * WRTI
	 APLV= APLVI
	 APST= APSTI
	 APRT= APRTI
	 APSOI= 0.
	 APSO= APSOI
	 PLOSSL= 0.
	 PLOSSR= 0.
	 PLOSSS= 0.
	 PUPTT= 0.
	 RPMINS= 0.
	 RPMINT= 0.
	 KMAXLVI= LINT (KMXLV, ILKMXL, DVSI)
       KMAXSTI= LSKR * KMAXLVI
	 KMAXRTI= LRKR * KMAXLVI
       AKLVI= KMAXLVI * WLVGI
	 AKSTI= KMAXSTI * WSTI
	 AKRTI= KMAXRTI * WRTI
	 AKLV= AKLVI
	 AKST= AKSTI
	 AKRT= AKRTI
	 AKSOI= 0.
	 AKSO= AKSOI
	 KLOSSL= 0.
	 KLOSSR= 0.
	 KLOSSS= 0.
	 KUPTT= 0.
	 RKMINS= 0.
	 RKMINT= 0.
	 CTRAN= 0.
       CNNI= 0.
	 CPNI= 0.
	 CKNI= 0.
	 CNPKI= 0.

       DAYPL= REAL(IDPL)
	 TRANRF= 1.
	 NNI= 1.
       
       IF (.NOT. PL) THEN
          DAYEM= REAL(IDEM)
	    EMERG= .TRUE.
       ELSE
          DAYEM= 999.
       ENDIF


200     CONTINUE

        


        IF (TERMIN) GOTO 999


*---------------------------------------------------------
*      INTEGRATION
*---------------------------------------------------------

*----- Temperature sums (C.d) from sowing/planting (P) and from emergence with and without daylength effect 

       TSULP= INTGRL(TSULP, DTSULP, 1.)
       TSUM= INTGRL(TSUM, DTSUM, 1.)
	 TSUML= INTGRL(TSUML, DTSUML, 1.)


*----- Start of flowering
       IF (.NOT. FLOW .AND. TSUML .GE. TSUM1)  IDFLOW= IDAY
       IF (.NOT. FLOW .AND. TSUML .GE. TSUM1)  FLOW= .TRUE.

*---- Total photosynthetically active radiation (MJ/m2)
      TPAR= INTGRL(TPAR, PAR, 1.)

*---- Total intercepted radiation (MJ/m2)
      TPARINT= INTGRL(TPARINT, PARINT, 1.)

*---- Dry weights of total biomass, living crop organs, and total above-ground living biomass (kg DM/ha)
	GTSUM= INTGRL(GTSUM, GRT, 1.)
      WLVG= INTGRL(WLVG, RWLVG, 1.) 
	WST=  INTGRL(WST, RWST, 1.)
	WRT=  INTGRL(WRT,RWRT, 1.)
	WSO=  INTGRL(WSO, RWSO, 1.)
	TAGBG= WLVG+WST+WSO

*---- Development stage (-)
	DVS=  INTGRL(DVS, DVR, 1.)
      
*---- Dry weights of dead crop organs and total above-ground biomass incl. dead crop organs (kg DM/ha)
      WLVD= INTGRL(WLVD, DLV, 1.)
	WSTD= INTGRL(WSTD,DRST, 1.)
	WRTD= INTGRL(WRTD, DRRT, 1.)
	TAGB= TAGBG + WLVD + WSTD


*----- Rooting depth and Leaf area index
       RD= INTGRL(RD, RR, 1.)
	 LAI= INTGRL(LAI, RLAI, 1.)

*     Soil mineral N/P/K  and Total mineral N/P/K available from both fertiliser and soil (kg N/P/K ha-1)
       NMIN= INTGRL(NMIN, RNMINS, 1.)
	 NMINT= INTGRL(NMINT, RNMINT,1.) 
	 PMIN= INTGRL(PMIN, RPMINS, 1.)
	 PMINT= INTGRL(PMINT, RPMINT,1.) 
	 KMIN= INTGRL(KMIN, RKMINS, 1.)
	 KMINT= INTGRL(KMINT, RKMINT,1.) 


*----- Total N/P/K uptake by crop over time (kg N/P/K ha-1) from soil and by biological fixation
       NUPTT= INTGRL(NUPTT, NUPTR, 1.)
       PUPTT= INTGRL(PUPTT, PUPTR, 1.)
       KUPTT= INTGRL(KUPTT, KUPTR, 1.)
	 NFIXTT= INTGRL(NFIXTT, NFIXTR, 1.)

*-----Actual N/P/K amount in various living organs and total living N/P/K amount(kg N/P/K ha-1)
      ANLV =  INTGRL (ANLV,RNLV, 1.)
      ANST =  INTGRL (ANST,RNST, 1.)
      ANRT =  INTGRL (ANRT,RNRT, 1.)
      ANSO =  INTGRL (ANSO,RNSO, 1.)
	NLIVT= ANLV+ANST+ANRT+ANSO

	APLV =  INTGRL (APLV,RPLV, 1.)
      APST =  INTGRL (APST,RPST, 1.)
      APRT =  INTGRL (APRT,RPRT, 1.)
      APSO =  INTGRL (APSO,RPSO, 1.)
	PLIVT= APLV+APST+APRT+APSO

	AKLV =  INTGRL (AKLV,RKLV, 1.)
      AKST =  INTGRL (AKST,RKST, 1.)
      AKRT =  INTGRL (AKRT,RKRT, 1.)
      AKSO =  INTGRL (AKSO,RKSO, 1.)
	KLIVT= AKLV+AKST+AKRT+AKSO

*-----N/P/K losses from leaves, roots and stems due to senescence and total N/P/K loss (kg N/P/K ha-1)
      NLOSSL=  INTGRL(NLOSSL, RNLDLV, 1.)
	NLOSSR=  INTGRL(NLOSSR, RNLDRT, 1.)
	NLOSSS=  INTGRL(NLOSSS, RNLDST, 1.)
      NLOSST=  NLOSSL+NLOSSR+NLOSSS

	PLOSSL=  INTGRL(PLOSSL, RPLDLV, 1.)
	PLOSSR=  INTGRL(PLOSSR, RPLDRT, 1.)
	PLOSSS=  INTGRL(PLOSSS, RPLDST, 1.)
      PLOSST=  PLOSSL+PLOSSR+PLOSSS

	KLOSSL=  INTGRL(KLOSSL, RKLDLV, 1.)
	KLOSSR=  INTGRL(KLOSSR, RKLDRT, 1.)
	KLOSSS=  INTGRL(KLOSSS, RKLDST, 1.)
      KLOSST=  KLOSSL+KLOSSR+KLOSSS

*---- total N/P/K in living and dead roots
      NROOT= ANRT + NLOSSR
	PROOT= APRT + PLOSSR
	KROOT= AKRT + KLOSSR
*-----------------------------------------------------------------------------
*      RATE CALCULATIONS
*-----------------------------------------------------------------------------

*------ Weather calculations

*------ Daily photosynthetically active radiation (PAR, MJ/m2)
        PAR= AVRAD/1.0E6 * 0.50

*------ Average daily temperature (C)
        TMPA= 0.5 *(TMIN + TMAX)
        DAY= REAL(IDAY)

*------ date of planting/sowing
        IF (PL) THEN
          IF (.NOT. YCH) PUSHPL= INSW(DAY - DAYPL, 0., 1.)
	    IF (YCH) PUSHPL= INSW(365. + DAY - DAYPL, 0., 1.)
        ELSE
          PUSHPL= 0.
        ENDIF

        IF (PL .AND. (TSULP .GE. TSUMEM) .AND. (.NOT. EMERG)) THEN
           IDEMERG= IDAY
           DAYEM= REAL(IDEMERG)
           EMERG= .TRUE.
        ENDIF 

*-----  Reduction of development rate until flowering by  day length

        DVRED= LINT(PHOTTB, ILPHOT, DAYLP)

        IF (IDSL .EQ. 1 .AND. .NOT. FLOW) THEN
           RDAYL= DVRED
        ELSE
           RDAYL= 1.
        ENDIF

*       emergence date set
        IF (.NOT. PL) IDEMERG= IDEM

*------ Change in temperature sums from sowing/planting (P) and from emergence for crop development without and with day length effect
        DTSULP= LIMIT(0., TEFFMX-TBASEM, TMPA-TBASEM) * PUSHPL
        IF (.NOT. YCH) PUSHEM= INSW(DAY - DAYEM, 0., 1.)
        IF (YCH) PUSHEM=INSW(365. + DAY-DAYEM, 0., 1.)
        DTSU=  MAX(0.,LINT(DTSMTB, ILDTSM, TMPA))
	  DTSUM= DTSU * PUSHEM
        DTSUML= DTSU * PUSHEM * RDAYL

*       Calculation of development stage        
        IF (DVS .LT. 1.0) THEN
*       effects of daylength and tenmperature on development during vegetative phase
	    DVR= DTSUML/TSUM1 
	  ELSE
*       development during generative phase
          DVR= DTSUML/TSUM2
	  END IF
		   

*------ Plant growth
        
*------ Radiation use efficiency as dependent on development stage (g DM MJ-1)
        RUE= LINT(RUETB,ILRUE,DVS)
        
*------ Correction of radiation use efficiency for change in atmospheric CO2 concentration (-)
        RCO= LINT(COTB,ILCO,CO)
        
*------ Reduction of radiation use efficiency for non-optimal day-time temperatures and for low minimum temperature
        DTEMP= TMAX - 0.25*(TMAX-TMIN)
        RTMP= LINT(TMPFTB,ILTMPF,DTEMP) * LINT(TMNFTB,ILTMNF,TMIN)
*------ Correction of RUE for both non-optimal temperatures and atmospheric CO2
	  RTMCO= RTMP * RCO
        
!      Calling the subroutine for translocatable N/P/K in leaves, stem, roots and
!      storage organs (kg N/P/K ha-1)

       CALL NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT,RNFLV,RNFST,RNFRT,
     $                  FNTRT,ATNLV,ATNST,ATNRT,ATN, 
     $                  APLV,APST,APRT,AKLV,AKST,AKRT, 
     $                  RPFLV,RPFST,RPFRT,RKFLV,RKFST,RKFRT, 
     $                  ATPLV,ATPST,ATPRT,ATP,ATKLV,ATKST,ATKRT,
     $                  ATK)
     
!    * Total vegetative living above-ground biomass (kg DM ha-1)
        TBGMR =WLVG+WST
    
!      N/P/K concentrations (kg N/P/K kg-1 DM) in the living leaves, stem, roots and storage
!      organs
       NFLV = ANLV/ NOTNUL(WLVG)
       NFST = ANST/ NOTNUL(WST)
       NFRT = ANRT/ NOTNUL(WRT)
       NFSO = ANSO/ NOTNUL(WSO)
    
       PFLV = APLV/ NOTNUL(WLVG)
       PFST = APST/ NOTNUL(WST)
       PFRT = APRT/ NOTNUL(WRT)
       PFSO = APSO/ NOTNUL(WSO)

       KFLV = AKLV/ NOTNUL(WLVG)
       KFST = AKST/ NOTNUL(WST)
       KFRT = AKRT/ NOTNUL(WRT)
       KFSO = AKSO/ NOTNUL(WSO)

!    Total N/P/K in vegetative living above-ground biomass (kg N/P/K ha-1)
       NUPGMR = ANLV + ANST
	 PUPGMR = APLV + APST
	 KUPGMR = AKLV + AKST
     
!    Fertilizer N/P/K application (kg N/P/K ha-1 d-1) and its recovery fraction (-)--------------*
       FERTN  = LINT (FERNTAB,ILFERN, DAY)
       NRF    = LINT (NRFTAB,ILNRFT, DAY)
       FERTNS = FERTN * NRF
       FERTP  = LINT (FERPTAB,ILFERP, DAY)
       PRF    = LINT (PRFTAB,ILPRFT, DAY)
       FERTPS = FERTP * PRF
       FERTK  = LINT (FERKTAB,ILFERK, DAY)
       KRF    = LINT (KRFTAB,ILKRFT, DAY)
       FERTKS = FERTK * KRF

!      Check on N balance
       NBALAN = ABS(NUPTT+NFIXTT+(ANLVI+ANSTI+ANRTI+ANSOI)-(ANLV
     $          +ANST+ANRT+ANSO+NLOSSL+NLOSSR+NLOSSS))  

      IF (NBALAN .GE. 1.) STOP
     $    ' nitrogen balance NBALAN not 0, program aborted'
   
!      Check on P balance
       PBALAN = ABS(PUPTT+(APLVI+APSTI+APRTI+APSOI)-(APLV
     $          +APST+APRT+APSO+PLOSSL+PLOSSR+PLOSSS))  

      IF (PBALAN .GE. 1.) STOP
     $    ' phosphorus balance PBALAN not 0, program aborted'

!      Check on K balance
       KBALAN = ABS(KUPTT+(AKLVI+AKSTI+AKRTI+AKSOI)-(AKLV
     $          +AKST+AKRT+AKSO+KLOSSL+KLOSSR+KLOSSS))  

      IF (KBALAN .GE. 1.) STOP
     $    ' potassium balance KBALAN not 0, program aborted'

!    Total leaf weight, both green and dead (kg DM ha-1)
       WLV    = WLVG + WLVD
  
!     Relative death rate of roots (d-1)
       RDRRT = LINT(RDRRTB, ILRDRR, DVS)
	 RDRST = LINT(RDRSTB, ILRDRS, DVS)
   
!     Total N/P/K in living above-ground crop organs (kg N/P/K ha-1)   
       NTAG   = ANLV   + ANST   + ANSO
	 PTAG   = APLV   + APST   + APSO
	 KTAG   = AKLV   + AKST   + AKSO
   
!     Relative death rate of leaves due to senescence/ageing as dependent on mean daily temperature (d-1)
       RDRTMP = LINT(RDRLTB,ILRDRL,TMPA)
  
!     Maximum N/P/K concentration in the leaves, from which the N/P/K conc. in the
!     stem and roots are derived, as a function of development stage (kg N/P/K kg-1 DM)  
       NMAXLV = LINT (NMXLV, ILNMXL, DVS)
       PMAXLV = LINT (PMXLV, ILPMXL, DVS)
       KMAXLV = LINT (KMXLV, ILKMXL, DVS)

!     N/P/K concentration in above-ground living biomass (kg N/P/K kg-1 DM)
	 NTAC  = NTAG/TAGBG
	 PTAC  = PTAG/TAGBG
	 KTAC  = KTAG/TAGBG
    
!     N/P/K supply to the storage organs (kg N/P/K ha-1 d-1)
       NSUPSO = INSW (DVS-DVSNT,0.,ATN/TCNT)
	 PSUPSO = INSW (DVS-DVSNT,0.,ATP/TCPT)
	 KSUPSO = INSW (DVS-DVSNT,0.,ATK/TCKT)

!      N/P/K concentrations in total vegetative living above-ground biomass  (kg N/P/K kg-1 DM) 
       NFGMR  = NUPGMR/NOTNUL(TBGMR)
	 PFGMR  = PUPGMR/NOTNUL(TBGMR)
	 KFGMR  = KUPGMR/NOTNUL(TBGMR)
     
!    * Residual N/P/K concentrations in total vegetative living above-ground biomass  (kg N/P/K kg-1 DM) 
       NRMR   = (WLVG*RNFLV+WST*RNFST)/NOTNUL(TBGMR)
	 PRMR   = (WLVG*RPFLV+WST*RPFST)/NOTNUL(TBGMR)
	 KRMR   = (WLVG*RKFLV+WST*RKFST)/NOTNUL(TBGMR)
  
!     Nutrient uptake limiting factor (-) at low moisture conditions in the
!     rooted soil layer before anthesis. After DVSNLT, there is no
!     nutrient uptake from the soil
       NLIMIT = INSW(DVS-DVSNLT, INSW(TRANRF-0.01,0.,1.) , 0.0)

!     Biomass partitioning functions under non-stressed situations (-)
       FRTWET = LINT(FRTB,ILFR, DVS )
       FLVT   = LINT(FLTB,ILFL, DVS )
       FSTT   = LINT(FSTB,ILFS, DVS)
       FSOT   = LINT(FOTB,ILFO, DVS )
 
*----  Carbon balance check  
       CBALAN = ABS(GTSUM + (WRTI+WLVGI+WSTI+WSOI)-(WLVG+WST+WSO+WRT
     $           +WLVD+WRTD+WSTD))	

       IF (CBALAN .GE. 1.) STOP
     $    ' carbon balance CBALAN not 0, program aborted'
    
*----  Maximum N/P/K concentrations in stems and roots (kg N/P/K kg-1 DM)
       NMAXST = LSNR * NMAXLV
       NMAXRT = LRNR * NMAXLV
	 PMAXST = LSPR * PMAXLV
       PMAXRT = LRPR * PMAXLV
	 KMAXST = LSKR * KMAXLV
       KMAXRT = LRKR * KMAXLV

*----- Root growth (cm d-1)
       IF (EMERG) RR = MIN(RRI * INSW( TRANRF-0.01, 0., 1. ),  RDM-RD)
        
!    Calling the subroutine for calculating optimal N/P/K concentrations in leaves
!    and stems
       CALL NOPTM(FRNX,NMAXLV,NMAXST,NOPTLV,NOPTST,FRPX,PMAXLV,
     $  PMAXST,POPTLV,POPTST,FRKX,KMAXLV,KMAXST,KOPTLV,KOPTST)
    
!     Optimal amount of N/P/K in vegetative above-ground living biomass and its N/P/K concentration
       NOPTS = NOPTST* WST
       NOPTL = NOPTLV* WLVG 
       NOPTMR = (NOPTL+ NOPTS)/NOTNUL(TBGMR)
   
       POPTS = POPTST* WST
       POPTL = POPTLV* WLVG 
       POPTMR = (POPTL+ POPTS)/NOTNUL(TBGMR)

       KOPTS = KOPTST* WST
       KOPTL = KOPTLV* WLVG 
       KOPTMR = (KOPTL+ KOPTS)/NOTNUL(TBGMR)
    
!     Calling the subroutine for calculating the NPK and N/P/K Nutrition Indices (NPKI & NNI etc.)
       CALL NNINDX(DAY,DAYEM,EMERG,NFGMR,NRMR,NOPTMR,NNI,
     $           	PFGMR,PRMR,POPTMR,PNI,KFGMR,KRMR,KOPTMR,KNI,NPKI)

!     With potential and water limited conditions there is no N/P/K stress and NNI etc. and NPKI are set to 1   
       IF (IOPT .EQ. 1 .or. IOPT .EQ. 2) THEN
	    NNI= 1.
	    PNI= 1.
	    KNI= 1.
	    NPKI= 1.
       ELSE IF (IOPT .EQ. 3) THEN
	    PNI= 1.
	    KNI= 1.
	 ENDIF


!    Calling the subroutine for calculating the relative modification for root and shoot
!    allocation.
       CALL SUBPAR (NPART,TRANRF,NNI,FRTWET,FLVT,FSTT,FSOT,
     $              FSHMOD,FLVMOD,FRT,FLV,FST,FSO)

       FCHECK  = ABS(FRT + (FLV+ FST+ FSO) * (1.- FRT) -1.)

       IF (FCHECK .GE. 0.05) STOP
     $ ' assimilate allocation check over crop organs FCHECK '
     $'not 0, program aborted'

    
!    Calling the subroutine for total growth rate (kg DM ha-1 d-1)
       KDIF= LINT(KDIFTB,ILKDIF,DVS)
       CALL GROWTH(DAY,EMERG,PAR,KDIF,NLUE,LAI,RUE,RTMCO,TRANRF,FINT,
     $	 FINTT,NPKI,PARINT,GRT)
  
!     Water-Nutrient stress factor
       RNW = MIN(TRANRF,NPKI)
     
!     cumulative values for TRANRF, NNI etc. and NPKI over growth period
       IF (EMERG) CTRAN= CTRAN + TRANRF
	 IF (EMERG) CNNI= CNNI + NNI
	 IF (EMERG) CPNI= CPNI + PNI
	 IF (EMERG) CKNI= CKNI + KNI	
	 IF (EMERG) CNPKI= CNPKI + NPKI

!      Specific Leaf area(ha/kg), as dependent on NPK stress.
       SLA = LINT (SLATB,ILSLA,DVS)*EXP(-NSLA * (1.-NPKI))
    
!    Calling the subroutine for relative death rate of leaves.
       CALL DEATHL(DAY,EMERG,DVS,DVSDLT,RDRTMP,RDRSHM,RDRL,TRANRF,LAI,
     $	 LAICR,WLVG,RDRNS,NPKI,
     $     SLA,RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI)
     
!    ** Leaf growth
       GLV    = FLV * GRT * (1-FRT)
    
!    Calling the subroutine for calculating the daily increase of leaf area index (m2 m-2 d-1).
       DTEFF= MAX(0., TMPA-TBASE)
       CALL GLA(DAY,EMERG,DTEFF,LAII,RGRLAI,DELT,SLA,LAI,GLV,NLAI,
     $      DVS,TRANRF,NPKI,GLAI)
    
!    Calling the subroutine for N/P/K losses due to death of leaves,stems and roots (kg N/P/K ha-1 d-1)
       CALL RNLD(DVS,WRT,WST,RDRRT,RDRST,RNFLV,DLV,RNFRT,
     $ RPFLV,RPFRT,RKFLV,RKFRT,RNFST,RPFST,RKFST,DVSDR,DRRT,DRST,
     $ RNLDLV,RNLDRT,RNLDST,RPLDLV,RPLDRT,RPLDST,RKLDLV,RKLDRT,RKLDST)
  
!    Net rate of change of Leaf area (m2 leaf area m-2 d-1)
       RLAI   = GLAI - DLAI
    
!    Calling the subroutine for calculating the relative growth rate of roots, leaves, stem
!    and storage organs (kg ha-1 d-1)
       CALL RELGR(DAY,DAYEM,EMERG,GRT,FLV,FRT,FST,
     $ FSO,DLV,DRRT,DRST,RWLVG,RWRT,RWST,RWSO)
    
!    Calling the subroutine for N/P/K demand of leaves, roots and stem storage
!    organs (kg N/P/K ha-1 d-1)
       CALL NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG,WST,WRT,WSO,
     $        PMAXLV,PMAXST,PMAXRT,PMAXSO,KMAXLV,KMAXST,KMAXRT,KMAXSO,
     $        ANLV,ANST,ANRT,ANSO,TCNT, APLV,APST,APRT,APSO,TCPT,
     $        AKLV,AKST,AKRT,AKSO,TCKT,
     $        NDEML,NDEMS,NDEMR,NDEMSO,PDEML,PDEMS,PDEMR,PDEMSO,
     $        KDEML,KDEMS,KDEMR,KDEMSO)
     
!        Total N/P/K demand (kg N/P/K ha-1)
       NDEMTO = MAX (0.0,(NDEML + NDEMS + NDEMR))
	 PDEMTO = MAX (0.0,(PDEML + PDEMS + PDEMR))
	 KDEMTO = MAX (0.0,(KDEML + KDEMS + KDEMR))
    
!    Rate of N/P/K uptake in grains (kg N/P/K ha-1 d-1)
       RNSO =  AMIN1 (NDEMSO,NSUPSO)
       RPSO =  AMIN1 (PDEMSO,PSUPSO)
       RKSO =  AMIN1 (KDEMSO,KSUPSO)

	IF (EMERG) THEN
!        Total N/P/K uptake (kg N/P/K ha-1 d-1) from soil and by biological N fixation
         NUPTR = (MAX (0., MIN ((1.-NFIXF)*NDEMTO, NMINT))* NLIMIT)/DELT
	   NFIXTR= MAX (0.,NUPTR * NFIXF / MAX (0.02, 1.-NFIXF) )
!	   Total P/K uptake (kg P/K ha-1 d-1) from soil 
         PUPTR = (MAX (0., MIN (PDEMTO, PMINT))* NLIMIT)/DELT
	   KUPTR = (MAX (0., MIN (KDEMTO, KMINT))* NLIMIT)/DELT

!        No N/P/K limitation for optimal and water limited production
         IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
	     NUPTR=(MAX (0., (1.-NFIXF)*NDEMTO)* NLIMIT ) / DELT 
           NFIXTR=(MAX (0., NFIXF*NDEMTO)* NLIMIT ) / DELT 
	     PUPTR=(MAX (0., PDEMTO) * NLIMIT ) / DELT 
	     KUPTR=(MAX (0., KDEMTO) * NLIMIT ) / DELT 
         ELSEIF (IOPT .EQ. 3) THEN
!        No P/K limitation for nitrogen limited production
	     PUPTR=(MAX (0., PDEMTO) * NLIMIT ) / DELT 
	     KUPTR=(MAX (0., KDEMTO) * NLIMIT ) / DELT 
	   END IF
	ELSE
	  NUPTR= 0.
	  NFIXTR= 0.
	  PUPTR= 0.
	  KUPTR= 0.
	END IF    

!    Calling the subroutine for calculating N/P/K translocated from leaves, stem, and roots (kg N/P/K ha-1 d-1)
       CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT,
     $   RPSO,ATPLV,ATPST,ATPRT,ATP,RPTLV,RPTST,RPTRT,	
     $   RKSO,ATKLV,ATKST,ATKRT,ATK,RKTLV,RKTST,RKTRT)
!    Calling the subroutine to compute the partitioning of the total
!    N/P/K uptake rates (NUPTR,PUPTR,KUPTR) over the leaves, stem and roots (kg N/P/K ha-1 d-1)
       CALL RNUSUB(DAY,DAYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR,
     $          PDEML,PDEMS,PDEMR,PUPTR,KDEML,KDEMS,KDEMR,KUPTR,
     $          NFIXTR,NDEMTO, RNULV,RNUST,RNURT,
     $          PDEMTO, RPULV,RPUST,RPURT,KDEMTO,RKULV,RKUST,RKURT)
!     Soil N/P/K supply (g N m-2 d-1) through mineralization during crop growth
       IF (EMERG) RNMINS= -MAX(0.,MIN( RTNMINS * NMINI * NLIMIT, NMIN))
       IF (EMERG) RPMINS= -MAX(0.,MIN( RTPMINS * PMINI * NLIMIT, PMIN))
       IF (EMERG) RKMINS= -MAX(0.,MIN( RTKMINS * KMINI * NLIMIT, KMIN))

!     Change in total inorganic N/P/K in soil as function of fertilizer
!     input, soil N/P/K mineralization and crop uptake.
       RNMINT = FERTNS/DELT -NUPTR - RNMINS
       RPMINT = FERTPS/DELT -PUPTR - RPMINS
       RKMINT = FERTKS/DELT -KUPTR - RKMINS
     
*----Rate of change of N/P/K in crop organs   
       RNST = RNUST-RNTST-RNLDST
       RNRT = RNURT-RNTRT-RNLDRT
       RNLV = RNULV-RNTLV-RNLDLV   
	 RPST = RPUST-RPTST-RPLDST
       RPRT = RPURT-RPTRT-RPLDRT
       RPLV = RPULV-RPTLV-RPLDLV   
	 RKST = RKUST-RKTST-RKLDST
       RKRT = RKURT-RKTRT-RKLDRT
       RKLV = RKULV-RKTLV-RKLDLV    

999      CONTINUE

* ----- Data for summary output

        IF (TERMIN) THEN
        WRITE (*, '(//2A)') ' Crop development completed'

        TAGB1= TAGB
        WSO1= WSO
        TPARINT1= TPARINT
        HI1= WSO/TAGB
        IDHALT= IDAY
        NUPTT1= NUPTT
	  PUPTT1= PUPTT
	  KUPTT1= KUPTT
	  NFIXTT1= NFIXTT
	  NLIV1= NLIVT
	  NLOSS1= NLOSST
	  NROOT1= NROOT
	  PLIV1= PLIVT
	  PLOSS1= PLOSST
	  PROOT1= PROOT
	  KLIV1= KLIVT
	  KLOSS1= KLOSST
	  KROOT1= KROOT
	    IF (IDAY .GT. IDEMERG) THEN
	    TRANRF1= CTRAN/(IDAY -IDEMERG)
	    NNI1= CNNI/(IDAY-IDEMERG)
	    NPKI1= CNPKI/(IDAY-IDEMERG)
	    ELSE
	    TRANRF1= CTRAN/ (IDAY+365-IDEMERG)
          NNI1= CNNI/(IDAY+365-IDEMERG)
	    NPKI1= CNPKI/(IDAY+365-IDEMERG)
	    END IF

*----- calculated radiation use efficiency in g above-ground D.M./radiation intercepted in MJ PAR
        RUEC= (TAGB/10.)/TPARINT
        END IF


*---- Finish conditions
      IF (FINT .LT. 0.05 .AND.  TAGB .GT. 200.) TERMIN= .TRUE.


         RETURN
         END
