************************************************************************
*
*  fddd.f
*
*  constants for F77-DDD (must correspond to ddd.h in actual version)
*
*  970212 kb  made consistent with ddd-1.8
*
*
************************************************************************


************************************************************************
*
*  constants for DDD_TypeDefine
*
************************************************************************

       INTEGER EL_GDATA, EL_LDATA, EL_GBITS, EL_OBJPTR, EL_CONTINUE, EL_END
       PARAMETER(EL_GDATA=-1)
       PARAMETER(EL_LDATA=-2)
       PARAMETER(EL_GBITS=-3)
       PARAMETER(EL_OBJPTR=-5)
       PARAMETER(EL_CONTINUE=-6)
       PARAMETER(EL_END=-7)


************************************************************************
*
*  constants for DDD-Options 
*
************************************************************************

       INTEGER OPT_IDENTIFY_MODE
       INTEGER OPT_WARNING_VARSIZE_OBJ, OPT_WARNING_SMALLSIZE
       INTEGER OPT_WARNING_PRIOCHANGE, OPT_WARNING_DESTRUCT_HDR
       INTEGER OPT_WARNING_REF_COLLISION, OPT_WARNING_OLDSTYLE
       INTEGER OPT_QUIET_CONSCHECK, OPT_DEBUG_XFERMESGS
       INTEGER OPT_INFO_XFER

       PARAMETER(OPT_IDENTIFY_MODE=0)

       PARAMETER(OPT_WARNING_VARSIZE_OBJ=8)
       PARAMETER(OPT_WARNING_SMALLSIZE=9)
       PARAMETER(OPT_WARNING_PRIOCHANGE=10)
       PARAMETER(OPT_WARNING_DESTRUCT_HDR=11)
       PARAMETER(OPT_WARNING_REF_COLLISION=12)
       PARAMETER(OPT_WARNING_OLDSTYLE=13)

       PARAMETER(OPT_QUIET_CONSCHECK=16)
       PARAMETER(OPT_DEBUG_XFERMESGS=17)
       PARAMETER(OPT_INFO_XFER=18)

       INTEGER OPT_ON, OPT_OFF
       INTEGER IDMODE_LISTS, IDMODE_SETS
       PARAMETER(OPT_OFF=0)
       PARAMETER(OPT_ON=1)
       PARAMETER(IDMODE_LISTS=10)
       PARAMETER(IDMODE_SETS=11)


************************************************************************
*
*  constants for DDD-HANDLERs
*
************************************************************************

       INTEGER HANDLER_LDATACONSTRUCTOR, HANDLER_DESTRUCTOR
       INTEGER HANDLER_DELETE, HANDLER_UPDATE
       INTEGER HANDLER_OBJMKCONS, HANDLER_SETPRIORITY
       INTEGER HANDLER_XFERCOPY, HANDLER_XFERDELETE
       INTEGER HANDLER_XFERGATHER, HANDLER_XFERSCATTER
       INTEGER HANDLER_XFERGATHERX, HANDLER_XFERSCATTERX
       INTEGER HANDLER_ALLOCOBJ, HANDLER_FREEOBJ, HANDLER_END
       PARAMETER(HANDLER_LDATACONSTRUCTOR=0)
       PARAMETER(HANDLER_DESTRUCTOR=1)
       PARAMETER(HANDLER_DELETE=2)
       PARAMETER(HANDLER_UPDATE=3)
       PARAMETER(HANDLER_OBJMKCONS=4)
       PARAMETER(HANDLER_SETPRIORITY=5)
       PARAMETER(HANDLER_XFERCOPY=6)
       PARAMETER(HANDLER_XFERDELETE=7)
       PARAMETER(HANDLER_XFERGATHER=8)
       PARAMETER(HANDLER_XFERSCATTER=9)
       PARAMETER(HANDLER_XFERGATHERX=10)
       PARAMETER(HANDLER_XFERSCATTERX=11)

       PARAMETER(HANDLER_ALLOCOBJ=12)
       PARAMETER(HANDLER_FREEOBJ=13)
       PARAMETER(HANDLER_END=999)
 

************************************************************************
*
*  constants for DDD_IFOneway  
*
************************************************************************


       INTEGER IF_FORWARD, IF_BACKWARD
       PARAMETER(IF_FORWARD=1)
       PARAMETER(IF_BACKWARD=2)


************************************************************************

