# $Id$

SOURCES = dgssv.c dgssvx.c dsp_blas2.c dsp_blas3.c dgscon.c dlacon.c dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c dgsrfs.c dgstrf.c dgstrs.c dcopy_to_ucol.c dsnode_dfs.c dsnode_bmod.c dpanel_dfs.c dpanel_bmod.c dreadhb.c dcolumn_dfs.c dcolumn_bmod.c dpivotL.c dpruneL.c dmemory.c dutil.c superlu_timer.c lsame.c util.c memory.c get_perm_c.c mmd.c sp_coletree.c sp_preorder.c sp_ienv.c relax_snode.c xerbla.c colamd.c dlamch.c

if BUILD2D
  LIB_2D = librealslu2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = librealslu3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

librealslu2_la_SOURCES = $(SOURCES)
librealslu2_la_CPPFLAGS = $(UG_2DFLAGS) -DAdd_

librealslu3_la_SOURCES = $(SOURCES)
librealslu3_la_CPPFLAGS = $(UG_3DFLAGS) -DAdd_
