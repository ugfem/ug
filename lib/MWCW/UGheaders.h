// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "compiler.h"

/* low */
#include "bio.h"
#include "debug.h"
#include "defaults.h"
#include "fifo.h"
#include "fileopen.h"
#include "general.h"
#include "heaps.h"
#include "initlow.h"
#include "misc.h"
#include "pfile.h"
#include "scan.h"
//#include "tree.h"
#include "ugenv.h"
#include "ugstruct.h"

/* dev */
#include "initdev.h"
#include "ugdevices.h"

/* dev/mif */
#include "MacGraph.h"
#include "MacMain.h"
#include "MacShell.h"
#include "MacSurface.h"

/* dom */
#include "domain.h"

/* conflicts between domain modules don't allow to iclude those headers
   // dom/std
   #include "std_domain.h"

   // dom/lgm
   #include "ansys2lgm.h"
   #include "lgm_domain.h"
   #include "lgm_load.h"
   #include "lgm_macros.h"
   #include "lgm_transfer.h"
 */

/* gm */
#include "algebra.h"
#include "cw.h"
#include "dlmgr.h"
#include "elements.h"
#include "enrol.h"
#include "er.h"
#include "evm.h"
#include "gm.h"
#include "initgm.h"
#include "mgio.h"
#include "pargm.h"
#include "refine.h"
#include "rm.h"
#include "shapes.h"
#include "ugio.h"
#include "ugm.h"

/* gm/gg2 */
#include "ggaccel.h"
#include "ggm.h"
#include "ggmain.h"

/* np */
#include "initnp.h"
#include "initnumerics.h"
#include "np.h"
#include "num.h"

/* np/udm */
#include "data_io.h"
#include "dio.h"
#include "disctools.h"
#include "formats.h"
#include "npscan.h"
#include "numproc.h"
#include "pcr.h"
#include "udm.h"

/* np/algebra */
#include "amgtools.h"
#include "block.h"
#include "fegeom.h"
#include "ff.h"
#include "ff_gen.h"
#include "fvgeom.h"
#include "npcheck.h"
#include "quadrature.h"
#include "transgrid.h"
#include "ugblas.h"

/* np/amglib */
#include "amg_blas.h"
#include "amg_coarsen.h"
#include "amg_header.h"
#include "amg_iter.h"
#include "amg_low.h"
#include "amg_solve.h"
#include "amg_sp.h"
#include "amg_ug.h"

/* np/field */
#include "field.h"

/* np/procs */
#include "amgtransfer.h"
#include "assemble.h"
#include "basics.h"
#include "bdf.h"
#include "db.h"
#include "error.h"
#include "ew.h"
#include "fas.h"
#include "freebnd.h"
#include "iter.h"
#include "ls.h"
#include "newton.h"
#include "nliter.h"
#include "nls.h"
#include "order.h"
#include "project.h"
#include "transfer.h"
#include "ts.h"
#include "tstep.h"

/* graphics */
#include "graphics.h"

/* graphics/uggraph */
#include "graph.h"
#include "initgraph.h"
#include "plotproc.h"
#include "wop.h"
#include "wpm.h"

/* ui */
#include "avs.h"
#include "cmdint.h"
#include "cmdline.h"
#include "commands.h"
#include "helpmsg.h"
#include "initui.h"
#include "tecplot.h"
#include "uginterface.h"
