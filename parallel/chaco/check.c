// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* debug_graph */
{
  extern int DEBUG_GRAPH;

  /* debug_graph */
  if (DEBUG_GRAPH>0)
  {
    int flag;
    char buf[100];
    flag = check_graph(graph,nvtxs,nedges);
    if (flag == FALSE)
    {
      sprintf(buf,"FATAL: check_graph returned with FALSE\n");
      UserWrite(buf);
      sprintf(buf,"nvtxs=%d,nedges=%d,ndims=%d\n",nvtxs,nedges,ndims);
      UserWrite(buf);
    }
    else
    {
      sprintf(buf,"OK: check_graph = TRUE in eigensolve() before pert,nvtxs=%d,nedges=%d\n",nvtxs,nedges);
      UserWrite(buf);
    }

  }
}
