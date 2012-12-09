AC_DEFUN([DDD_PARAMETERS],[
  AC_ARG_WITH([ddd_maxprocbits],
    AS_HELP_STRING([--with-ddd-maxprocbits],
	  [Set number of bits of an unsigned int used to store the process number, 
       the remaining bits are used to store the local entity id]),[],[with_ddd_maxprocbits=no])
  AS_IF([test "x$with_ddd_maxprocbits" = "xno"],[
    with_ddd_maxprocbits=24;
  ])
  AC_DEFINE_UNQUOTED(DDD_MAX_PROCBITS_IN_GID, $with_ddd_maxprocbits,
                      [see parallel/ddd/dddi.h])
])
