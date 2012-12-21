dnl This macro introduces a configure flag --enable-system-heap.
dnl If set UG takes it memory through malloc/free instead of using
dnl the custom heap data structures in low/heaps.c.
dnl
dnl Using the system heap has several advantages:
dnl - there is no artificial upper bound to the amount of memory UG can use
dnl - valgrind may see more bugs
dnl
dnl The system heap may also be faster, but it may also be slower.
dnl No way to know without measuring.
dnl
dnl Since this is an experimental feature it is switched off by default

AC_DEFUN([UG_ENABLE_SYSTEM_HEAP],[
  AC_ARG_ENABLE(system-heap,
    AS_HELP_STRING([--enable-system-heap],[If this is set, the operating system heap is used instead of UG's own heap data structure.]))

  AS_IF([test "x$enable_system_heap" = "xyes"],
    AC_DEFINE(UG_USE_SYSTEM_HEAP, 1, [If this is set, the operating system heap is used instead of UG's own heap data structure.]))
])
