dnl This macro introduces a configure flag --disable-system-heap.
dnl If set, UG uses a custom heap data structure in low/heaps.c instead
dnl of taking its memory through malloc/free.
dnl
dnl Using the system heap has several advantages:
dnl - there is no artificial upper bound to the amount of memory UG can use
dnl - valgrind may see more bugs
dnl
dnl The system heap may also be faster, but it may also be slower.
dnl No way to know without measuring.

AC_DEFUN([UG_ENABLE_SYSTEM_HEAP],[
  AC_ARG_ENABLE(system-heap,
    AS_HELP_STRING([--disable-system-heap],[If this is set, UG's own heap data structure is used instead of the operating system heap.]), [], [enable_system_heap=yes])

  AS_IF([test "x$enable_system_heap" = "xyes"],
    AC_DEFINE(UG_USE_SYSTEM_HEAP, 1, [If this is set, the operating system heap is used instead of UG's own heap data structure.]))
])
