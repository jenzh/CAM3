#ifdef ESMC_RCS_HEADER
"$Id: ESMF_conf.h,v 1.1.8.1 2004/01/02 18:50:57 mvr Exp $"
"Defines the configuration for this machine"
#endif

#if !defined(INCLUDED_CONF_H)
#define INCLUDED_CONF_H
 
#define PARCH_IRIX64 
#define ESMF_ARCH_NAME "IRIX64"
#define ESMC_HAVE_LIMITS_H
#define ESMC_HAVE_PWD_H 
#define ESMC_HAVE_STRING_H 
#define ESMC_HAVE_STROPTS_H 
#define ESMC_HAVE_MALLOC_H 
#define ESMC_HAVE_DRAND48 
#define ESMC_HAVE_GETDOMAINNAME 
#define ESMC_HAVE_UNAME 
#define ESMC_HAVE_UNISTD_H 
#define ESMC_HAVE_STDLIB_H
#define ESMC_HAVE_SYS_TIME_H 
#define ESMC_HAVE_SYS_UTSNAME_H
#define ESMC_USE_SHARED_MEMORY

#undef ESMC_HAVE_OMP_THREADS
#undef ESMC_HAVE_PCL

#undef ESMC_HAVE_MPI

#define ESMC_POINTER_SIZE 8

#define ESMC_SUBSTITUTE_CTRL_CHARS 1

#define ESMC_HAVE_FORTRAN_UNDERSCORE
#define ESMC_SIZEOF_VOIDP 8
#define ESMC_SIZEOF_INT 4
#define ESMC_SIZEOF_DOUBLE 8

#define ESMC_HAVE_IRIXF90

#define ESMC_WORDS_BIGENDIAN 1

#define ESMC_HAVE_MEMMOVE

#define ESMC_HAVE_DOUBLE_ALIGN
#define ESMC_HAVE_DOUBLE_ALIGN_MALLOC

#define ESMC_HAVE_MEMALIGN

#define ESMC_HAVE_FAST_MPI_WTIME

#define ESMC_USE_DBX_DEBUGGER
#define ESMC_HAVE_SYS_RESOURCE_H

#define ESMC_HAVE_RTLD_GLOBAL 1

#define ESMC_CAN_SLEEP_AFTER_ERROR

#define ESMC_HAVE_4ARG_SIGNAL_HANDLER

#define ESMC_USE_KBYTES_FOR_SIZE
#define ESMC_USE_P_FOR_DEBUGGER

#endif
