/* config/LSMLIB_config.h.  Generated from LSMLIB_config.h.in by configure.  */
/*
 * File:        LSMLIB_config.h.in
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: LSMLIB configuration header file
 */


#ifndef included_LSMLIB_config_h
#define included_LSMLIB_config_h

/* Debugging macro */
#ifndef LSMLIB_DEBUG_NO_INLINE
/* #undef LSMLIB_DEBUG_NO_INLINE */
#endif

/* Macro defined if double precision library is being built. */
#ifndef LSMLIB_DOUBLE_PRECISION
#define LSMLIB_DOUBLE_PRECISION 1
#endif

/* Floating-point precision for LSMLIB_REAL */
#ifndef LSMLIB_REAL
#define LSMLIB_REAL double
#endif

/* Zero tolerance */
#ifndef LSMLIB_ZERO_TOL
#define LSMLIB_ZERO_TOL 1.e-11
#endif

/* Maximum value for LSMLIB_REAL */
#ifndef LSMLIB_REAL_MAX
#define LSMLIB_REAL_MAX DBL_MAX
#endif

/* Minimum value for LSMLIB_REAL */
#ifndef LSMLIB_REAL_MIN
#define LSMLIB_REAL_MIN DBL_MIN
#endif

/* Machine epsilon value for LSMLIB_REAL */
#ifndef LSMLIB_REAL_EPSILON
#define LSMLIB_REAL_EPSILON DBL_EPSILON
#endif

#endif

