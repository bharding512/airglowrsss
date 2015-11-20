/*****************************************************************************/
/*        Copyright (C) Roper Scientific, Inc. 2002 All rights reserved.     */
/*****************************************************************************/
#ifndef _MASTER_H
#define _MASTER_H
static const char *_master_h_="$Header: /PVCAM/SourceLinux/master.h 1     7/18/02 8:24a Dtrent $";


#ifndef linux
    #error OS Not Supported
#endif


#if defined  __cplusplus                                                     
  #define PV_C_PLUS_PLUS
#endif                                  

/**************************** Calling Conventions ****************************/
#define PV_CDECL 
#define PV_DECL PV_CDECL


/**************************** PVCAM Pointer Types ****************************/
#define PV_PTR_DECL  *
#define PV_BUFP_DECL *

/******************************** PVCAM Types ********************************/
enum { PV_FAIL, PV_OK };

typedef int            rs_bool, PV_PTR_DECL  rs_bool_ptr;
typedef char                    PV_PTR_DECL  char_ptr;
typedef signed char    int8,    PV_PTR_DECL  int8_ptr;
typedef unsigned char  uns8,    PV_PTR_DECL  uns8_ptr;
typedef short          int16,   PV_PTR_DECL  int16_ptr;
typedef unsigned short uns16,   PV_PTR_DECL  uns16_ptr;
typedef long           int32,   PV_PTR_DECL  int32_ptr;
typedef unsigned long  uns32,   PV_PTR_DECL  uns32_ptr;
typedef double         flt64,   PV_PTR_DECL  flt64_ptr;
typedef void                    PV_BUFP_DECL void_ptr;
typedef void_ptr                PV_BUFP_DECL void_ptr_ptr;

typedef const rs_bool PV_PTR_DECL rs_bool_const_ptr;
typedef const char    PV_PTR_DECL char_const_ptr;
typedef const int8    PV_PTR_DECL int8_const_ptr;
typedef const uns8    PV_PTR_DECL uns8_const_ptr;
typedef const int16   PV_PTR_DECL int16_const_ptr;
typedef const uns16   PV_PTR_DECL uns16_const_ptr;
typedef const int32   PV_PTR_DECL int32_const_ptr;
typedef const uns32   PV_PTR_DECL uns32_const_ptr;
typedef const flt64   PV_PTR_DECL flt64_const_ptr;

/* deprecated types */
typedef int boolean;
typedef boolean PV_PTR_DECL  boolean_ptr;
typedef const boolean PV_PTR_DECL boolean_const_ptr;

/****************************** PVCAM Constants ******************************/
#ifndef FALSE
  #define FALSE  PV_FAIL      /* FALSE == 0                                  */
#endif

#ifndef TRUE
  #define TRUE   PV_OK        /* TRUE  == 1                                  */
#endif

#define BIG_ENDIAN    FALSE /* TRUE for Motorola byte order, FALSE for Intel */
#define CAM_NAME_LEN     32 /* Max length of a cam name (includes null term) */

/************************ PVCAM-Specific Definitions *************************/
#define MAX_CAM          16 /* Maximum number of cameras on this system.     */

#endif /* _MASTER_H */
