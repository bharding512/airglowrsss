//========================================================================
// File name  : ksbiglpt.h
// Description: Function prototypes for kernel LPT code.
// Author     : Jan Soldan - Linux
// Author     : Matt Longmire (SBIG) - Windows
// Copyright (C) 1999 - 2003 Matt Longmire, Jan Soldan
// All rights reserved.
//========================================================================
#ifndef _KSBIG_LPT_H_
#define _KSBIG_LPT_H_

#include "sbigdef.h"

// Kernel functions prototypes
void      KLptInitPort  (
                         struct private_data  *pd
                        );

int       KLptCameraOutWrapper
                        (
                         struct private_data  *pd, 
                         LinuxCameraOutParams *lcop
			);
void      KLptCameraOut           
                        (
			 struct private_data  *pd, 
                         unsigned char         reg, 
			 unsigned char         val   
			);

int       KLptSendMicroBlock      
                        (
			 struct private_data  *pd, 
			 LinuxMicroblock      *lmb
			);

int       KLptGetMicroBlock       
                        (
			 struct private_data  *pd, 
			 LinuxMicroblock      *lmb
			);

int       KLptSetVdd              
                        (
			 struct private_data  *pd, 
			 IocSetVdd            *svdd
			);

int       KLptClearImagingArray   
                        (
			 struct private_data  *pd, 
			 IOC_CLEAR_CCD_PARAMS *cp
			);

int       KLptClearTrackingArray  
                        (
			 struct private_data  *pd, 
			 IOC_CLEAR_CCD_PARAMS *cp
			);

int       KLptGetPixels           
                        (
			 struct private_data  *pd, 
			 LinuxGetPixelsParams *lgpp
			);

int       KLptGetArea             
                        (
			 struct private_data  *pd, 
			 LinuxGetAreaParams   *lgpp
			);

int       KLptGetDriverInfo       
                        (
			 GetDriverInfoResults0 *gdir0
			);

int       KLptRVClockImagingCCD   
                        (
			 struct private_data   *pd, 
			 IOC_VCLOCK_CCD_PARAMS *cp
			);

int       KLptRVClockTrackingCCD  
                        (
			 struct private_data   *pd, 
			 IOC_VCLOCK_CCD_PARAMS *cp
			);

int       KLptRVClockST5CCCD      
                        (
			 struct private_data   *pd, 
			 IOC_VCLOCK_CCD_PARAMS *cp
			);

unsigned char     KLptCameraIn            
                        (
			 struct private_data   *pd, 
			 unsigned char         reg              
			);

void      KLptIoDelay             
                        (
			 struct private_data   *pd, 
			 short                 i
			);

unsigned char     KLptMicroIn             
                        (
			 struct private_data   *pd, 
			 unsigned char         ackIt
			);

void      KLptMicroOut            
                        (
			 struct private_data   *pd, 
			 unsigned char         val
			);

int       KLptHClear              
                        (
			 struct private_data   *pd, 
			 short                 times            
			);

int       KLptVClockImagingCCD    
                        (
			 struct private_data   *pd, 
			 CAMERA_TYPE           cameraID,	
			 unsigned char         baseClks,	  
			 short                 hClears          
			);

void      KDisable                
                        (
			 struct private_data *pd
			);

void      KEnable        
                        (
			 struct private_data *pd
			);

void      KLptDisableLptInterrupts
                        (
			 struct private_data *pd
			);

void      KLptReadyToRx           
                        (
			 struct private_data *pd
			);

void      KLptForceMicroIdle      
                        (
			 struct private_data *pd
			);

unsigned char     KLptMicroStat           
                        (
			 struct private_data *pd
			);

int       KLptWaitForPLD          
                        (
			 struct private_data *pd
			);

int       KLptWaitForAD           
                        (
			 struct private_data *pd
			);

int       KLptGetJiffies          
                        (
			 unsigned long *arg
			);

int       KLptGetHz          
                        (
			 unsigned long *arg
			);

int       KLptDumpImagingLines    
                        (
			 struct private_data   *pd, 
			 IOC_DUMP_LINES_PARAMS *dlp
			);

int       KLptDumpTrackingLines   
                        (
			 struct private_data   *pd, 
			 IOC_DUMP_LINES_PARAMS *dlp
			);

int       KLptDumpST5CLines       
                        (
			 struct private_data   *pd, 
			 IOC_DUMP_LINES_PARAMS *dlp
			);

int       KLptClockAD             
                        (
			 struct private_data *pd, 
			 short               *len
			);

int       KLptSetBufferSize       
                        (
			 struct private_data *pd, 
			 spinlock_t          *lock, 
                         unsigned short      *new_size
			);

int       KLptGetBufferSize       
                        (
			 struct private_data *pd
			);

int       KLptTestCommand
                        (
			 void
			);

int       KDevOpen      
                        (
			 struct inode *inode,
                         struct file  *filp,
	                 int           lpt_base,
	                 int 	       lpt_span,
	                 int           buffer_size
			);

int       KDevRelease             
                        (
			 struct inode *inode,
                         struct file  *filp
			);

int       KDevIoctl     
                        (
                         struct inode  *inode,
                         struct file   *filp,
                         unsigned int   cmd,
                         unsigned long  arg,
			 spinlock_t    *spin_lock
			);
//========================================================================
#endif  // _KSBIG_LPT_H_

