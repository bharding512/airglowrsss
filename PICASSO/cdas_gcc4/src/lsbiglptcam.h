//========================================================================
// File name  : lsbiglptcam.h
// Description: Header file of the SbigLptCam class.
// Copyright (C) 2002 - 2003 Jan Soldan, 251 65 Ondrejov-236, Czech Republic
// All rights reserved.
//========================================================================
#ifndef _LSBIG_LPT_CAM_H_
#define _LSBIG_LPT_CAM_H_

#include "lsbigcam.h"
//========================================================================
class SbigLptCam : public SbigCam{
public:
                 SbigLptCam();
                 SbigLptCam(char *device_name);
                ~SbigLptCam();
 virtual int     GetDriverInfo(GetDriverInfoParams  *params,
	                       void                 *results);


 int             ReallocateLptPorts(LptPortParams   *params);
 int             StartReadout (StartReadoutParams   *params);
 /*
 int             ReadoutArea(LinuxReadoutAreaParams *params,
                             unsigned long          *results,
                             bool                   subtract);
 */
 int             EndReadout(EndReadoutParams *params);
 int             TestCommand();

 protected:

 virtual void    VirtBuildMicroCommand1(unsigned long  &exposureTime){}

 virtual void    VirtBuildMicroCommand2(unsigned short &len){len = 2;}

 virtual void    VirtBuildMicroCommand3(unsigned char *p,
                                        EEPROMParams  *params){}

 virtual int     FreezeTEControl(bool freezeIt);

 virtual bool    VirtSetVdd(bool raiseIt);

 virtual int     CCDDigitizeLine(CCD_REQUEST     ccd,
		   	         short           left,
                                 short           len,
                                 short           right,
                                 short           windowHeight,
			         short           onHorzBin,
                                 short           offHorzBin,
			         short           onVertBin,
                                 short           offVertBin,
                                 unsigned short *dest,
			         bool            subtract,
                                 short           clearWidth);

 virtual int     CCDDumpLines(CCD_REQUEST ccd,
                              short       width,
	        	      short       len,
                              short       vertBin);

 virtual int     ClearITArray(CCD_REQUEST ccd,
                              short       height,
                              short       width,
                              short       times,
                              short       left);

 virtual int     OffsetITArray(CCD_REQUEST     ccd,
                               short           height,
                               short           width,
                               unsigned short *offset,
                               short           left);

 virtual int     OffsetST5CArray(short           width, 
	         		 unsigned short *offset,
	                         unsigned short  mask);

 virtual int     CCDMeasureBias(CCD_REQUEST ccd,
	        		short       clearWidth,
                                short       left);

 virtual int     SendMicroBlock(unsigned long length);

 virtual int     GetMicroBlock(unsigned char *pBuf, 
                               unsigned long length);
 
 virtual void    VirtEndExposure();

 virtual void    SetRxLen(unsigned long &rxLen){rxLen = 1;}

 virtual int     VirtValGetMicroBlock(MICRO_COMMAND command,
                                      unsigned char subcommand,
				      unsigned long cmp_len,
				      unsigned char *pRxData,
                                      unsigned long &rx_len);

 virtual void    InitPort();
 virtual void    SetTrackerIsSt237();
 virtual void    VirtReadOffset(){}

 virtual void    CameraOut(unsigned char reg,  unsigned char val);

 int             LptClearImagingArray(short height,
                                      short times);

 int             LptClearTrackingArray(short height,
                                       short times);

 int             LptGetPixels(CCD_REQUEST  ccd,
                              unsigned short *dest,
	       		      short           left,
                              short           len,
                              short           right,
                              short           horzBin,
			      short           vertBin,
                              short           clearWidth);

 int             CCDLMeasureBias(CCD_REQUEST  ccd,
		 	         unsigned short  width,
                                 unsigned short *pTemp,
			         short           clearWidth,
                                 short           left);


 /*
 int             CCDDigitizeArea(CCD_REQUEST ccd,
		  	         short           left,
                                 short           len,
                                 short           right,
                                 short           windowHeight,
			         short           onHorzBin,
                                 short           offHorzBin,
			         short           onVertBin,
                                 short           offVertBin,
                                 unsigned short *dest,
			         bool            subtract,
                                 short           clearWidth);
 */				           

 int             LptDumpImagingLines(short width,
                                     short len,
                                     short vertBin,
				     short vToHRatio);

 int             LptDumpTrackingLines(short width,
                                      short len,
                                      short vertBin);

 int             LptDumpST5CLines(short width,
                                  short len,
                                  short vertBin);
};
//========================================================================
#endif //_LSBIG_LPT_CAM_H_
