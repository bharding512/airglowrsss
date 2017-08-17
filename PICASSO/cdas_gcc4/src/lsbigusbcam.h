//========================================================================
// File name  : lsbigusbcam.h
// Description: Header file of the SbigUsbCam class.
// Copyright (C) 2002 - 2003 Jan Soldan, 251 65 Ondrejov-236, Czech Republic
// All rights reserved.
//========================================================================
#ifndef _LSBIG_USB_CAM_H_
#define _LSBIG_USB_CAM_H_

#include "lsbigcam.h"
//========================================================================
class SbigUsbCam : public SbigCam{
public:
               SbigUsbCam();
               SbigUsbCam(char *name);
              ~SbigUsbCam();
 int           GetDriverInfo(GetDriverInfoParams *, void *);
 int           TestCommand();

 protected:
 int           UsbDumpLines(CCD_REQUEST ccd, 
                            short       width, 
			    short       len, 
			    short       vertBin);

 int           UsbGetPixels(CCD_REQUEST     ccd, 
                            unsigned short *dest,
		            short           left, 
			    short           len, 
			    short           right, 
			    short           horzBin);

 int           MicroInitPixelReadout(CCD_REQUEST ccd, 
                                     short       left, 
				     short       noPixels,
                                     short       right, 
				     short       windowHeight, 
				     short       horzBin, 
				     short       vertBin);

 int           MicroGetPixels(unsigned short *dest);

 virtual void  VirtBuildMicroCommand1(unsigned long &exposureTime);

 virtual void  VirtBuildMicroCommand2(unsigned short &len){ len = 3;}

 virtual void  VirtBuildMicroCommand3(unsigned char *p,
                                      EEPROMParams  *eepp)
	       {
	         *p++ = eepp->deviceAddress; 
               }

 virtual int   FreezeTEControl(bool freezeIt);

 virtual bool  VirtSetVdd(bool raiseIt);

 virtual void  SetRxLen(unsigned long &rxLen){ rxLen = 2;}

 virtual int   CCDDigitizeLine(CCD_REQUEST     ccd,
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

 virtual int   CCDDumpLines(CCD_REQUEST ccd,
                            short       width,
	        	    short       len,
                            short       vertBin);
			    
 virtual int   ClearITArray(CCD_REQUEST ccd,
                            short       height,
                            short       width,
                            short       times,
                            short       left);

 virtual int   OffsetITArray(CCD_REQUEST     ccd,
                             short           height,
                             short           width,
                             unsigned short *offset,
                             short           left);

 virtual int   OffsetST5CArray(short           width, 
	         	       unsigned short *offset,
	                       unsigned short  mask);

 virtual int   CCDMeasureBias(CCD_REQUEST ccd,
	        	      short       clearWidth,
                              short       left);

 virtual int   SendMicroBlock(unsigned long length);

 virtual int   GetMicroBlock(unsigned char *pBuf, 
                             unsigned long length);
 
 virtual void  VirtEndExposure();

 virtual int   VirtValGetMicroBlock(MICRO_COMMAND command,
                                    unsigned char subcommand,
				    unsigned long cmp_len,
				    unsigned char *pRxData,
                                    unsigned long &rx_len);

 virtual void  InitPort(){}
 virtual void  SetTrackerIsSt237();
 virtual void  VirtReadOffset();

 virtual void  CameraOut(unsigned char reg,  unsigned char val){}
};
//========================================================================
#endif //_LSBIG_USB_CAM_H_
