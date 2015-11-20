//=============================================================================
// File name  : lsbigcam.h
// Description: Header file of the SbigCam class.
// Copyright (C) 2002 - 2003 Jan Soldan, 251 65 Ondrejov-236, Czech Republic
// All rights reserved.
//=============================================================================
#ifndef _LSBIG_CAM_H_
#define _LSBIG_CAM_H_

#include "base.h"
#include "sbigdef.h"
#include "ksbigictl.h"
#include "ksbiglptmain.h"
//=============================================================================
class SbigCam{
protected:
 int                        m_status;
 int                        m_fd;
 CAMERA_TYPE                m_camera_id;
 char                       m_device_name[MAX_DEVICE_NAME_SIZE];
 bool                       m_vdd_set;
 bool                       m_imaging_rip;
 bool                       m_tracking_rip;
 bool                       m_vdd_override;
 bool                       m_st237A;
 StartReadoutParams         m_start_readout_params;
 unsigned char              m_control_out;
 unsigned char              m_imaging_clocks_out;
 unsigned short             m_imaging_bias;
 unsigned short             m_tracking_bias;
 unsigned short             m_st5c_bias;		
 unsigned short             m_st237_bias;	
 short                      m_command_status[CC_LAST_COMMAND];
 EXPOSURE_STATE             m_exposure_state[2];
 bool                       m_link_established;
 CCD_INFO                   m_ccd_info[2];
 MiscellaneousControlParams m_last_control;
 EEPROMContents             m_eeprom;
 bool                       m_shutter_open;
 bool                       m_wait_for_trigger;
 CCD                        m_triggered_ccd;
 StartExposureParams        m_trigger_exp_params;
 unsigned long              m_last_tracking_time;
 bool                       m_te_frozen;
 bool                       m_te_frozen_off;
 unsigned short             m_last_te_setpoint;
 unsigned short             m_last_te_power;
 bool                       m_te_auto_freeze;
 bool                       m_st253;
 bool                       m_kai311;
 bool                       m_close_shutter_on_end_exposure;
 CCD_REQUEST                m_max_ccd;
 short                      m_has_tracking_ccd;
 short                      m_base_is_st7;
 short                      m_base_is_st5c;
 short                      m_tracker_is_st237;
 unsigned char              m_shutter_edge;
 unsigned char              m_diddle_line_counter;
 unsigned short             m_threshold;
 unsigned short             m_saturation;
 unsigned short             m_offset;
 unsigned short             m_hot_threshold;
 long                       m_hot_count;
 TXSerialBytesParams        m_tx_serial_bytes_params;
 UDRV_OPTIONS               m_udrv_opt;
 CAMERA_TYPE                m_romMSNtoID[6];
 unsigned short             m_temp_video [MAX_DIG_WIDTH];
 unsigned short             m_temp_video2[MAX_DIG_WIDTH];
 unsigned short             m_last_line1 [MAX_HP_WIDTH];
 unsigned short             m_last_line2 [MAX_HP_WIDTH];
 unsigned char              m_microblock_in_buf[BUFSIZE];
 unsigned char              m_microblock_out_buf[512+4];
 unsigned short             m_pixelFifo[USB_FIFO_SIZE];
 FifoInfo                   m_fifoInfo;
 USBIGA                     m_usbIGA;
 unsigned long              m_driverControlParams[DCP_LAST];

public:
                 SbigCam();
 virtual        ~SbigCam();
 unsigned long   GetJiffies();
 unsigned long   GetHz();
 virtual int     GetDriverInfo(GetDriverInfoParams *, void *)             = 0;
 virtual int     TestCommand()                                            = 0;
 int             CalcImageParams(GrabImageParams *);

 /*
 int             Dimm(DimmParams *, unsigned short **);
 int             Dimm(DimmParams *, IMG_DATA_TYPE  **);
 int             CheckDriftScanMode(DriftScanIo *);	      
 int             GrabDriftScanImage(DriftScanIo *, unsigned short **);
 int             GrabDriftScanImage(DriftScanIo *, IMG_DATA_TYPE  **);
 int             SetFilter(int filter);
 int             GrabImage(GrabImageParams *, unsigned short **);
 int             GrabImage(GrabImageParams *, IMG_DATA_TYPE  **);
 */

 int             QueryCommandStatus(QueryCommandStatusParams  *,
	                            QueryCommandStatusResults *);
 int             EstablishLink();
 int             StartExposure(StartExposureParams *);
 int             EndExposure  (EndExposureParams   *);
 int             DumpLines    (DumpLinesParams     *);
 int             ReadoutLine  (ReadoutLineParams   *,
                               unsigned short      *results,
                               bool                 subtract);
 int             SetTemperatureRegulation(SetTemperatureRegulationParams *);
 int             QueryTemperatureStatus  (QueryTemperatureStatusResults  *);
 bool            DetectSt237A();
 int             OpenDevice(char *name);
 int             CloseDevice();
 int             ActivateRelay(ActivateRelayParams *);
 int             PulseOut(PulseOutParams *);
 int             TXSerialBytes(TXSerialBytesParams *, TXSerialBytesResults *);
 int             GetSerialStatus(GetSerialStatusResults *);
 int             AOTipTilt(AOTipTiltParams   *);
 int             AOSetFocus(AOSetFocusParams *);
 void            AODelay(AODelayParams       *);
 void            GetLinkStatus(LinuxGetLinkStatusResults *);
 unsigned short  CalcSetpoint(double temperature);
 double          CalcCcdTemperature(QueryTemperatureStatusResults *);
 double          CalcAmbTemperature(QueryTemperatureStatusResults *);
 int             CalcPower(QueryTemperatureStatusResults *);
 char           *GetStatusString();
 void            InitCCDInfo();
 void            InitTrackingCCDInfo();
 int             Bcd2Int(unsigned short bcd);
 void            TimerDelay(unsigned long us);
 int             MiscellaneousControl(MiscellaneousControlParams *);
 int             GetCCDInfo(GetCCDInfoParams *, void *);
 inline char    *GetDeviceName()      {return(m_device_name);      }
 inline int      GetCameraType()      {return(m_camera_id);        }
 inline bool     GetLinkEstablished() {return(m_link_established); }
 inline int      GetStatus()          {return(m_status);           }
 inline void     ClearStatus()        {m_status = CE_NO_ERROR;     }
 inline bool     IsDeviceOpen(){return((m_fd) ? true : false);     } 
   
protected:
 int             ReadOffset(ReadOffsetParams *, ReadOffsetResults *);
 
 int             GetEEPROM(EEPROMContents *);

 int             PutEEPROM(EEPROMContents *);

 unsigned short  CalcEEPROMChecksum(EEPROMContents *);

 int             MicroCommand(MICRO_COMMAND command,
                              void          *txDataPtr,
                              void          *rxDataPtr);

 int             ValidateMicroResponse(MICRO_COMMAND command,
         	                       void          *dataPtr,
                                       void          *txDataPtr,
                                       unsigned long rxLen);

 int             ValGetMicroAck();

 int             ValGetMicroBlock(MICRO_COMMAND command,
                                  void          *rxDataPtr,
                                  void          *txDataPtr,
                                  unsigned long rxLen);

 void            SetVdd(bool raiseIt);

 void            BuildMicroCommand(MICRO_COMMAND command,
                                   unsigned char *dataPtr,
                                   unsigned long *pTxLen,
                                   unsigned long *pRxLen);

 int             CheckShutter();

 virtual void    CameraOut(unsigned char reg,  unsigned char val)         = 0;
 virtual void    VirtBuildMicroCommand1(unsigned long  &exposureTime)     = 0;
 virtual void    VirtBuildMicroCommand2(unsigned short &len)              = 0;
 virtual void    VirtBuildMicroCommand3(unsigned char *p, EEPROMParams *) = 0;
 virtual int     FreezeTEControl(bool freezeIt)        = 0;
 virtual void    VirtEndExposure()                     = 0;
 virtual void    InitPort()                            = 0;
 virtual void    SetTrackerIsSt237()                   = 0;
 virtual void    VirtReadOffset()                      = 0;
 virtual void    SetRxLen(unsigned long &rxLen)        = 0;
 virtual int     VirtValGetMicroBlock(
                 MICRO_COMMAND command,
                 unsigned char subcommand,
		 unsigned long cmp_len,
		 unsigned char *pRxData,
                 unsigned long &rx_len)                = 0;

 virtual int     SendMicroBlock(unsigned long length)  = 0;
     
 virtual bool    VirtSetVdd(bool raiseIt)              = 0;
 
 virtual int     GetMicroBlock(unsigned char *pBuf, 
                               unsigned long length)   = 0;

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
                                 short           clearWidth) = 0;

 virtual int     CCDDumpLines(CCD_REQUEST ccd,
                              short       width,
	           	      short       len,
                              short       vertBin)           = 0; 

 virtual int     ClearITArray(CCD_REQUEST ccd,
                              short       height,
                              short       width,
                              short       times,
                              short       left)              = 0;

 virtual int     OffsetITArray(CCD_REQUEST     ccd,
                               short           height,
                               short           width,
                               unsigned short *offset,
                               short           left)         = 0;

 virtual int     OffsetST5CArray(short           width, 
		   	         unsigned short *offset,
	                         unsigned short  mask)       = 0;
							   
 virtual int     CCDMeasureBias(CCD_REQUEST ccd,
	      		        short       clearWidth,
                                short       left)            = 0;

 void            make_mode(GetCCDInfoResults0 *r, 
                           short ccd, 
			   short position, 
                           short binning);

 void            RelayClick();
 int             GetBufferSize();
 int             SetBufferSize(unsigned short new_size);
 int             GetFileDescriptor()      {return(m_fd);}
 bool            GetVdd()                 {return(m_vdd_set);}
 bool            GetVddOverride()         {return(m_vdd_override);}
 bool            GetSt237A()              {return(m_st237A);}
 CCD_REQUEST     GetMaxCcd()              {return(m_max_ccd);}
 int             GetImagingBias()         {return(m_imaging_bias);}
 int             GetTrackingBias()        {return(m_tracking_bias);}
 int             GetSt5cBias()            {return(m_st5c_bias);}		
 int             GetSt237Bias()           {return(m_st237_bias);}		
 int             GetEepromChecksum()      {return(m_eeprom.checksum);}
 int             GetEepromVersion()       {return(m_eeprom.version);}
 int             GetEepromModel()         {return(m_eeprom.model);}
 int             GetEepromAbgType()       {return(m_eeprom.abgType);}
 int             GetEepromBadColumns()    {return(m_eeprom.badColumns);}
 int             GetEepromTrackingOffset(){return(m_eeprom.trackingOffset);}
 int             GetEepromTrackingGain()  {return(m_eeprom.trackingGain);}
 int             GetEepromImagingOffset() {return(m_eeprom.imagingOffset);}
 int             GetEepromImagingGain()   {return(m_eeprom.imagingGain);}
 char           *GetEepromSerialNumber()  {return((char*)m_eeprom.serialNumber);}
 void            InitVars();
 double          CalcTemperature(short thermistorType, short ccdSetpoint);
 void            make_n_modes(GetCCDInfoResults0 *r, short ccd, 
                              short start, short n);
 void            make_by_n_modes(GetCCDInfoResults0 *r, short 
                                 start, short n);
 void            clear  (char *dest, short len);
 void            scpy   (char *dest, char *src);
 void            cpy    (char *dest, char *src, short len);
 void            swapcpy(char *dest, char *src, short len);
 unsigned short  bcd_nx(unsigned short bcd, short n);
 void            OffVertBinPixels(unsigned short *dest,
                                  unsigned short *src,
		    	          short           len,
                                  unsigned short  threshold,
                                  unsigned short  saturation);

 void            OffHorzBinPixels(unsigned short *dest,
                                  unsigned short *src,
	    		          short           len,
                                  short           offHorzBin,
                                  unsigned short  threshold,
			          unsigned short  saturation);

 void            HFlipPixels(unsigned short *dest, short len);

 void            SubtractPixels(unsigned short *dest,
                                unsigned short *src,
                                short           len,
	   		        unsigned short  threshold,
                                unsigned short  saturation);

 void            OffsetPixels(unsigned short *dest,
                              short           len,
                              unsigned short  offset,
	    		      unsigned short  threshold,
                              unsigned short  saturation);

 short           HotPatchPixels(unsigned short *temp_video,
                                unsigned short *last_line1,
	   		        unsigned short *last_line2,
                                short           len,
	    		        unsigned short  hot_threshold);

 void            SetSatThr(CCD_REQUEST ccd, short offFactor);

 unsigned char   CalculateSerialChecksum(void);

};
//=============================================================================
#endif //_LSBIG_CAM_H_



