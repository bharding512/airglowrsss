#if !defined(__atmcdLXd_h)
#define __atmcdLXd_h

#ifdef __cplusplus
extern "C" {
#endif

#define at_u16 unsigned short
#ifdef _LP64
#define at_32 int
#define at_u32 unsigned int
#else
#define at_32 long
#define at_u32 unsigned long
#endif


typedef struct ANDORCAPS
{
  at_u32 ulSize;
  at_u32 ulAcqModes;
  at_u32 ulReadModes;
  at_u32 ulTriggerModes;
  at_u32 ulCameraType;
  at_u32 ulPixelMode;
  at_u32 ulSetFunctions;
  at_u32 ulGetFunctions;
  at_u32 ulFeatures;
  at_u32 ulPCICard;
  at_u32 ulEMGainCapability;
} AndorCapabilities;

typedef struct COLORDEMOSAICINFO
{
  int iX; // Number of X pixels. Must be >2.
  int iY; // Number of Y pixels. Must be >2.
  int iAlgorithm;  // Algorithm to demosaic image.
  int iXPhase;   // First pixel in data (Cyan or Yellow/Magenta or Green). For Luca/SurCam = 1.
  int iYPhase;   // First pixel in data (Cyan or Yellow/Magenta or Green). For Luca/SurCam = 0.
  int iBackground;  // Background to remove from raw data when demosaicing. Internal Fixed Baseline Clamp for Luca/SurCam is 512
} ColorDemosaicInfo;

typedef struct WHITEBALANCEINFO
{
  int iSize;  // Structure size
  int iX; // Number of X pixels. Must be >2.
  int iY; // Number of Y pixels. Must be >2.
  int iAlgorithm;  // Algorithm to calculate white balance.
  int iROI_left;   // Region Of Interest from which the White Balance is calculated
  int iROI_right;  // Region Of Interest from which the White Balance is calculated
  int iROI_top;    // Region Of Interest from which the White Balance is calculated
  int iROI_bottom; // Region Of Interest from which the White Balance is calculated
  int iOperation; // 0: do NOT calculate white balance, apply given relative factors to colour fields
                  // 1: calculate white balance, do NOT apply calculated relative factors to colour fields
                  // 2: calculate white balance, apply calculated relative factors to colour fields
} WhiteBalanceInfo;

unsigned int AbortAcquisition(void);
unsigned int CancelWait();
unsigned int CoolerOFF(void);
unsigned int CoolerON(void);
unsigned int DemosaicImage(unsigned short* _input, unsigned short* _red, unsigned short* _green, unsigned short* _blue, ColorDemosaicInfo* _info);
unsigned int FreeInternalMemory();
unsigned int GetAcquiredData16(unsigned short* array, at_u32 size);
unsigned int GetAcquiredData(at_32 * array, at_u32 size);
unsigned int GetAcquiredFloatData(float* array, at_u32 size);
unsigned int GetAcquisitionProgress(at_32* acc, at_32* series);
unsigned int GetAcquisitionTimings(float* exposure, float* accumulate, float* kinetic);
unsigned int GetAdjustedRingExposureTimes(int _inumTimes, float *_fptimes);
unsigned int GetAllDMAData(at_32* array, at_u32 size);
unsigned int GetAmpDesc(int index, char* name, int len);
unsigned int GetAmpMaxSpeed(int index, float* speed);
unsigned int GetAvailableCameras(at_32* totalCameras);
unsigned int GetBackground(at_32* array, at_u32 size);
unsigned int GetBitDepth(int channel, int* depth);
unsigned int GetCameraEventStatus(at_u32 *cam_status);
unsigned int GetCameraHandle(at_32 cameraIndex, at_32* cameraHandle);
unsigned int GetCameraInformation(int index, at_32 *information);
unsigned int GetCameraSerialNumber(int* number);
unsigned int GetCapabilities(AndorCapabilities* caps);
unsigned int GetControllerCardModel(char* controllerCardModel);
unsigned int GetCurrentCamera(at_32* cameraHandle);
unsigned int GetCYMGShift(int * _iXshift, int * _iYshift);
unsigned int GetDDGIOCFrequency(double* frequency);
unsigned int GetDDGIOCNumber(unsigned long* number_pulses);
unsigned int GetDDGIOCPulses(int* pulses);
unsigned int GetDDGPulse(double width, double resolution, double* Delay, double *Width);
unsigned int GetDetector(int* xpixels, int* ypixels);
unsigned int GetDICameraInfo(void *info);
unsigned int GetEMCCDGain(int* gain);
unsigned int GetEMGainRange(int* low, int* high);
unsigned int GetFastestRecommendedVSSpeed(int *index, float* speed);
unsigned int GetFIFOUsage(int* FIFOusage);
unsigned int GetFilterMode(int* mode);
unsigned int GetFKExposureTime(float* time);
unsigned int GetFKVShiftSpeedF(int index, float* speed);
unsigned int GetFKVShiftSpeed(int index, int* speed);
unsigned int GetHardwareVersion(unsigned int* PCB, unsigned int* Decode, unsigned int* SerPar, unsigned int* Clocks, unsigned int* dummy1, unsigned int* dummy2);
unsigned int GetHeadModel(char* name);
unsigned int GetHorizontalSpeed(int index, int* speed);
unsigned int GetHSSpeed(int channel, int type, int index, float* speed);
unsigned int GetHVflag(int *bFlag);
unsigned int GetID(int devNum, int* id);
unsigned int GetImages16 (at_32 first, at_32 last, unsigned short* array, at_u32 size, at_32* validfirst, at_32* validlast);
unsigned int GetImages (at_32 first, at_32 last, at_32* array, at_u32 size, at_32* validfirst, at_32* validlast);
unsigned int GetImagesPerDMA(at_u32* images);
unsigned int GetIRQ(int* IRQ);
unsigned int GetMaximumBinning(int ReadMode, int HorzVert, int* MaxBinning);
unsigned int GetMaximumExposure(float* MaxExp);
unsigned int GetMCPGain(int iNum, int *iGain, float *iPhotoepc);
unsigned int GetMCPVoltage(int *iVoltage);
unsigned int GetMinimumImageLength(int* MinImageLength);
unsigned int GetMostRecentColorImage16 (at_u32 size, int algorithm, unsigned short* red, unsigned short* green, unsigned short* blue);
unsigned int GetMostRecentImage16 (unsigned short* array, at_u32 size);
unsigned int GetMostRecentImage (at_32* array, at_u32 size);
//unsigned int GetMSTimingsEnabled(void);
//unsigned int GetMSTimingsData(void *TimeOfStart,float *Differences, int noOfimages);
unsigned int GetNewData16(unsigned short* array, at_u32 size);
unsigned int GetNewData8(unsigned char* array, at_u32 size);
unsigned int GetNewData(at_32* array, at_u32 size);
unsigned int GetNewFloatData(float* array, at_u32 size);
unsigned int GetNumberADChannels(int* channels);
unsigned int GetNumberAmp(int* amp);
unsigned int GetNumberAvailableImages (at_32* first, at_32* last);
unsigned int GetNumberDevices(int* numDevs);
unsigned int GetNumberFKVShiftSpeeds(int* number);
unsigned int GetNumberHorizontalSpeeds(int* number);
unsigned int GetNumberHSSpeeds(int channel, int type, int* speeds);
unsigned int GetNumberNewImages (at_32* first, at_32* last);
unsigned int GetNumberPreAmpGains(int* noGains);
unsigned int GetNumberRingExposureTimes(int *_ipnumTimes);
unsigned int GetNumberVerticalSpeeds(int* number);
unsigned int GetNumberVSAmplitudes(int* number);
unsigned int GetNumberVSSpeeds(int* speeds);
unsigned int GetOldestImage16 (unsigned short* array, at_u32 size);
unsigned int GetOldestImage (at_32* array, at_u32 size);
unsigned int GetPhysicalDMAAddress(at_u32* Address1, at_u32* Address2);
unsigned int GetPixelSize(float* xSize, float* ySize);
unsigned int GetPreAmpGain(int index, float* gain);
unsigned int GetRegisterDump(int* mode);
unsigned int GetRingExposureRange(float *_fpMin, float *_fpMax);
unsigned int GetSizeOfCircularBuffer (at_32* index);
unsigned int GetSlotBusDeviceFunction(at_u32 *dwSlot, at_u32 *dwBus, at_u32 *dwDevice, at_u32 *dwFunction);
unsigned int GetSoftwareVersion(unsigned int* eprom, unsigned int* coffile, unsigned int* vxdrev, unsigned int* vxdver, unsigned int* dllrev, unsigned int* dllver);
unsigned int GetSpoolProgress(at_32* index);
unsigned int GetStatus(int* status);
unsigned int GetTemperatureF(float* temperature);
unsigned int GetTemperature(int* temperature);
unsigned int GetTemperatureRange(int* mintemp,int* maxtemp);
unsigned int GetTemperatureStatus(float *SensorTemp, float *TargetTemp, float *AmbientTemp, float *CoolerVolts);
unsigned int GetTotalNumberImagesAcquired (at_32* index);
unsigned int GetVerticalSpeed(int index, int* speed);
unsigned int GetVirtualDMAAddress(void** Address1, void** Address2);
unsigned int GetVSSpeed(int index, float* speed);
unsigned int GPIBReceive(int id, short address, char* text, int size);
unsigned int GPIBSend(int id, short address, char* text);
unsigned int I2CBurstRead(unsigned char i2cAddress, at_32 nBytes, unsigned char* data);
unsigned int I2CBurstWrite(unsigned char i2cAddress, at_32 nBytes, unsigned char* data);
unsigned int I2CRead(unsigned char deviceID, unsigned char intAddress, unsigned char* pdata);
unsigned int I2CReset(void);
unsigned int I2CWrite(unsigned char deviceID, unsigned char intAddress, unsigned char data);
unsigned int IdAndorDll();
unsigned int InAuxPort(int port, int* state);
unsigned int Initialize(char * dir);	 //	read ini file to get head and card
unsigned int InitializeDevice(char * dir);
unsigned int IsInternalMechanicalShutter(int *InternalShutter);
unsigned int IsPreAmpGainAvailable(int channel, int amplifier, int index, int pa, int* status);
unsigned int IsTriggerModeAvailable(int _itriggerMode);
unsigned int Merge(const at_32* array, at_32 nOrder, at_32 nPoint, at_32 nPixel, float* coeff, at_32 fit, at_32* output, float* start, float* step);
unsigned int OutAuxPort(int port, int state);
unsigned int PrepareAcquisition();
unsigned int SaveAsBmp(char* path, char* palette, at_32 ymin, at_32 ymax);
unsigned int SaveAsCommentedSif(char* path, char* comment);
unsigned int SaveAsEDF(char* _szPath, int _iMode);
unsigned int SaveAsFITS(char* szFileTitle, int type);
unsigned int SaveAsRaw(char* szFileTitle, int type);
unsigned int SaveAsSif(char* path);
//unsigned int SaveAsSPC(char* path);
unsigned int SaveAsTiff(char* path, char* palette, int position, int type);
unsigned int SaveEEPROMToFile(char *cFileName);
unsigned int SaveToClipBoard(char* palette);
unsigned int SelectDevice(int devNum);
unsigned int SendSoftwareTrigger(void);
unsigned int SetAccumulationCycleTime(float time);
//unsigned int SetAcqStatusEvent(void* event);
unsigned int SetAcquisitionMode(int mode);
unsigned int SetAcquisitionType(int type);
unsigned int SetADChannel(int channel);
unsigned int SetBackground(at_32* array, at_u32 size);
unsigned int SetBaselineClamp(int state);
unsigned int SetBaselineOffset(int offset);
unsigned int SetCameraStatusEnable(unsigned long Enable);
unsigned int SetComplexImage(int numAreas, int* areas);
unsigned int SetCoolerMode(int mode);
unsigned int SetCropMode(int active, int cropheight, int reserved);
unsigned int SetCurrentCamera(at_32 cameraHandle);
unsigned int SetCustomTrackHBin(int bin);
unsigned int SetDataType(int type);
unsigned int SetDDGAddress(unsigned char t0, unsigned char t1, unsigned char t2, unsigned char tt, unsigned char address);
unsigned int SetDDGGain(int gain);
unsigned int SetDDGGateStep(double step);
unsigned int SetDDGInsertionDelay(int state);
unsigned int SetDDGIntelligate(int state);
unsigned int SetDDGIOC(int state);
unsigned int SetDDGIOCFrequency(double frequency);
unsigned int SetDDGIOCNumber(unsigned long number_pulses);
unsigned int SetDDGTimes(double t0, double t1, double t2);
unsigned int SetDDGTriggerMode(int mode);
unsigned int SetDDGVariableGateStep(int mode, double p1, double p2);
unsigned int SetDelayGenerator(int board, short address, int type);
unsigned int SetDMAParameters(int MaxImagesPerDMA, float SecondsPerDMA);
//unsigned int SetDriverEvent(void* event);
unsigned int SetEMAdvanced(int state);
unsigned int SetEMCCDGain(int gain);
unsigned int SetEMClockCompensation(int EMClockCompensationFlag);
unsigned int SetEMGainMode(int mode);
unsigned int SetExposureTime(float time);
unsigned int SetFanMode(int mode);
unsigned int SetFastExtTrigger(int mode);
unsigned int SetFastKineticsEx(int exposedRows, int seriesLength, float time, int mode, int hbin, int vbin, int offset);
unsigned int SetFastKinetics(int exposedRows, int seriesLength, float time, int mode, int hbin, int vbin);
unsigned int SetFilterMode(int mode);
unsigned int SetFilterParameters(int width, float sensitivity, int range, float accept, int smooth, int noise);
unsigned int SetFKVShiftSpeed(int index);
unsigned int SetFPDP(int state);
unsigned int SetFrameTransferMode(int mode);
unsigned int SetFullImage(int hbin, int vbin);
unsigned int SetFVBHBin(int bin);
unsigned int SetGain(int gain);
unsigned int SetGate(float delay, float width, float step);
unsigned int SetGateMode(int gatemode);
unsigned int SetHighCapacity(int state);
unsigned int SetHorizontalSpeed(int index);
unsigned int SetHSSpeed(int type, int index);
unsigned int SetImage(int hbin, int vbin, int hstart, int hend, int vstart, int vend);
unsigned int SetImageFlip(int iHFlip, int iVFlip);
unsigned int SetImageRotate(int iRotate);
unsigned int SetIsolatedCropMode (int active, int cropheight, int cropwidth, int vbin, int hbin);
unsigned int SetKineticCycleTime(float time);
unsigned int SetMCPGating(int gating);
unsigned int SetMessageWindow(int wnd);
unsigned int SetMultiTrackHBin(int bin);
unsigned int SetMultiTrack(int number, int height, int offset,int* bottom,int* gap);
unsigned int SetNextAddress16(at_32* data, at_32 lowAdd, at_32 highAdd, at_32 len, at_32 physical);
unsigned int SetNextAddress(at_32* data, at_32 lowAdd, at_32 highAdd, at_32 len, at_32 physical);
unsigned int SetNumberAccumulations(int number);
unsigned int SetNumberKinetics(int number);
unsigned int SetOutputAmplifier(int type);
unsigned int SetPCIMode(int mode,int value);
unsigned int SetPhotonCounting(int state);
unsigned int SetPhotonCountingThreshold(at_32 counts);
unsigned int SetPixelMode(int bitdepth, int colormode);
unsigned int SetPreAmpGain(int index);
unsigned int SetRandomTracks(int numTracks, int* areas);
unsigned int SetReadMode(int mode);
unsigned int SetRegisterDump(int mode);
unsigned int SetRingExposureTimes(int numTimes, float *times);
//unsigned int SetSaturationEvent(HANDLE event);
unsigned int SetShutter(int type, int mode, int closingtime, int openingtime);
unsigned int SetShutterEx(int type, int mode, int closingtime, int openingtime, int ext_mode);
unsigned int SetSifComment(char* comment);
unsigned int SetSingleTrackHBin(int bin);
unsigned int SetSingleTrack(int centre, int height);
unsigned int SetSpool(int active, int method, char* path, int framebuffersize);
unsigned int SetStorageMode(at_32 mode);
unsigned int SetTemperature(int temperature);
unsigned int SetTriggerMode(int mode);
//unsigned int SetUserEvent(at_u32 event);
unsigned int SetUSGenomics(at_32 width, at_32 height);
unsigned int SetVerticalRowBuffer(int rows);
unsigned int SetVerticalSpeed(int index);
unsigned int SetVirtualChip(int state);
unsigned int SetVSAmplitude(int index);
unsigned int SetVSSpeed(int index);
unsigned int ShutDown(void);
unsigned int StartAcquisition(void);
unsigned int UnMapPhysicalAddress();
unsigned int WaitForAcquisition();
unsigned int WaitForAcquisitionByHandle(at_32 cameraHandle);
unsigned int WaitForAcquisitionByHandleTimeOut(long cameraHandle, int _iTimeOutMs);
unsigned int WaitForAcquisitionTimeOut(int _iTimeOutMs);
unsigned int WhiteBalance(unsigned short* _wRed, unsigned short* _wGreen, unsigned short* _wBlue, float *_fRelR, float *_fRelB, WhiteBalanceInfo *_info);

#define DRV_ERROR_CODES	20001
#define DRV_SUCCESS	20002
#define DRV_VXDNOTINSTALLED	20003
#define DRV_ERROR_SCAN	20004
#define DRV_ERROR_CHECK_SUM	20005
#define DRV_ERROR_FILELOAD	20006
#define DRV_UNKNOWN_FUNCTION	20007
#define DRV_ERROR_VXD_INIT	20008
#define DRV_ERROR_ADDRESS	20009
#define DRV_ERROR_PAGELOCK	20010
#define DRV_ERROR_PAGEUNLOCK	20011
#define DRV_ERROR_BOARDTEST	20012
#define DRV_ERROR_ACK	20013
#define DRV_ERROR_UP_FIFO	20014
#define DRV_ERROR_PATTERN	20015

#define DRV_ACQUISITION_ERRORS	20017
#define DRV_ACQ_BUFFER	20018
#define DRV_ACQ_DOWNFIFO_FULL	20019
#define DRV_PROC_UNKONWN_INSTRUCTION	20020
#define DRV_ILLEGAL_OP_CODE	20021
#define DRV_KINETIC_TIME_NOT_MET	20022
#define DRV_ACCUM_TIME_NOT_MET	20023
#define DRV_NO_NEW_DATA	20024
#define DRV_SPOOLERROR 20026
#define KERN_MEM_ERROR					20025
#define DRV_SPOOLSETUPERROR 20027
#define DRV_FILESIZELIMITERROR 20028
#define DRV_ERROR_FILESAVE 20029

#define DRV_TEMPERATURE_CODES	20033
#define DRV_TEMPERATURE_OFF	20034
#define DRV_TEMPERATURE_NOT_STABILIZED	20035
#define DRV_TEMPERATURE_STABILIZED	20036
#define DRV_TEMPERATURE_NOT_REACHED	20037
#define DRV_TEMPERATURE_OUT_RANGE	20038
#define DRV_TEMPERATURE_NOT_SUPPORTED	20039
#define DRV_TEMPERATURE_DRIFT	20040

#define DRV_TEMP_CODES	20033
#define DRV_TEMP_OFF	20034
#define DRV_TEMP_NOT_STABILIZED	20035
#define DRV_TEMP_STABILIZED	20036
#define DRV_TEMP_NOT_REACHED	20037
#define DRV_TEMP_OUT_RANGE	20038
#define DRV_TEMP_NOT_SUPPORTED	20039
#define DRV_TEMP_DRIFT	20040

#define DRV_GENERAL_ERRORS	20049
#define DRV_INVALID_AUX	20050
#define DRV_COF_NOTLOADED	20051
#define DRV_FPGAPROG 20052
#define DRV_FLEXERROR 20053
#define DRV_GPIBERROR 20054
#define DRV_EEPROMVERSIONERROR 20055

#define DRV_DATATYPE	20064
#define DRV_DRIVER_ERRORS	20065
#define DRV_P1INVALID	20066
#define DRV_P2INVALID	20067
#define DRV_P3INVALID	20068
#define DRV_P4INVALID	20069
#define DRV_INIERROR	20070
#define DRV_COFERROR	20071
#define DRV_ACQUIRING	20072
#define DRV_IDLE	20073
#define DRV_TEMPCYCLE	20074
#define DRV_NOT_INITIALIZED 20075
#define DRV_P5INVALID	20076
#define DRV_P6INVALID	20077
#define DRV_INVALID_MODE	20078
#define DRV_INVALID_FILTER 20079

#define DRV_I2CERRORS	20080
#define DRV_I2CDEVNOTFOUND	20081
#define DRV_I2CTIMEOUT	20082
#define DRV_P7INVALID	20083
#define DRV_USBERROR	20089
#define DRV_IOCERROR 20090
#define DRV_VRMVERSIONERROR 20091
#define DRV_USB_INTERRUPT_ENDPOINT_ERROR 20093
#define DRV_RANDOM_TRACK_ERROR 20094
#define DRV_INVALID_TRIGGER_MODE 20095
#define DRV_LOAD_FIRMWARE_ERROR 20096
#define DRV_DIVIDE_BY_ZERO_ERROR 20097
#define DRV_INVALID_RINGEXPOSURES 20098

#define DRV_ERROR_NOCAMERA 20990
#define DRV_NOT_SUPPORTED 20991
#define DRV_NOT_AVAILABLE 20992

#define DRV_ERROR_MAP 20115
#define DRV_ERROR_UNMAP 20116
#define DRV_ERROR_MDL 20117
#define DRV_ERROR_UNMDL 20118
#define DRV_ERROR_BUFFSIZE 20119
#define DRV_ERROR_NOHANDLE 20121

#define DRV_GATING_NOT_AVAILABLE 20130
#define DRV_FPGA_VOLTAGE_ERROR   20131

#define DRV_MSTIMINGS_ERROR 20156


#define AC_ACQMODE_SINGLE 1
#define AC_ACQMODE_VIDEO 2
#define AC_ACQMODE_ACCUMULATE 4
#define AC_ACQMODE_KINETIC 8
#define AC_ACQMODE_FRAMETRANSFER 16
#define AC_ACQMODE_FASTKINETICS 32

#define AC_READMODE_FULLIMAGE 1
#define AC_READMODE_SUBIMAGE 2
#define AC_READMODE_SINGLETRACK 4
#define AC_READMODE_FVB 8
#define AC_READMODE_MULTITRACK 16
#define AC_READMODE_RANDOMTRACK 32

#define AC_TRIGGERMODE_INTERNAL 1
#define AC_TRIGGERMODE_EXTERNAL 2
#define AC_TRIGGERMODE_EXTERNAL_FVB_EM 4
#define AC_TRIGGERMODE_CONTINUOUS 8
#define AC_TRIGGERMODE_EXTERNALSTART 16
#define AC_TRIGGERMODE_BULB 32

#define AC_CAMERATYPE_PDA 0
#define AC_CAMERATYPE_IXON 1
#define AC_CAMERATYPE_ICCD 2
#define AC_CAMERATYPE_EMCCD 3
#define AC_CAMERATYPE_CCD 4
#define AC_CAMERATYPE_ISTAR 5
#define AC_CAMERATYPE_VIDEO 6
#define AC_CAMERATYPE_IDUS 7
#define AC_CAMERATYPE_NEWTON 8
#define AC_CAMERATYPE_SURCAM 9
#define AC_CAMERATYPE_USBISTAR 10
#define AC_CAMERATYPE_LUCA 11
#define AC_CAMERATYPE_RESERVED 12
#define AC_CAMERATYPE_IKON 13
#define AC_CAMERATYPE_INGAAS 14

#define AC_PIXELMODE_8BIT  1
#define AC_PIXELMODE_14BIT 2
#define AC_PIXELMODE_16BIT 4
#define AC_PIXELMODE_32BIT 8

#define AC_PIXELMODE_MONO 0
#define AC_PIXELMODE_RGB (1 << 16)
#define AC_PIXELMODE_CMY (2 << 16)

#define AC_SETFUNCTION_VREADOUT 1
#define AC_SETFUNCTION_HREADOUT 2
#define AC_SETFUNCTION_TEMPERATURE 4
#define AC_SETFUNCTION_GAIN 8
#define AC_SETFUNCTION_ICCDGAIN 8
#define AC_SETFUNCTION_EMCCDGAIN 16
#define AC_SETFUNCTION_BASELINECLAMP 32
#define AC_SETFUNCTION_VSAMPLITUDE 64
#define AC_SETFUNCTION_HIGHCAPACITY 128
#define AC_SETFUNCTION_BASELINEOFFSET 256
#define AC_SETFUNCTION_PREAMPGAIN 512
#define AC_SETFUNCTION_CROPMODE 1024
#define AC_SETFUNCTION_DMAPARAMETERS 2048

#define AC_GETFUNCTION_TEMPERATURE 1
#define AC_GETFUNCTION_TARGETTEMPERATURE 2
#define AC_GETFUNCTION_TEMPERATURERANGE 4
#define AC_GETFUNCTION_DETECTORSIZE 8
#define AC_GETFUNCTION_GAIN 16
#define AC_GETFUNCTION_ICCDGAIN 16
#define AC_GETFUNCTION_EMCCDGAIN 32

#define AC_FEATURES_POLLING 1
#define AC_FEATURES_EVENTS 2
#define AC_FEATURES_SPOOLING 4
#define AC_FEATURES_SHUTTER 8
#define AC_FEATURES_SHUTTEREX 16
#define AC_FEATURES_EXTERNAL_I2C 32
#define AC_FEATURES_SATURATIONEVENT 64
#define AC_FEATURES_FANCONTROL 128
#define AC_FEATURES_MIDFANCONTROL 256
#define AC_FEATURES_TEMPERATUREDURINGACQUISITION 512

#define AC_EMGAIN_8BIT 1
#define AC_EMGAIN_12BIT 2
#define AC_EMGAIN_LINEAR12 4
#define AC_EMGAIN_REAL12 8

#ifdef __cplusplus
}
#endif

#endif
