//========================================================================
// SBIGDEF.H
// Contains the function prototypes and enumerated constants for the 
// linux driver.
//
// This supports the following devices:
// ST-7E/8E/9E/10E
// ST-5C/237/237A (PixCel255/237)
// ST-1K, AO-7
//
// Version 4.0 - September 17, 2001
// Copyright (C) 1995-2001 - Santa Barbara Instrument Group
//
// History:
//
// 2002-03-31:
// Adapted for Linux character device drivers by Jan Soldan, 
// soldan@obs.unige.ch
//
// 2003-12-22: 
// Changed values like FULL_CLEAR_TIME, NIBBLE_TIMEOUT, IDLE_STATE_DELAY, 
// etc. for Linux timing. These values are mostly defined in milliseconds 
// while under Windows in number clock ticks. See lines around 1247... 

//========================================================================

#ifndef _SBIGDEF_H_	      
#define _SBIGDEF_H_

typedef unsigned short MY_LOGICAL;

/*
	Supported Camera Commands

	These are the commands supported by the driver.
	They are prefixed by CC_ to designate them as
	camera commands and avoid conflicts with other
	enums.

	Some of the commands are marked as SBIG use only
	and have been included to enhance testability
	of the driver for SBIG.
*/

typedef enum {
 CC_START_EXPOSURE = 1,
 CC_END_EXPOSURE,
 CC_READOUT_LINE,
 CC_DUMP_LINES,
 CC_SET_TEMPERATURE_REGULATION,
 CC_QUERY_TEMPERATURE_STATUS,
 CC_ACTIVATE_RELAY,
 CC_PULSE_OUT,
 CC_ESTABLISH_LINK,
 CC_GET_DRIVER_INFO,

 CC_GET_CCD_INFO,
 CC_QUERY_COMMAND_STATUS,
 CC_MISCELLANEOUS_CONTROL,
 CC_READ_SUBTRACT_LINE,
 CC_UPDATE_CLOCK,
 CC_READ_OFFSET,
 CC_OPEN_DRIVER,
 CC_CLOSE_DRIVER,
 CC_TX_SERIAL_BYTES,
 CC_GET_SERIAL_STATUS,
 
 CC_AO_TIP_TILT,
 CC_AO_SET_FOCUS,
 CC_AO_DELAY,
 CC_GET_TURBO_STATUS,
 CC_END_READOUT,
 CC_GET_US_TIMER,
 CC_OPEN_DEVICE,
 CC_CLOSE_DEVICE,
 CC_SET_IRQL,
 CC_GET_IRQL,
 
 CC_GET_LINE,
 CC_GET_LINK_STATUS,
 CC_GET_DRIVER_HANDLE,
 CC_SET_DRIVER_HANDLE,
 CC_START_READOUT,
 CC_GET_ERROR_STRING,
 CC_SET_DRIVER_CONTROL, 
 CC_GET_DRIVER_CONTROL,
 CC_USB_AD_CONTROL,

 /* SBIG internal commands */
 CC_SEND_BLOCK = 90,
 CC_SEND_BYTE,
 CC_GET_BYTE,
 CC_SEND_AD,
 CC_GET_AD,
 CC_CLOCK_AD,
 CC_SYSTEM_TEST,
 CC_GET_DRIVER_OPTIONS, 
 CC_SET_DRIVER_OPTIONS,
 CC_LAST_COMMAND,

 /* JSIMAGE internal commands */
 CC_START_GRAB_MODE = 1000,
 CC_STOP_GRAB_MODE,
 CC_FAN_CHANGED,
 CC_TERMINATE,
 CC_UNDEFINED
}PAR_COMMAND;

/*
	Return Error Codes

	These are the error codes returned by the driver
	function.  They are prefixed with CE_ to designate
	them as camera errors.
*/

#ifndef CE_ERROR_BASE
#define CE_ERROR_BASE 1
#endif

typedef enum {
 CE_NO_ERROR,
 CE_CAMERA_NOT_FOUND = CE_ERROR_BASE,
 CE_EXPOSURE_IN_PROGRESS,
 CE_NO_EXPOSURE_IN_PROGRESS,
 CE_UNKNOWN_COMMAND,
 CE_BAD_CAMERA_COMMAND,
 CE_BAD_PARAMETER,
 CE_TX_TIMEOUT,
 CE_RX_TIMEOUT,
 CE_NAK_RECEIVED,
 CE_CAN_RECEIVED,

 CE_UNKNOWN_RESPONSE,
 CE_BAD_LENGTH,
 CE_AD_TIMEOUT,
 CE_KBD_ESC,
 CE_CHECKSUM_ERROR,
 CE_EEPROM_ERROR,
 CE_SHUTTER_ERROR,
 CE_UNKNOWN_CAMERA,
 CE_DRIVER_NOT_FOUND,
 CE_DRIVER_NOT_OPEN,

 CE_DRIVER_NOT_CLOSED,
 CE_SHARE_ERROR,
 CE_TCE_NOT_FOUND,
 CE_AO_ERROR,
 CE_ECP_ERROR,
 CE_MEMORY_ERROR,
 CE_DEVICE_NOT_FOUND,
 CE_DEVICE_NOT_OPEN,
 CE_DEVICE_NOT_CLOSED,
 CE_DEVICE_NOT_IMPLEMENTED,

 CE_DEVICE_DISABLED,
 CE_OS_ERROR,
 CE_SOCK_ERROR,
 CE_SERVER_NOT_FOUND,
 CE_NEXT_ERROR,

 CE_IOCTL = 1000,
 CE_STOP_GRAB_MODE

}PAR_ERROR;

/*
	Camera Command State Codes

	These are the return status codes for the Query
	Command Status command.  They are prefixed with
	CS_ to designate them as camera status.

*/
typedef enum{
 CS_IDLE,
 CS_IN_PROGRESS,
 CS_INTEGRATING,
 CS_INTEGRATION_COMPLETE
}PAR_COMMAND_STATUS;

#define CS_PULSE_IN_ACTIVE       0x8000
#define CS_WAITING_FOR_TRIGGER   0x8000

#define CONFIG_LOW_IMAGER_IS_311 2
/*
 Miscellaneous Enumerated Constants

 ABG_STATE7             - Passed to Start Exposure Command
 unsigned char          - General purpose type
 DRIVER_REQUEST         - Used with Get Driver Info command
 CCD_REQUEST            - Used with Imaging commands to specify CCD
 CCD_INFO_REQUEST       - Used with Get CCD Info Command
 PORT                   - Used with Establish Link Command
 CAMERA_TYPE            - Returned by Establish Link and Get CCD Info commands
 SHUTTER_COMMAND        - Used with Start Exposure and Misc. Control Commands
 SHUTTER_STATE7         - Used with Start Exposure and Misc. Control Commands
 TEMPERATURE_REGULATION - Used with Enable Temperature Regulation
 LED_STATE              - Used with the Miscellaneous Control Command
 FILTER_COMMAND,        - Used with the Miscellaneous	Control Command
 FILTER_STATE           - Used with the Miscellaneous	Control Command
 AD_SIZE, FILTER_TYPE   - Used with the GetCCDInfo3 Command
 AO_FOCUS_COMMAND       - Used with the AO Set Focus Command
 SBIG_DEVICE_TYPE       - Used with Open Device Command
*/

typedef enum{ 
 MC_START_EXPOSURE=0, 
 MC_END_EXPOSURE, 
 MC_REGULATE_TEMP,
 MC_TEMP_STATUS, 
 MC_RELAY, 
 MC_PULSE, 
 MC_GET_VERSION, 
 MC_EEPROM,
 MC_MISC_CONTROL, 
 MC_STATUS,
 MC_SYSTEM_TEST, 
 MC_TX_BYTES,
 MC_CONTROL_CCD, 
 MC_COMMAND_13, 
 MC_SYSTEM,
 MC_READOUT, 
 MC_REGULATE_TEMP2 = 128 + MC_REGULATE_TEMP 
}MICRO_COMMAND;

typedef enum{ 
 SYS_READ_INT, 
 SYS_WRITE_INT, 
 SYS_READ_EXT, 
 SYS_WRITE_EXT,
 SYS_GET_ROM_SUM 
}SYSTEM_SUBCOMMAND;

typedef enum{
 RS_DIG_ROW, 
 RS_DLP_ROW, 
 RS_DL_ROW, 
 RS_DLP_ROWS,
 RS_DUMP_FIFO, 
 RS_DL_SETUP, 
 RS_DUMP_ROWS, 
 RS_CLEAR_CCD,
 RS_SET_VDD, 
 RS_WRITE_AD, 
 RS_DLPP_ROWS 
}READOUT_SUBCOMMAND;

typedef enum{
 ABG_LOW7,
 ABG_CLK_LOW7,
 ABG_CLK_MED7,
 ABG_CLK_HI7
}ABG_STATE7;

typedef enum{
 DRIVER_STD,
 DRIVER_EXTENDED
}DRIVER_REQUEST;

typedef enum{
 CCD_IMAGING,
 CCD_TRACKING
}CCD_REQUEST;

typedef enum{
 CCD_INFO_IMAGING,
 CCD_INFO_TRACKING,
 CCD_INFO_EXTENDED,
 CCD_INFO_EXTENDED_5C,
 CCD_INFO_EXTENDED2_IMAGING,
 CCD_INFO_EXTENDED2_TRACKING
}CCD_INFO_REQUEST;

typedef enum{
 ABG_NOT_PRESENT,
 ABG_PRESENT
}IMAGING_ABG;

typedef enum{
 BR_AUTO,
 BR_9600,
 BR_19K,
 BR_38K,
 BR_57K,
 BR_115K
}PORT_RATE;

typedef enum{
 ST7_CAMERA = 4,
 ST8_CAMERA,
 ST5C_CAMERA,
 TCE_CONTROLLER,
 ST237_CAMERA,
 STK_CAMERA,
 ST9_CAMERA,
 STV_CAMERA,
 ST10_CAMERA,
 ST1K_CAMERA,
 ST2K_CAMERA
}CAMERA_TYPE;

typedef enum{
 SC_LEAVE_SHUTTER,
 SC_OPEN_SHUTTER,
 SC_CLOSE_SHUTTER,
 SC_INITIALIZE_SHUTTER
}SHUTTER_COMMAND;

typedef enum{
 SS_OPEN,
 SS_CLOSED,
 SS_OPENING,
 SS_CLOSING
}SHUTTER_STATE7;

typedef enum{
 REGULATION_OFF,
 REGULATION_ON,
 REGULATION_OVERRIDE,
 REGULATION_FREEZE,
 REGULATION_UNFREEZE,
 REGULATION_ENABLE_AUTOFREEZE,
 REGULATION_DISABLE_AUTOFREEZE
}TEMPERATURE_REGULATION;

#define REGULATION_FROZEN_MASK 0x8000

typedef enum{
 LED_OFF,
 LED_ON,
 LED_BLINK_LOW,
 LED_BLINK_HIGH
}LED_STATE;

typedef enum{
 FILTER_LEAVE,
 FILTER_SET_1,
 FILTER_SET_2,
 FILTER_SET_3,
 FILTER_SET_4,
 FILTER_SET_5,
 FILTER_STOP,
 FILTER_INIT
}FILTER_COMMAND;

typedef enum{
 FS_MOVING,
 FS_AT_1,
 FS_AT_2,
 FS_AT_3,
 FS_AT_4,
 FS_AT_5,
 FS_UNKNOWN
}FILTER_STATE;

typedef enum{
 AD_UNKNOWN,
 AD_12_BITS,
 AD_16_BITS
}AD_SIZE;

typedef enum{
 FW_UNKNOWN,
 FW_EXTERNAL,
 FW_VANE,
 FW_FILTER_WHEEL
}FILTER_TYPE;

typedef enum{
 AOF_HARD_CENTER,
 AOF_SOFT_CENTER,
 AOF_STEP_IN,
 AOF_STEP_OUT
}AO_FOCUS_COMMAND;

typedef enum{
 DEV_NONE,
 DEV_LPT1,
 DEV_LPT2,
 DEV_LPT3,
 DEV_USB = 0x7F00,
 DEV_ETH
}SBIG_DEVICE_TYPE;

typedef enum{ 
 DCP_USB_FIFO_ENABLE, 
 DCP_CALL_JOURNAL_ENABLE,
 DCP_IVTOH_RATIO, 
 DCP_USB_FIFO_SIZE, 
 DCP_USB_DRIVER, 
 DCP_LAST 
}DRIVER_CONTROL_PARAM;

typedef enum{ 
 USB_AD_IMAGING_GAIN, 
 USB_AD_IMAGING_OFFSET, 
 USB_AD_TRACKING_GAIN,
 USB_AD_TRACKING_OFFSET 
}USB_AD_CONTROL_COMMAND;

typedef enum{ 
 USBD_SBIGE, 
 USBD_SBIGI, 
 USBD_SBIGJ, 
 USBD_NEXT 
}ENUM_USB_DRIVER;


typedef enum{
 EXP_IDLE,
 EXP_IN_PROGRESS = 2,
 EXP_COMPLETE
}EXPOSURE_STATE;

typedef enum{
 IMAGING_CCD,
 TRACKING_CCD
}CCD;

typedef enum{
 CCD_THERMISTOR,
 AMBIENT_THERMISTOR
}THERMISTOR_TYPE;

/*
	General Purpose Flags
*/

#define END_SKIP_DELAY 0x8000
/*
set in ccd parameter of EndExposure command to skip synchronization
delay - Use this to increase the rep rate when taking darks to later
be subtracted from SC_LEAVE_SHUTTER exposures such as when tracking
and imaging
*/

#define START_SKIP_VDD	0x8000		
/* set in ccd parameter of StartExposure command to skip lowering
Imaging CCD Vdd during integration. - Use this to increase the rep
rate when you don't care about glow in the upper-left corner of the
imaging CCD */

#define EXP_WAIT_FOR_TRIGGER_IN 0x80000000
/* set in exposureTime to wait for trigger in pulse */

#define EXP_SEND_TRIGGER_OUT    0x40000000
/* set in exposureTime to send  trigger out Y- */

#define EXP_LIGHT_CLEAR         0x20000000	
/* set to do light clear of the CCD (Kaiser only) */

#define EXP_TIME_MASK           0x00FFFFFF
/* mask with exposure time to remove flags */


/*	Defines */
#define MIN_ST7_EXPOSURE    12  /* Minimum exposure is 12/100ths second */

/*
	Command Parameter and Results Structs

	Make sure you set your compiler for byte structure alignment
	as that is how the driver was built.
*/
typedef struct{
 unsigned char subCommand;
 unsigned char len;
 short         destAddress;
 unsigned char data[8];
}SystemMemoryParams;

typedef struct{
 unsigned char subCommand;
 unsigned char byte1;
 unsigned char byte2;
}ReadoutWriteADParams;

typedef struct{
 unsigned char data[8];
}SystemMemoryResults;

typedef struct{
 unsigned short ccd;          /* CCD_REQUEST     */
 unsigned long  exposureTime;
 unsigned short abgState;     /* ABG_STATE7      */
 unsigned short openShutter;  /* SHUTTER_COMMAND */
}StartExposureParams;

typedef struct{
 unsigned short ccd;
}EndExposureParams;

typedef struct{
 unsigned short ccd;
 unsigned short readoutMode;
 unsigned short pixelStart;
 unsigned short pixelLength;
}ReadoutLineParams;

typedef struct{
 unsigned short ccd;
 unsigned short readoutMode;
 unsigned short lineLength;
}DumpLinesParams;

typedef struct{
 unsigned short ccd;
}EndReadoutParams;

typedef struct{
 unsigned short ccd;
 unsigned short readoutMode;
 unsigned short top;
 unsigned short left;
 unsigned short height;
 unsigned short width;
}StartReadoutParams;

typedef struct{
 unsigned short regulation; /* TEMPERATURE_REGULATION */
 unsigned short ccdSetpoint;
}SetTemperatureRegulationParams;

typedef struct{
 MY_LOGICAL     enabled;
 unsigned short ccdSetpoint;
 unsigned short power;
 unsigned short ccdThermistor;
 unsigned short ambientThermistor;
}QueryTemperatureStatusResults;

typedef struct{
 unsigned short tXPlus;
 unsigned short tXMinus;
 unsigned short tYPlus;
 unsigned short tYMinus;
}ActivateRelayParams;

typedef struct{
 unsigned short numberPulses;
 unsigned short pulseWidth;
 unsigned short pulsePeriod;
}PulseOutParams;

typedef struct{
 unsigned short dataLength;
 unsigned char  data[256];
}TXSerialBytesParams;

typedef struct{
 unsigned short bytesSent;
}TXSerialBytesResults;

typedef struct{
 MY_LOGICAL clearToCOM;
}GetSerialStatusResults;

typedef struct{
 unsigned short sbigUseOnly;
}EstablishLinkParams;

typedef struct{
 unsigned short cameraType; /* CAMERA_TYPE */
}EstablishLinkResults;

typedef struct {
 unsigned short request; /* DRIVER_REQUEST */
}GetDriverInfoParams;

typedef struct{
 unsigned short version;
 char           name[64];
 unsigned short maxRequest;
}GetDriverInfoResults0;

typedef struct{
 unsigned short request; /* CCD_INFO_REQUEST */
}GetCCDInfoParams;

typedef struct{
 unsigned short mode;
 unsigned short width;
 unsigned short height;
 unsigned short gain;
 unsigned long  pixel_width;
 unsigned long  pixel_height;
}READOUT_INFO;

typedef struct{
 unsigned short firmwareVersion;
 unsigned short cameraType; /* CAMERA_TYPE */
 char           name[64];
 unsigned short readoutModes;
 struct{
   unsigned short mode;
   unsigned short width;
   unsigned short height;
   unsigned short gain;
   unsigned long  pixelWidth;
   unsigned long  pixelHeight;
 }readoutInfo[20];
}GetCCDInfoResults0;

typedef struct{
 unsigned short badColumns;
 unsigned short columns[4];
 unsigned short imagingABG; /* IMAGING_ABG */
 char           serialNumber[10];
}GetCCDInfoResults2;

typedef struct{
 unsigned short adSize;     /* AD_SIZE     */
 unsigned short filterType; /* FILTER_TYPE */
}GetCCDInfoResults3;

typedef struct{
 unsigned short capabilitiesBits;
 unsigned short dumpExtra;
}GetCCDInfoResults4;

typedef struct{
 unsigned short command;
}QueryCommandStatusParams;

typedef struct{
 unsigned short status;
}QueryCommandStatusResults;

typedef struct{
 MY_LOGICAL     fanEnable;
 unsigned short shutterCommand; /* SHUTTER_COMMAND */
 unsigned short ledState;       /* LED_STATE       */
}MiscellaneousControlParams;

typedef struct{
 unsigned short ccd;
}ReadOffsetParams;

typedef struct{
 unsigned short offset;
}ReadOffsetResults;

typedef struct{
 unsigned short xDeflection;
 unsigned short yDeflection;
}AOTipTiltParams;

typedef struct{
 unsigned short focusCommand; /* AO_FOCUS_COMMAND */
}AOSetFocusParams;

typedef struct{
 unsigned long  delay;
}AODelayParams;

typedef struct{
 MY_LOGICAL turboDetected;
}GetTurboStatusResults;

typedef struct{
 unsigned short	deviceType;    /* SBIG_DEVICE_TYPE, specifies LPT, Ethernet, etc */
 unsigned short	lptBaseAddress;/* DEV_LPTN: Windows 9x Only, Win NT uses deviceSelect */
 unsigned long	ipAddress;     /* DEV_ETH:  Ethernet address */
}OpenDeviceParams;

typedef struct{
 unsigned short level;
}SetIRQLParams;

typedef struct{
 unsigned short level;
}GetIRQLResults;

typedef struct{
 MY_LOGICAL	linkEstablished;
 unsigned short baseAddress;
 unsigned short cameraType;  /* CAMERA_TYPE */
 unsigned long  comTotal;
 unsigned long  comFailed;
}GetLinkStatusResults;

typedef struct{
 unsigned long count;
}GetUSTimerResults;

typedef struct{
 unsigned short port;
 unsigned short length;
 unsigned char  *source;
}SendBlockParams;

typedef struct{
 unsigned short port;
 unsigned short data;
}SendByteParams;

typedef struct{
 unsigned short ccd;
 unsigned short readoutMode;
 unsigned short pixelStart;
 unsigned short pixelLength;
}ClockADParams;

typedef struct{
 unsigned short testClocks;
 unsigned short testMotor;
 unsigned short test5800;
}SystemTestParams;

typedef struct{
 unsigned short outLength;
 unsigned char  *outPtr;
 unsigned short inLength;
 unsigned char  *inPtr;
}SendSTVBlockParams;

typedef struct{
 unsigned short errorNo;
}GetErrorStringParams;

typedef struct{
 char errorString[64];
}GetErrorStringResults;

typedef struct{
 short handle;
}SetDriverHandleParams;

typedef struct{
 short handle;
}GetDriverHandleResults;

typedef struct{
 unsigned short /* DRIVER_CONTROL_PARAM */ controlParameter;
 unsigned long controlValue;
}SetDriverControlParams;

typedef struct{
 unsigned short /* DRIVER_CONTROL_PARAM */ controlParameter;
}GetDriverControlParams; 

typedef struct{
 unsigned long controlValue;
}GetDriverControlResults;

typedef struct{
 unsigned short /* USB_AD_CONTROL_COMMAND */ command;
 short data;
}USBADControlParams;

typedef struct{
 unsigned short  raiseIt;
 unsigned short  vddWasLow;
}IocSetVdd;

typedef struct{
 short cameraID; // CAMERA_TYPE
 short ccd;      // CCD_REQUEST
 short left;
 short len;
 short right;
 short horzBin;
 short vertBin;
 short clearWidth;
 short st237A;
 short st253;
}IOC_GET_PIXELS_PARAMS;

typedef struct{
 short cameraID;
 short onVertBin;
 short clearWidth;
}IOC_VCLOCK_CCD_PARAMS;

typedef struct{
 short cameraID;
 short width;
 short len;
 short vertBin;
 short vToHRatio;
 short st253;
}IOC_DUMP_LINES_PARAMS;

typedef struct{
 short cameraID;
 short height;
 short times;
}IOC_CLEAR_CCD_PARAMS;

typedef struct{
 short cameraID;
 short ccd;
 short left;
 short len;
 short right;
 short horzBin;
 short vertBin;
 short clearWidth;
 short st237A;
 short height;
}IOC_GET_AREA_PARAMS;

typedef struct {
 short baseAddress;
}IOC_SET_PORT_ADDRESS_PARAMS;

#define HSI   0x10	/*  hand shake input bit */
#define CIP   0x10	/*  conversion in progress bit */
#define CAN   0x18	/*  CAN response from Micro */
#define NAK   0x15	/*  NAK response from Micro */
#define ACK   0x06	/*  ACK response from Micro */

typedef enum{
 TRACKING_CLOCKS = 0x00,
 IMAGING_CLOCKS  = 0x10,
 MICRO_OUT       = 0x20,
 CONTROL_OUT     = 0x30,
 READOUT_CONTROL = 0x00,
 DEVICE_SELECT   = 0x60
}OUTPUT_REGISTER;

typedef enum{
 AD0     = 0x00,
 AD1     = 0x10,
 AD2     = 0x20,
 AD3_MDI = 0x30
}INPUT_REGISTER;

typedef enum{
 HSO             = 0x01,
 MICRO_SELECT    = 0x02,
 IMAGING_SELECT  = 0x00,
 TRACKING_SELECT = 0x04,
 MICRO_SYNC      = 0x04,
 AD_TRIGGER      = 0x08
}CONTROL_BITS;

typedef enum{
 V1_H   = 1,
 V2_H   = 2,
 TRG_H  = 4,
 IABG_M = 8
}IMAGING_CLOCK_BITS;

typedef enum{
 P1V_H = 1,
 P2V_H = 2,
 P1H_L = 8
}KI_CLOCK_BITS;

typedef enum{
 IAG_H  = 1,
 TABG_M = 2,
 BIN    = 4,
 CLR    = 8
}TRACKING_CLOCK_BITS;

typedef enum{
 KCLR  = 0,
 KBIN1 = 4,
 KBIN2 = 8,
 KBIN3 = 12
}KT_CLOCK_BITS;

typedef enum{
 /* IAG_H=1 */
 SAG_H = 2,
 SRG_H = 4,
 R1_D3 = 8
}CCD_CLOCK_BITS;

typedef enum{
 CLR_SELECT       = 1,
 R3_D1            = 2,
 MICRO_CLK_SELECT = 4,
 PLD_TRIGGER      = 8
}READOUT_CONTROL_BITS;

#define IMAGING_GAIN  0x0230			// e-/ADU in XX.XX BCD
#define IMAGING_GAIN9 0x0280			// ST-9 e-/ADU in XX.XX BCD
#define TRACKING_GAIN 0x0130			// e-/ADU in XX.XX BCD
#define ST5C_GAIN     0x0125			// e-/ADU in XX.XX BCD
#define ST237_GAIN    0x0420			// e-/ADU in XX.XX BCD
#define ST237A_GAIN   0x0072			// e-/ADU in XX.XX BCD
#define STK_HIGH_GAIN 0x0250			// e-/ADU in XX.XX BCD
#define STK_LOW_GAIN  0x1000			// e-/ADU in XX.XX BCD

typedef	struct{
 short  width;
 short  height;
 short  top;
 short  left;
 short  bottom;
 short  right;
 short  full_width;
 short  full_height;
 short  binning[10];
 unsigned short gain;
 unsigned short offset;
 unsigned short bcd_pixel_width;
 unsigned short bcd_pixel_height;
 unsigned short dump_extra;
 unsigned short storage_height;
}CCD_INFO;
//========================================================================
#define MAX_DIG_WIDTH   4096  //  maximum pixels to digitize, was 2184
#define SAT_THR_4K      4094  //  threshold above which pixels forced to sat
#define SAT_VAL_4K      4095  //  value sat pixels forced to
#define SAT_THR_8K      8190  //  threshold above which pixels forced to sat
#define SAT_VAL_8K      8191  //  value sat pixels forced to
#define SAT_THR_12K     12285 //  threshold above which pixels forced to sat
#define SAT_VAL_12K     12287 //  value sat pixels forced to
#define SAT_THR_16K     16380 //  threshold above which pixels forced to sat
#define SAT_VAL_16K     16383 //  value sat pixels forced to
#define SAT_THR_32K     32760 //  threshold above which pixels forced to sat
#define SAT_VAL_32K     32767 //  value sat pixels forced to
#define SAT_THR_48K     49140 //  threshold above which pixels forced to sat
#define SAT_VAL_48K     49151 //  value sat pixels forced to
#define SAT_THR_65K     63000 //  threshold above which pixels forced to sat
#define SAT_VAL_65K     65535 //  value sat pixels forced to
#define HOT_THR_4K	600   //  if hot is this high jimmy pixel
#define HOT_THR_8K	3500  //  if hot is this high jimmy pixel
#define HOT_THR_12K	5000  //  if hot is this high jimmy pixel
#define HOT_THR_32K	8000  //  if hot is this high jimmy pixel
#define HOT_THR_65K	8000  //  if hot is this high jimmy pixel
#define MAX_HP_WIDTH	657   //  max width to do hot pixels over

typedef struct{
 unsigned short checksum;
 unsigned char  version;
 unsigned char  model;
 unsigned char  abgType;
 unsigned char  badColumns;
 unsigned short trackingOffset;
 unsigned short trackingGain;
 unsigned short imagingOffset;
 unsigned short imagingGain;
 unsigned short columns[3];
 unsigned short configWord;
 unsigned char  serialNumber[10];
}EEPROMContents;

typedef enum{
 SHUTTER_ASIS,
 OPEN_SHUTTER,
 CLOSE_SHUTTER
}SHUTTER_CONTROL;

typedef enum{
 SHUTTER_OPEN,
 SHUTTER_CLOSED,
 SHUTTER_OPENING,
 SHUTTER_CLOSING
}SHUTTER_STATUS;

typedef enum{
 CCD_IDLE,
 PRE_SHUTTER,
 INTEGRTAING,
 POST_SHUTTER
}CCD_STATE;

typedef enum{
 LINK_NONE,
 LINK_LPT,
 LINK_USB,
 LINK_ETH
}LINK_TYPE;

typedef struct{
 unsigned short sum;
}GetROMSumResults;

typedef struct{
 LINK_TYPE        linkType;
 SBIG_DEVICE_TYPE deviceType;
 unsigned char  open;
 unsigned long  comTotal;
 unsigned long  comFailed;
 unsigned long  comPassed;
}LINK_INFO;

typedef struct{
 TEMPERATURE_REGULATION regulation;
 unsigned short ccdSetpoint;
 unsigned short preload;
}MicroTemperatureRegulationParams;

typedef struct{
 unsigned char  freezeTE;
 unsigned char  lowerPGain;
}MicroTemperatureSpecialParams;

typedef struct{
 unsigned char  address;
 unsigned char  data;
 unsigned short write;
 unsigned char  deviceAddress;
}EEPROMParams;

typedef struct{
 unsigned char  data;
}EEPROMResults;

typedef struct{
 unsigned short firmwareVersion;
 unsigned short cameraID;
}GetVersionResults;

typedef struct{
 unsigned char  subCommand;
}ReadoutParams;

typedef struct{
 unsigned char  subCommand;
 unsigned char  ccd;
 unsigned short columns;
 unsigned short rows;
}ReadoutAreaParams;

typedef struct{
 unsigned char  subCommand;
 unsigned char  ccd;
 unsigned char  horzBin;
 unsigned char  vertBin;
 unsigned short left;
 unsigned short right;
}ReadoutSetupParams;

typedef struct{
 unsigned char  subCommand;
 unsigned char  ccd;
 unsigned char  vertBin;
 unsigned short rowWidth;
 unsigned short rows;
}ReadoutDumpParams;

typedef struct{
 unsigned char  subCommand;
 unsigned char  raiseIt;
}ReadoutVddParams;

typedef struct{
 unsigned char subCommand;
 unsigned char vddWasLow;
}ReadoutVddResults;

typedef struct {
 short rowWidth;       /* rows to download from camera    */
 short reqRowWidth;    /* pixels asked for by user        */
 short rowsPerPass;    /* rows to to per USB transactions */
 short rowsInFifo;     /* no. rows in fifo                */
 short bPipelineFull;  /* TRUE when already pipelined     */
}FifoInfo;

typedef struct{
 unsigned short imagingState;
 unsigned short trackingState;
 unsigned short shutterStatus;
 unsigned short ledStatus;
 unsigned short xPlusActive;
 unsigned short xMinusActive;
 unsigned short yPlusActive;
 unsigned short yMinusActive;
 unsigned short fanEnabled;
 unsigned short CFW6Active;
 unsigned short CFW6Input;
 unsigned char  shutterEdge;
 unsigned short filterState;
 unsigned short adSize;
 unsigned short filterType;
}StatusResults;

#define DO_VERSION_1		1  // initial version
#define DO_VERSION_2		2  // added something
#define DO_VERSION_3		3  // added doReportShutterError
#define DO_VERSION_4		4  // added doUseImagerDSM
#define DO_VERSION_5		5  // removed doUseImageDSM, doUSBGARev

#define DO_CURRENT_VERSION	DO_VERSION_5

typedef enum{
 DO_PER_DRIVER, 
 DO_TRUE, 
 DO_FALSE 
}DRIVER_OPTION_BOOL;

typedef struct{
 short doVersion;
 short doSize;
 short doBiasSubtraction;      // TRUE when bias subtraction enabled
 short doVddLowForIntegration; // TRUE when lower Vdd for integration
 short doAutoFreezeTE;	       // TRUE when auto freeze TE cooler for readout
 short doReportShutterErrors;  // TRUE when shutter errors reported 
 short doUSBFifoSize;	       // Size of FIFO in USB cameras
}UDRV_OPTIONS;

typedef struct {
 short         desiredPixels;
 unsigned char ccd;
}USBIGA;

#define USB_FIFO_SIZE           (4*4096) // size of Fifo in camera
#define USB_PIPELINE_SIZE        4       // size of AD pipeline register
//========================================================================
// Linux Data Structures
//========================================================================
//#define _USB_CHATTY_
//#define _SBIG_DEBUG_

typedef unsigned char  uchar;
//typedef unsigned short ushort;
//typedef unsigned long  ulong;
typedef float          IMG_DATA_TYPE;

#define SBIG_DRIVER_VERSION_BCD         0x0109
#define SBIG_LPT_DRIVER_VERSION_STR    "modsbiglpt.o Ver 1.12"
#define SBIG_USB_DRIVER_VERSION_STR    "modsbigusb.o Ver 1.12"

#define SBIG_LPT_DRV_DESCRIPTION       "SBIG LPT Linux Kernel Driver"
#define SBIG_USB_DRV_DESCRIPTION       "SBIG USB Linux Kernel Driver"
#define SBIG_DRV_AUTHOR                "Jan Soldan, jsoldan@asu.cas.cz"

#define SBIG_USB_DEVICE_NAME           "sbigusb"
#define SBIG_USB_VENDOR_ID		0x0D97
#define SBIG_USB_UNBOOTED_PRODUCT_ID	0x0001
#define SBIG_USB_BOOTED_PRODUCT_ID	0x0101
#define SBIG_USB_MINOR_BASE		192
#define SBIG_USB_MAX_DEVICES		8

#define MAX_DEVICE_NAME_SIZE            64  // device name string
#define BUFSIZE                         64  // i/o buffer
#define CURRENT_EEPROM_VERSION          1   // version of EEPROM contents
#define NEXT_ROM_MSN                    5        

#define STR_SBIG_LPT0                   "sbiglpt0"
#define STR_SBIG_LPT1                   "sbiglpt1"
#define STR_SBIG_LPT2                   "sbiglpt2"
#define STR_SBIG_USB0                   "sbigusb0"
#define STR_SBIG_USB1                   "sbigusb1"
#define STR_SBIG_USB2                   "sbigusb2"

#define	STR_IMAGING_CCD                 "Imaging"
#define	STR_TRACKING_CCD                "Tracking"

#define	STR_BIN1x1_ON                   "1x1"
#define	STR_BIN2x2_ON                   "2x2"
#define	STR_BIN3x3_ON                   "3x3"
#define	STR_BIN1x1_OFF                  "1x1 -"
#define	STR_BIN2x2_OFF                  "2x2 -"
#define	STR_BIN3x3_OFF                  "3x3 -"

#define	STR_FAN_ON                      "On"
#define	STR_FAN_OFF                     "Off"

#define	STR_ABG_OFF                     "Off"
#define	STR_ABG_LOW                     "Low"
#define	STR_ABG_MED                     "Med"
#define	STR_ABG_HIGH                    "High"
			      
#define	STR_WND_1                       "Wnd 1"
#define	STR_WND_2                       "Wnd 2"
#define	STR_WND_12                      "Wnd 1-2"

#define STR_SB_SELECT_FILE              "Select new SBIG or FITS file..."
#define STR_SB_DONE                     "Done."
#define STR_SB_READY                    "Ready."

#define STR_NONE_CALIBRATION            "None"
#define STR_BASIC_CALIBRATION           "LF = LF - DF"
#define STR_STANDARD_CALIBRATION        "LF = (LF - DF)/FF"
#define STR_MASTER_DF                   "DF = (DF1 + DF2 + ... + DFn) / n"

#define	DEFAULT_CAM1_CAPTION		"Camera 1"
#define	DEFAULT_CAM2_CAPTION		"Camera 2"
#define	DEFAULT_CAM3_CAPTION		"Camera 3"

#define	DEFAULT_FILE_NAME               "test.fits"
#define	DEFAULT_SELECT_DEVICE_MESSAGE   "Select Camera->Device, than press Link button."

#define	UNDEF_VALUE  -999

// Setpoint
#define	DEFAULT_SETPOINT_VALUE          -5.0
#define	DEFAULT_SETPOINT_TOLERANCE       0.5   
#define	DEFAULT_TEMP_STATUS_UPDATE       1500L

// Exposures
#define	DEFAULT_F1_EXPOSURE              5.0
#define	DEFAULT_F2_EXPOSURE              5.0
#define	DEFAULT_F3_EXPOSURE              5.0
#define	DEFAULT_F4_EXPOSURE              5.0
#define	DEFAULT_F5_EXPOSURE              5.0

// Filters' positions
#define	DEFAULT_F1_POSITION              1
#define	DEFAULT_F2_POSITION              2
#define	DEFAULT_F3_POSITION              3
#define	DEFAULT_F4_POSITION              4
#define	DEFAULT_F5_POSITION              5

// Filters' names
#define	DEFAULT_F1_NAME                  "None"
#define	DEFAULT_F2_NAME                  "Red"
#define	DEFAULT_F3_NAME                  "Green"
#define	DEFAULT_F4_NAME                  "Blue"
#define	DEFAULT_F5_NAME                  "H Alpha"


#define	DEFAULT_LIGHT_PATH		 "/home/soldan/wfoms/2003mmdd/raw/lf/"
#define	DEFAULT_DARK_PATH                "/home/soldan/wfoms/2003mmdd/raw/df/"
#define	DEFAULT_FILENAME_DIR             "/home/soldan/wfoms/"
#define	DEFAULT_FILENAME_FIL             "Images (*.fits *.FITS *.FTS *.fts *.st* *.ST*)"

// Exposure mode
#define	DEFAULT_SERIES                   1
#define	DEFAULT_LIGHT_SUB_SERIES         1
#define	DEFAULT_DARK_SUB_SERIES          1

// FOV 
#define	DEFAULT_DLE                      48.0
#define	DEFAULT_FLLE                     135.0
#define	DEFAULT_CCDWLE                   765
#define	DEFAULT_CCDHLE                   510
#define	DEFAULT_PIXELWLE                 9
#define	DEFAULT_PIXELHLE                 9

// Power
#define	DEFAULT_POWER_LIMIT              95

// Main window size
#define	DEFAULT_MAIN_WINDOW_WIDTH        10
#define	DEFAULT_MAIN_WINDOW_HEIGHT       10

#define BUF_LENGTH                       16384
#define	SBIG_FILE_HEADER_LENGTH          2048
#define	SBIG_STR_LENGTH                  64  
#define	MAX_PATH_LENGTH                  512  


// Under Linux we use loop around kernel function udelay(1000)
// which busy waits for one millisecond.

// define max USB timeout in jiffies for 1 second
#define	USB_MAX_TIMEOUT   100  
 
#define NO_CLEARS         4             // number of times to clear CCD at
				        // ticks before use full clear
// SBIG: FULL_CLEAR_TIME  6 ticks, ie. 330 ms 
#define FULL_CLEAR_TIME   330           // Linux: value in ms, ie. 330

// SBIG: NIBBLE_TIMEOUT = 5 ticks between chars received from camera,
// ie. 5 * 55 = 275 ms
#define NIBBLE_TIMEOUT    275	       // Linux: value in ms, ie. 275

// SBIG: IDLE_STATE_DELAY  3 ticks, ie. 3x55 = 165 ms
#define IDLE_STATE_DELAY  165          // Linux: value in ms, ie. 165

// SBIG: CONVERSION_DELAY = 500  queries (approx. 1us each) before
// A/D must signal not busy
#define CONVERSION_DELAY  500	       // Linux: value 500 queries

// SBIG: VCLOCK_DELAY 10 loops between vertical clocks on the imaging CCD
#define VCLOCK_DELAY      10	       // Linux: value 10 loops
#define ST1K_VCLOCK_X     10           // multiplier for VClock an ST1K

// number of pixels cleared in a block
#define CLEAR_BLOCK       10		

// after every N vertical clocks do a full horizontal clock on
#define DUMP_RATIO        5		

// number of pixels to average for bias
#define BIAS_WIDTH        128	

// CCD/Ambient thermistor constant
#define T0                    25.000
#define MAX_AD              4096.000
#define R_RATIO_CCD            2.570
#define R_BRIDGE_CCD          10.000
#define DT_CCD                25.000
#define R0                     3.000
#define R_RATIO_AMBIENT        7.791
#define R_BRIDGE_AMBIENT       3.000
#define DT_AMBIENT            45.000
#define SIDEREAL_DRIFT        15.04108444
#define PI                     3.1415926
#define	TEST_NUM_DRIFT_ROWS  100 
//------------------------------------------------------------------------
typedef enum{
	CAM_TEMP_EVENT = 65000,
	CAM_STATUS_TEXT_EVENT,
	CAM_TYPE_TEXT_EVENT,
        CAM_CLEAR_LABELS_EVENT,
	CAM_GRAB_DONE_EVENT,
	CAM_PB_SET_TOTAL_STEPS_EVENT,
	CAM_PB_SET_PROGRESS_EVENT,
	CAM_PB_RESET_EVENT
}CAMERA_EVENT;
//------------------------------------------------------------------------
typedef enum{
 INIT_EXPOSURE,
 EXPOSURE,
 INIT_READOUT_LINES,
 READOUT_LINES
}ACTION;
//------------------------------------------------------------------------
typedef struct{
 unsigned short action;
 unsigned long  start;
 unsigned long  position;
 unsigned long  end;
}LinuxCallbackParams;
//------------------------------------------------------------------------
typedef struct{
 char           pathName  [MAX_PATH_LENGTH];
 char           timeStamp [SBIG_STR_LENGTH];
 char           cameraName[SBIG_STR_LENGTH];
 char           filterName[SBIG_STR_LENGTH];
 int            fileSaveFlag;
 unsigned short ccd;          // CCD_REQUEST
 unsigned short abgState;     // ABG_STATE7
 unsigned short openShutter;  // SHUTTER_COMMAND
 double         exposureTime;
 unsigned short readoutMode;	    
 unsigned short top;
 unsigned short left;
 unsigned short width;
 unsigned short height;   
 unsigned short subtract;
}GrabImageParams;
//------------------------------------------------------------------------
typedef struct{
 unsigned short ccd;          // CCD_REQUEST
 unsigned short abgState;     // ABG_STATE7
 unsigned short openShutter;  // SHUTTER_COMMAND
 unsigned short readoutMode;
 float          exposureTime;
 unsigned short top;
 unsigned short left;
 unsigned short width;
 unsigned short height;   
 unsigned short subtract;
}GrabAreaParams;
//------------------------------------------------------------------------
typedef struct{
 // input 
 unsigned short driftImageWidth;         // lenght of drift scan image in pixels 
 unsigned short readoutMode;	         // 0=1x1, 1=2x2, 2=3x3 only
 double         focalLength;             // telescope's focal length [m]
 double         dec;                     // declination of center of the CCD  [degr]

 // output 
 int            ccdWidth;                // [pixels]
 int            ccdHeight;               // [pixels]

 double         pixelWidth;              // [um]
 double         pixelHeight;             // [um]
 double	        pixelFovWidth;           // [arcsec/pixel]
 double         pixelFovHeight;          // [arcsec/pixel]
 
 double         driftFovWidth;	         // [arcsec]
 double         driftFovHeight;          // [arcsec]
 double         driftRateDifference;     // [arcsec/sec]

 double         crossingExpTime;         // [s]
 double         driftExpTime;	         // [s]
 double         totalExpTime;            // [s]

 double         maxReadoutTime;          // [ms]
 double         rowReadoutTime;		 // [ms]
 double         success;                 // 0-100% 

}DriftScanIo;
//------------------------------------------------------------------------
typedef struct{
 unsigned short ccd;          // CCD_REQUEST
 unsigned short readoutMode;
 unsigned short height;
}DimmParams;
//------------------------------------------------------------------------
typedef struct{
 unsigned short ccd;          // CCD_REQUEST
 unsigned short readoutMode;
 unsigned short top;
 unsigned short left;
 unsigned short width;
 unsigned short height;   
}LinuxReadoutAreaParams;
//------------------------------------------------------------------------
typedef struct{
 unsigned char  reg;
 unsigned char  value;
}LinuxCameraOutParams;
//------------------------------------------------------------------------
typedef struct{
 unsigned char  *pBuffer;
 unsigned long  length;
}LinuxMicroblock;
//------------------------------------------------------------------------
typedef struct{
 IOC_GET_PIXELS_PARAMS  gpp;
 unsigned short *dest;
 unsigned long  length;
}LinuxGetPixelsParams;
//------------------------------------------------------------------------
typedef struct{
 IOC_GET_AREA_PARAMS  gap;
 unsigned short *dest;
 unsigned long  length;
}LinuxGetAreaParams;
//------------------------------------------------------------------------
typedef struct{
 unsigned short portBase;
 unsigned short portSpan;
}LptPortParams;
//------------------------------------------------------------------------
typedef struct{
 MY_LOGICAL	linkEstablished;
 unsigned short cameraType;
 char           deviceName[MAX_DEVICE_NAME_SIZE];
}LinuxGetLinkStatusResults;
//------------------------------------------------------------------------
#define CB_CCD_TYPE_FRAME_TRANSFER   0x0001 //b0=1 is frame transfer CCD
#define CONFIG_LOW_TRACKER_IS_237         1
//========================================================================

#endif // _SBIGDEF_H_
