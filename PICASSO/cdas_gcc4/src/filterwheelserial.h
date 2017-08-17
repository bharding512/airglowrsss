#ifndef FILTERWheelSERIAL_H
#define FILTERWheelSERIAL_H

#ifndef _LIBFLI_H_
  #include "libfli.h"
#endif

#include <filterwheel.h>
#include <termios.h>
#include <iostream>
using namespace std;
#include <fcntl.h>
#include <sys/ioctl.h>
//#include <linux/ioctl.h>
//#include <i386-linux-gnu/sys/ioctl.h>

#define FLISTRINGLENGTH 255

class FilterWheelSERIAL : public FilterWheel  {
public:
	FilterWheelSERIAL();
	~FilterWheelSERIAL();
  /** No descriptions */
  void WriteXML(ofstream & os);
  /** Write the filter class to config file */
  void Write(ofstream & os);
  /** Read filterwheel class from file */
  void Read(ifstream & is);
  /** Initializes the filterwheel. */
  bool SetUp_FilterWheel();
  /** Moves the filterwheel to the position given by pos. */
  int MoveToPosition(int pos);
  /** Returns whether the shutter is open or closed. */
  bool IsShutterOpen();
  /** Returns the current position of the filterwheel. */
  int GetCurrentPosition();
  /** Closes the filterwheel shutter. */
  void CloseShutter();
  /** Return a string describing the last error encountered */
  string GetErrorString();
  /** Sets the filter position containing the shutter. */
  void SetShutterPosition(int num);
  /** Returns the shutter position. */
  int GetShutterPosition();
    void OpenShutter();
  void SetNumFilters(int num);
  int GetNumFilters();
public: // Public attributes
  /**  */
  bool myFilterOK;
  /** The number of filters in the current filterwheel. */
  int myNumFilters;
  /** The location of the shutter position. */
    int myShutterPos;
    int myTempControlDevice;
    int myFWDevice;
private:
  char usbdevice[FLISTRINGLENGTH];
  flidev_t filterwheelDevice;
  int myFWControlFD;
  int myTempControlFD;
  bool myShutterOpen;
  bool SendByte(char c, int myFD);
  void millisleep(int ms);
  int ishome();
  int move_FWSerial(int k);
  int myCurrentPosition;
  int myPositionNumber;
 
struct termios newtio;
};

#endif

