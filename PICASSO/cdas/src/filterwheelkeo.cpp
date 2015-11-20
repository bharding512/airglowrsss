/***************************************************************************
                          filterWheelFLI.cpp  -  description
                             -------------------
    begin                : Fri Oct 21 2005
    copyright            : (C) 2004,5 by Jonathan Makela and Ethan Miller
    email                : jmakela@uiuc.edu, esmiller@uiuc.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "filterwheelkeo.h"
#include "system.h"
#include <sstream>
#include <strstream>
#include "libfli.h"

#ifndef INVALID_HANDLE_VALUE
  #define INVALID_HANDLE_VALUE -1
#endif

extern System mySystem;

FilterWheelKEO::FilterWheelKEO(){
	myFilterOK = true;
	myShutterOpen = false;
	myNumFilters = 5;
	myShutterPos = 0;	// We don't techinically use a shutter position here
				// but we need this variable to be set
	
	myFWDevice = 0;
	myTempControlDevice = 1;

	// Initialize COM ports
	myFilterOK = SetUp_FilterWheel();

	// Always start with a closed shutter
	CloseShutter();
}

FilterWheelKEO::~FilterWheelKEO(){
	// Clean up
	myFilterOK = false;
	CloseShutter();

	// Close the ports
	::close(myFWControlFD);
	::close(myTempControlFD);
}

/** Read filterwheel class from file */
void FilterWheelKEO::Read(ifstream & is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#FilterWheel" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #FilterWheel found");
    return;
  }

  // Read the data
  is >> mySystem.myConfiguration.myFilterWheelType;
  is >> myNumFilters;
  is >> myShutterPos;
  is >> myFWDevice;
  is >> myTempControlDevice;

  // Close the shutter now that we've read in the position
  CloseShutter();
}

/** Write the filter class to config file */
void FilterWheelKEO::Write(ofstream & os){
  // Put the header label
  os << "#FilterWheel" << endl;

  // Write the CCD settemp
  os << mySystem.myConfiguration.myFilterWheelType << endl;
  os << myNumFilters << endl;
  os << myShutterPos << endl;
  os << myFWDevice << endl;
  os << myTempControlDevice << endl;
}

/** No descriptions */
void FilterWheelKEO::WriteXML(ofstream & os){
}

/** Closes the filterwheel shutter. */
void FilterWheelKEO::CloseShutter(){
	int i;
	int len;

	char *szOut = "d=0\r";
	DWORD out;

	if(myFWControlFD > 0)
	{
		// Create string to tell filterwheel where to go
		len = strlen(szOut);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut, myFWControlFD))
			{
				// There was an error
				string str("Cannot close shutter: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return;
			}
			*szOut++;
		}

		// Create string to tell filterwheel to move
	        char *szOut2 = "GOSUB1\r";
		len = strlen(szOut2);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut2, myFWControlFD))
			{
				// There was an error
				string str("Cannot change close shutter: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return;
			}
			szOut2++;
		}
		
	}

	// Set the variable containing the shutter state
	myShutterOpen = false;

	return;
}

/** Returns the current position of the filterwheel. */
int FilterWheelKEO::GetCurrentPosition() {
	int i;
	int myPos;
	int len;

	char *szOut="g=-1\r";
	DWORD out;
	char myBuffer[4096];

	if(myFWControlFD > 0)
	{
		// Create string to tell filterwheel where to go
		len = strlen(szOut);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut, myFWControlFD))
			{
				// There was an error
				string str("Cannot querry filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return -999;
			}
			*szOut++;
		}

		// Create string to tell filterwheel to move
                char *szOut2 = "GOSUB4\r";
		//sprintf(szOut, "GOSUB4\r");
		len = strlen(szOut2);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut2, myFWControlFD))
			{
				// There was an error
				string str("Cannot querry filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return -999;
			}
			szOut2++;
		}
		
		// Read the string on the line
		int bytesRead =::read(myFWControlFD,myBuffer,4096);
		if(bytesRead == 0)
		{
			// There was an error
			string str("Cannot querry filter position: ");
     			mySystem.myConfiguration.myLogger.WriteErrorString(str);
      			myFilterOK = false;
      			return -999;
		}
	
		// Parse the return string
		sscanf(myBuffer,"FILT:%d\r", &myPos);
	}

	return myPos;
}

/** Returns whether the shutter is open or closed. */
bool FilterWheelKEO::IsShutterOpen(){
  return myShutterOpen;
}

/** Moves the filterwheel to the position given by pos. */
int FilterWheelKEO::MoveToPosition(int pos){
	int i;
	int len;

	char szOut[4];// = "g=1\r";
	DWORD out;

	if(myFWControlFD > 0)
	{
		// Create string to tell filterwheel where to go
		char *szOut1 = "g=3\r";
		sprintf(szOut, "g=%d\r", pos);
		len = strlen(szOut);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(szOut[i], myFWControlFD))
			{
				// There was an error
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return -999;
			}
//			*szOut++;
		}

		sleep(1);
		
		// Create string to tell filterwheel to move
		char *szOut2 = "GOSUB4\r";
		//sprintf(szOut, "GOSUB4\r");
		len = strlen(szOut2);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut2, myFWControlFD))
			{
				// There was an error
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return -999;
			}
			szOut2++;
		}
		
	}

	// Sleep for a second to let the filterwheel do it's thing
 	sleep(1);

	// We were successful
	return pos;
}

/** Initializes the filterwheel. */
bool FilterWheelKEO::SetUp_FilterWheel(){
	// Assume the FW is OK
	myFilterOK = true;
	ostrstream dev;
	dev << "/dev/ttyS" << myFWDevice << ends;

	// Set up the COM ports (the settings will be the same for both)
	int flags = 0;

	///////////////////////////////////
	// Open the FW device as read/write
	flags = O_RDWR;
//	myFWControlFD = open(myFWDevice, flags | O_NDELAY);
	myFWControlFD = open(dev.str(), flags | O_NDELAY);


	// Check to make sure we were able to open the files
	if(myFWControlFD < 0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice;
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		return false;
	}

	// Flush the port
	tcflush(myFWControlFD, TCIOFLUSH);
	int n = fcntl(myFWControlFD, F_GETFL, 0);
	fcntl(myFWControlFD, F_SETFL, n & ~O_NDELAY);

	// The Terminal I/O structure and required variables
	speed_t _baud;
	bzero(&newtio, sizeof(newtio));
	
	// Get the current attributes
	if(tcgetattr(myFWControlFD, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice << ", Cannot get attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		return false;
	}

	// Set the speed of the serial port
	_baud = B9600;
	cfsetospeed(&newtio, (speed_t)_baud);
	cfsetispeed(&newtio, (speed_t)_baud);

	// Set 8 data bits, N parity, no handshaking, 1 stop bit
	newtio.c_cflag = (newtio.c_cflag & ~CSIZE) | CS8;
	newtio.c_cflag |= CLOCAL | CREAD;
	newtio.c_cflag &= ~(PARENB | PARODD);
	newtio.c_cflag &= ~CRTSCTS;
	newtio.c_cflag &= ~CSTOPB;

	// Other flags and turnoff software handshake
	newtio.c_iflag = IGNBRK;
	newtio.c_iflag &= ~(IXON|IXOFF|IXANY);
	newtio.c_lflag = 0;
	newtio.c_oflag = 0;
	newtio.c_cc[VTIME] = 1;
	newtio.c_cc[VMIN] = 60;

	if(tcsetattr(myFWControlFD, TCSANOW, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice << ", Cannot set attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		return false;
	}

	int mcs = 0;
	ioctl(myFWControlFD, TIOCMGET, &mcs);
	mcs |= TIOCM_RTS;
	ioctl(myFWControlFD, TIOCMSET, &mcs);

	if(tcgetattr(myFWControlFD, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice << ", Cannot verify attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		return false;
	}

	newtio.c_cflag &= ~CRTSCTS;
	if(tcsetattr(myFWControlFD, TCSANOW, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice << ", Cannot set CRTSCTS";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		return false;
	}

/*	// Create string to tell filterwheel to initialize
        char *szOut2 = "GOSUB5\r";
        //sprintf(szOut, "GOSUB1\r");
        int len = strlen(szOut2);

        // Send this to the filterwheel port
        for(int i = 0; i < len; i++)
        {
                if(!SendByte(*szOut2, myFWControlFD))
                {
                        // There was an error
                        string str("Cannot initialize KEO filterwheel: ");
                        mySystem.myConfiguration.myLogger.WriteErrorString(str);
                        myFilterOK = false;
                        return false;
                }
                szOut2++;
        }
*/

	///////////////////////////////////
	// Open the devices as read/write
	flags = 0;
	flags = O_RDWR;
	dev << "/dev/ttyS" << myTempControlDevice << ends;
//	myTempControlFD = open(myTempControlDevice, flags | O_NDELAY);
//	myTempControlFD = open("/dev/ttyS1", flags | O_NDELAY);
 	myTempControlFD = open(dev.str(), flags | O_NDELAY);

	// Check to make sure we were able to open the files
	if(myTempControlFD < 0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open temperature port /dev/ttyS" << myTempControlDevice;
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		return false;
	}

	// Flush the port
	tcflush(myTempControlFD, TCIOFLUSH);
	n = fcntl(myTempControlFD, F_GETFL, 0);
	fcntl(myTempControlFD, F_SETFL, n & ~O_NDELAY);
	
	// Get the current attributes
	if(tcgetattr(myTempControlFD, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open temperature port /dev/ttyS" << myTempControlDevice << ", Cannot get attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		::close(myTempControlFD);
		return false;
	}

	// Set the speed of the serial port
	_baud = B9600;
	cfsetospeed(&newtio, (speed_t)_baud);
	cfsetispeed(&newtio, (speed_t)_baud);

	// MAY NEED TO BE 7 BIT, ODD PARITY
	
	// Set 8 data bits, N parity, no handshaking, 1 stop bit
	newtio.c_cflag = (newtio.c_cflag & ~CSIZE) | CS8;
	newtio.c_cflag |= CLOCAL | CREAD;
	newtio.c_cflag &= ~(PARENB | PARODD);
	newtio.c_cflag &= ~CRTSCTS;
	newtio.c_cflag &= ~CSTOPB;

	// Other flags and turnoff software handshake
	newtio.c_iflag = IGNBRK;
	newtio.c_iflag &= ~(IXON|IXOFF|IXANY);
	newtio.c_lflag = 0;
	newtio.c_oflag = 0;
	newtio.c_cc[VTIME] = 1;
	newtio.c_cc[VMIN] = 60;

	if(tcsetattr(myTempControlFD, TCSANOW, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open temperature port /dev/ttyS" << myTempControlDevice << ", Cannot set attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		::close(myTempControlFD);
		return false;
	}

	mcs = 0;
	ioctl(myTempControlFD, TIOCMGET, &mcs);
	mcs |= TIOCM_RTS;
	ioctl(myTempControlFD, TIOCMSET, &mcs);

	if(tcgetattr(myTempControlFD, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open temperature port /dev/ttyS" << myTempControlDevice << ", Cannot verify attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		::close(myTempControlFD);
		return false;
	}

	newtio.c_cflag &= ~CRTSCTS;
	if(tcsetattr(myTempControlFD, TCSANOW, &newtio)!=0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open temperature port /dev/ttyS" << myTempControlDevice << ", Cannot set CRTSCTS";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		::close(myFWControlFD);
		::close(myTempControlFD);
		return false;
	}

	return true;
}

/** Sets the filter position containing the shutter. */
void FilterWheelKEO::SetShutterPosition(int num){
  myShutterPos = num;
}
/** Returns the shutter position. */
int FilterWheelKEO::GetShutterPosition(){
  return myShutterPos;
}


/*!
    \fn FilterWheelFLI::OpenShutter()
 */
void FilterWheelKEO::OpenShutter()
{
	int i;
	int len;

	char *szOut="d=1\r";
	DWORD out;

	if(myFWControlFD > 0)
	{
		// Create string to tell filterwheel where to go
	//	sprintf(szOut, "d=1\r");
		len = strlen(szOut);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut, myFWControlFD))
			{
				// There was an error
				string str("Cannot open shutter: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return;
			}
			*szOut++;
		}

		// Create string to tell filterwheel to move
		char *szOut2 = "GOSUB1\r";
		//sprintf(szOut, "GOSUB1\r");
		len = strlen(szOut2);

		// Send this to the filterwheel port
		for(i = 0; i < len; i++)
		{
			if(!SendByte(*szOut2, myFWControlFD))
			{
				// There was an error
				string str("Cannot open shutter: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
      				return;
			}
			szOut2++;
		}
		
	}

	// Set the variable containing the shutter state
	myShutterOpen = true;

  return;
}

bool FilterWheelKEO::SendByte(char c, int myFD)
{
	// Check to make sure we have a valid file descriptor
	if(myFD==-1)
	{
		mySystem.myConfiguration.myLogger.WriteErrorString("Trying to write to unopened serial port");
		return false;
	}
	
	// Write the byte to the port
	int res=::write(myFD, &c, 1);
	
	// Check to make sure the operation worked
	if(res < 1)
	{
		mySystem.myConfiguration.myLogger.WriteErrorString("SendByte failed in FilterWheelKEO");
		return false;
	}

	// Sleep for a millisecond
	millisleep(1);

	return true; 
}

/*!
\fn FilterWheelKEO::millisleep(int ms)
*/
void FilterWheelKEO::millisleep(int ms)
{
	// Sleep for the requested number of milliseconds
	if(ms > 0)
	{
		struct timeval tv;
		tv.tv_sec=0;
		tv.tv_usec=ms*1000;
		select(0,0,0,0,&tv);
	}
}

void FilterWheelKEO::SetNumFilters(int num){
  myNumFilters = num;
}

int FilterWheelKEO::GetNumFilters(){
  return myNumFilters;
}
