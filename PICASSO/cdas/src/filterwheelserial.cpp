#include "filterwheelserial.h"
#include "system.h"
#include <sstream>
#include <strstream>
#include "libfli.h"

#ifndef INVALID_HANDLE_VALUE
  #define INVALID_HANDLE_VALUE -1
#endif

extern System mySystem;

FilterWheelSERIAL::FilterWheelSERIAL(){
	myFilterOK = true;
	myShutterOpen = false;
	myNumFilters = 5;
	myShutterPos = 0;	// We don't techinically use a shutter position here
				// but we need this variable to be set
	myFWDevice = 2;
	myTempControlDevice = 1;

}

FilterWheelSERIAL::~FilterWheelSERIAL(){
	// Clean up
	myFilterOK = false;
	CloseShutter();
	// Close the ports
	::close(myFWControlFD);
	::close(myTempControlFD);
}

/** Read filterwheel class from file */
void FilterWheelSERIAL::Read(ifstream & is){
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
void FilterWheelSERIAL::Write(ofstream & os){
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
void FilterWheelSERIAL::WriteXML(ofstream & os){
}

/** Closes the filterwheel shutter. */
void FilterWheelSERIAL::CloseShutter(){

	// Set the variable containing the shutter state
	myShutterOpen = false;

	return;
}

/** Returns the current position of the filterwheel. */
int FilterWheelSERIAL::GetCurrentPosition() {
		
	  sleep(3);
	  return myCurrentPosition;
	//}
}

/** Returns whether the shutter is open or closed. */
bool FilterWheelSERIAL::IsShutterOpen(){
  return myShutterOpen;
}

/** Moves the filterwheel to the position given by pos.
    Here pos is  */
int FilterWheelSERIAL::MoveToPosition(int pos){
	//Gohome then move to pos x 4000
	int h_pos=myPositionNumber;

//	move_FWSerial(h_pos);
//	sleep(5);

//	while(ishome()<=0)
//	{
//	move_FWSerial(h_pos);
//	h_pos+=100;
//	}	
	h_pos=h_pos+(pos-1)*4000;
	move_FWSerial(h_pos);
	// We were successful

	myCurrentPosition = pos;
	return myCurrentPosition;
}

/** Initializes the filterwheel,by moving it to the home position.
    Offset is the position offset from where the HallSensor is high so as to center the Filter. */ 
bool FilterWheelSERIAL::SetUp_FilterWheel(){
	// Initialize FW motion to OFF
	
	int offset = 200;
	
	int pos=0;

	move_FWSerial(pos);
	sleep(5);

	while(ishome()<=0)
	{
	move_FWSerial(pos);
	pos+=100;
	}
	
	move_FWSerial(pos-offset);
	myPositionNumber=pos-offset;
	return true;
}

/** Sets the filter position containing the shutter. */
void FilterWheelSERIAL::SetShutterPosition(int num){
  myShutterPos = num;
}
/** Returns the shutter position. */
int FilterWheelSERIAL::GetShutterPosition(){
  return myShutterPos;
}


/*!
    \fn FilterWheelFLI::OpenShutter()
 */
void FilterWheelSERIAL::OpenShutter()
{
	// Set the variable containing the shutter state
	myShutterOpen = true;

  return;
}

bool FilterWheelSERIAL::SendByte(char c, int myFD)
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


void FilterWheelSERIAL::millisleep(int ms)
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

int FilterWheelSERIAL::ishome()
{
	myFWDevice=0;
	myFWControlFD;
	myFilterOK=true;	
	//FW_ON=true;
	ostrstream dev;

	dev << "/dev/ttyUSB" << myFWDevice << ends;// Set up the COM ports (the settings will be the same for both)

    int flags = 0;

    ///////////////////////////////////

    // Open the FW device as read/write

    flags = O_RDWR; //unsure

// myFWControlFD = open(myFWDevice, flags | O_NDELAY);

    myFWControlFD = open(dev.str(), flags | O_NDELAY);

    // Check to make sure we were able to open the files

    if(myFWControlFD < 0)
    {

          // An error occured
	  ostringstream sout;
    	  sout << "Error 'fopen' Cannot open FW port /dev/ttyUSB" << myFWDevice<<endl;
	  mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
	  //FW_ON=false;
	  myFilterOK=false;	
          return -1;

    }

    // Flush the port

    tcflush(myFWControlFD, TCIOFLUSH);

    int n = fcntl(myFWControlFD, F_GETFL, 0);

    fcntl(myFWControlFD, F_SETFL, n & ~O_NDELAY);

    // The Terminal I/O structure and required variables

    speed_t _baud;

    //tcgetattr(myFWControlFD,&newtio);//maybe wrong

    bzero(&newtio, sizeof(newtio));//

    // Get the current attributes

    if(tcgetattr(myFWControlFD, &newtio)!=0)
    {

          // An error occured
	  ostringstream sout;
    	  sout << "Error Current Att. Cannot open FW port /dev/ttyUSB" << myFWDevice << ", Cannot get attributes";
	  mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
	  myFilterOK=false;
	  //FW_ON=false;
          ::close(myFWControlFD);
          return -1;

    }

    // Set the speed of the serial port

    _baud = B9600;

    cfsetospeed(&newtio, (speed_t)_baud);

    cfsetispeed(&newtio, (speed_t)_baud);

    // Set 8 data bits, N parity, no handshaking, 1 stop bit

    newtio.c_cflag = (newtio.c_cflag & ~CSIZE) | CS8;

    newtio.c_cflag |= CLOCAL | CREAD;

    newtio.c_cflag &= ~(PARENB | PARODD);

    newtio.c_cflag |= CRTSCTS;

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
                sout << "Error set Att. Cannot open FW port /dev/ttyUSB" << myFWDevice << ", Cannot set attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		//FW_ON=false;
		::close(myFWControlFD);
		return -1;

    }

if(tcgetattr(myFWControlFD, &newtio)!=0)
    {

          // An error occured
	  ostringstream sout;        
sout << "Error get Att.2 Cannot open FW port /dev/ttyUSB" << myFWDevice << ", Cannot verify attributes"<<endl;
	  mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
	  myFilterOK=false;
	  //FW_ON=false;
          ::close(myFWControlFD);
          return -1;

    }

    newtio.c_cflag |= CRTSCTS;

    if(tcsetattr(myFWControlFD, TCSANOW, &newtio)!=0)
    {

          // An error occured
		ostringstream sout;
		sout << "Error Set Att 2. Cannot open FW port /dev/ttyUSB" << myFWDevice << ", Cannot set CRTSCTS"<<endl;
          	mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK = false;
		//FW_ON=false;
		::close(myFWControlFD);
		return -1;         
    }

int mcs = 0;

if(ioctl(myFWControlFD, TIOCMGET, &mcs)==-1)
{
cout<<"TIOCMGET failed"<<endl;
//FW_ON=false;
return -1;
}
int status = mcs & TIOCM_CTS;

//FW_ON=false;
::close(myFWControlFD);

return status;
}


int FilterWheelSERIAL::move_FWSerial(int k)
{
	myFWDevice=0;
	myFWControlFD;
	//FW_ON=true;

	ostrstream dev;

	dev << "/dev/ttyS" << myFWDevice << ends;

	int flags=0;

	//Open the FW as read/write

	flags=O_RDWR;

	myFWControlFD = open(dev.str(), flags | O_NDELAY);

	//Check to make sure we are able to open the files

	if(myFWControlFD < 0)
	{
		// An error occured
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice;
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK=false;
		//FW_ON=false;
		return -1;
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
		myFilterOK=false;
		//FW_ON=false;
		::close(myFWControlFD);
		return -1;
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
		myFilterOK=false;
		//FW_ON=false;
		::close(myFWControlFD);
		return -1;
	}

	if(tcgetattr(myFWControlFD, &newtio)!=0)
	{
		// An error occured
	
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice << ", Cannot verify attributes";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK=false;
		//FW_ON=false;
		::close(myFWControlFD);
		return -1;
	}

	newtio.c_cflag &= ~CRTSCTS;

	if(tcsetattr(myFWControlFD, TCSANOW, &newtio)!=0)
	{
		// An error occured
	
		ostringstream sout;
		sout << "Cannot open FW port /dev/ttyS" << myFWDevice << ", Cannot set CRTSCTS";
		mySystem.myConfiguration.myLogger.WriteErrorString(sout.str());
		myFilterOK=false;
		//FW_ON=false;
		::close(myFWControlFD);
		return -1;
	}

	int i,len;
	char ch[6];
	char *szOut,*szOut2,*szOut3,*szOut4,*szOut5;

		if(myFWControlFD>0)
		{
			szOut="MP\r";
			len=strlen(szOut);
			for(i=0;i<len;i++)
			{
				if(!SendByte(*szOut, myFWControlFD))
				{
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
				//FW_ON=false;
					return -1;
				}
			szOut++;
			}

			szOut2="A=60\r";
			len=strlen(szOut2);
			for(i=0;i<len;i++)
			{
				if(!SendByte(*szOut2, myFWControlFD))
				{
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
				//FW_ON=false;
				return -1;
				}
			szOut2++;
			}

			szOut3="V=150000\r";
			len=strlen(szOut3);
			for(i=0;i<len;i++)
				{
					if(!SendByte(*szOut3, myFWControlFD))
					{
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
				//FW_ON=false;
					return -1;
					}
				szOut3++;
				}

			szOut4;
			sprintf(ch,"%d",k);
			szOut4 = (char *)malloc(5 + strlen(ch)*sizeof(char));
			//szOut4= new char[5+strlen(ch)];
			strcpy(szOut4, "P=");
			strcat(szOut4, ch);
			strcat(szOut4, "\r");

			len=strlen(szOut4);
			for(i=0;i<len;i++)
			{
				if(!SendByte(*szOut4, myFWControlFD))
				{
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
				//FW_ON=false;
				return -1;
				}
				szOut4++;
			}

			szOut5="G\r";
			len=strlen(szOut5);
			for(i=0;i<len;i++)
			{
				if(!SendByte(*szOut5, myFWControlFD))
				{
				string str("Cannot change filter position: ");
     				mySystem.myConfiguration.myLogger.WriteErrorString(str);
      				myFilterOK = false;
				//FW_ON=false;
				return -1;
				}
				szOut5++;
			}
		}
//FW_ON=false;
//printf("Done moving fw\n");
::close(myFWControlFD);
return 1;
}

void FilterWheelSERIAL::SetNumFilters(int num){
  myNumFilters = num;
}

int FilterWheelSERIAL::GetNumFilters(){
  return myNumFilters;
}
