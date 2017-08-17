/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Wed Aug  6 12:05:11 EDT 2003
    copyright            : (C) 2003 by Jonathan Makela
    email                : jmakela@nrl.navy.mil
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include "system.h"
#include <pthread.h>
#include "ServerSocket.h"
#include "SocketException.h"
#include "ServerCommands.h"
#include <time.h>

System mySystem;

using namespace std;

// The thread function definitions
void *server_function(void *threadid);
void *control_function(void *threadid);

// The command parser
void parse_command(string command, string *result);

// The mutex definition
pthread_mutex_t mutex;

int main(int argc, char *argv[])
{
  // String var
  string myString;

  // Return code variable
  int rc;

  // The threads;
  pthread_t server;
  pthread_t control;

  // Load the configuration file
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Starting CDAS v2.0");
  myString = "Reading in configuration file: " + string(mySystem.myConfFileName);
  mySystem.myConfiguration.myLogger.WriteSystemLogString(myString);
  mySystem.myConfiguration.LoadConfiguration(mySystem.myConfFileName);
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Read in configuration file");
  
  // Initialize the system
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Initializing System");

  // We are now running
  mySystem.myIsRunning = true;

  // Save the XML configuration file
  myString = "Writing XML file: " + mySystem.myConfiguration.myDirectories.GetQuickLookPath() + mySystem.myXMLFilename;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(myString);
  mySystem.myConfiguration.SaveXML(mySystem.myConfiguration.myDirectories.GetQuickLookPath() + mySystem.myXMLFilename);

  // Initialize the Mutex that helps the thread controls
  pthread_mutex_init(&mutex, NULL);

  // Initialize the two threads
  int theadid = 0;
  rc = pthread_create(&server, NULL, server_function, (void *) theadid);
  theadid = 1;
  rc = pthread_create(&control, NULL, control_function, (void *) theadid);

  // An infinite loop while we are running
  while(mySystem.myIsRunning) {
    sleep(1);
  }

  // The system is stopped, wait for the server and control threads to end
  while(mySystem.myServerRunning) {
    sleep(1);
  }
  while(mySystem.myControlRunning) {
    sleep(1);
  }

  // Save the XML configuration file
  myString = "Writing XML file: " + mySystem.myConfiguration.myDirectories.GetQuickLookPath() + mySystem.myXMLFilename;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(myString);
  mySystem.myConfiguration.SaveXML(mySystem.myConfiguration.myDirectories.GetQuickLookPath() + mySystem.myXMLFilename);
  
  // Save the configuration file
  mySystem.myConfiguration.SaveConfiguration(mySystem.myConfFileName);
  
  return EXIT_SUCCESS;
}

void *server_function(void *threadid)
{
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Starting server thread");
  mySystem.myServerRunning = true;

  // Create a System object to store information in before the update command
  // is received.  It's initial state is the global server.

  try
  {
    // Create the socket
    ServerSocket server(30000);

    while(mySystem.myIsRunning)
    {
      ServerSocket new_sock;
      server.accept(new_sock);

      try
      {
        while(true)
        {
          string data;
          string output;
          ostringstream sout;

          new_sock >> data;

          // We've received data, try and parse the command
          parse_command(data, &output);

          new_sock << output;
        }
      }
      catch(SocketException&)
      {
      }
    }
  }
  catch(SocketException &e)
  {
    cout << "Exception was caught: " << e.description() << "\nExiting\n";
  }

  mySystem.myConfiguration.myLogger.WriteSystemLogString("Stopping server thread");
  mySystem.myServerRunning = false;
  
  pthread_exit(NULL);
}

void parse_command(string command, string *result)
{
  // Parses a command of the form "action: value"
  int action;
  string str;
  string filename;
  ostringstream sout;
  char *temp;
  int tempBool;

  action = atoi(strtok((char *)command.c_str(), ":"));

  // Use a switch statement to go through all of the different possible
  // commands.  The commands are defined in ServerCommands.h
  switch(action)
  {
    case SERVER_STOP:
    {
      pthread_mutex_lock(&mutex);
      *result = "Stop command received.  Stopping program...";
      mySystem.myIsRunning = false;
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SCHEDULE_MODE:
    {
      int mode;
      temp = strtok(NULL, ":");
      mode = atoi(temp);
      pthread_mutex_lock(&mutex);
      switch(mode)
      {
        case SCHEDULE_OFF:
	{
          str = "OFF";
          break;
	}
        case SCHEDULE_AUTO:
	{
          str = "AUTO";
          break;
	}
        case SCHEDULE_MANUAL:
	{
          time_t start, stop;
          char time_string[64];
          str = "MANUAL";
          temp = strtok(NULL, ":");
          if(temp == NULL)
          {
            // There's an error in the string sent.
            // Keep the current scheduling mode.
            str = "ERROR: Not enough arguments to set MANUAL";
            mode = mySystem.myConfiguration.mySchedule.GetScheduleMode();
            break;
          }
          start = atol(temp);
          temp = strtok(NULL, ":");
          if(temp == NULL)
          {
            // There's an error in the string sent.
            // Keep the current scheduling mode.
            str = "ERROR: Not enough arguments to set MANUAL";
            mode = mySystem.myConfiguration.mySchedule.GetScheduleMode();
            break;
          }
          stop = atol(temp);
          strftime(time_string, 64, "%a %b %d %H:%M:%S %Y", localtime(&start));
          str = str + "(" + time_string + "-";
          strftime(time_string, 64, "%a %b %d %H:%M:%S %Y", localtime(&stop));
          str = str + time_string + ")";
          mySystem.myConfiguration.mySchedule.SetStartTime(*localtime(&start));
          mySystem.myConfiguration.mySchedule.SetStopTime(*localtime(&stop));
          break;
	}
        default:
	{
          str = "???";
          break;
	}
      }
      *result = "Switching schedule mode to " + str;
      mySystem.myConfiguration.mySchedule.SetScheduleMode(mode);
      pthread_mutex_unlock(&mutex);
      break;
    }
    case ANGLES_SET:
    {
      temp = strtok(NULL, ":");
      int cmd;
      double angle;
      bool set;
      cmd = atoi(temp);
      pthread_mutex_lock(&mutex);
      switch(cmd)
      {
        case ANGLE_MOON:
	{
          // Change the moon angle
          temp = strtok(NULL, ":");
          if(temp == NULL)
          {
            // There's an error in the string sent.
            // Keep the current scheduling mode.
            *result = "ERROR: Not enough arguments to set ANGLES_SET";
            break;
          }
          angle = atof(temp);
          sout << angle;
          *result = "Change moon angle to " + sout.str();
          mySystem.myConfiguration.mySchedule.SetMoonAngle(angle);
          break;
	}
        case ANGLE_SUN:
	{
          // Change the sun angle
          temp = strtok(NULL, ":");
          if(temp == NULL)
          {
            // There's an error in the string sent.
            // Keep the current scheduling mode.
            *result = "ERROR: Not enough arguments to set ANGLES_SET";
            break;
          }
          angle = atof(temp);
          sout << angle;
          *result = "Change sun angle to " + sout.str();
          mySystem.myConfiguration.mySchedule.SetSunAngle(angle);
          break;
	}
        case ANGLE_DOMOON:
	{
          // Change the myMoonSet flag
          temp = strtok(NULL, ":");
          if(temp == NULL)
          {
            // There's an error in the string sent.
            // Keep the current scheduling mode.
            *result = "ERROR: Not enough arguments to set ANGLES_SET";
            break;
          }
          set = (bool) atoi(temp);
          sout << set;
          *result = "Change moon set to " + sout.str();
          mySystem.myConfiguration.mySchedule.SetMoonSet(set);
          break;
	}
        default:
	{
          *result = "ANGLES_SET command not recognized";
	}
      }
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_SITENAME:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SITENAME";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySite.SetSiteName(sout.str());
      *result = "Change sitename to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_SITEABBR:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SITEABBR";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySite.SetAbbreviation(temp);
      *result = "Change site abbreviation to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_LATITUDE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_LATITUDE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySite.SetLatitude(atof(temp));
      *result = "Change site latitude to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_SCHEDULELATITUDE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SCHEDULELATITUDE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySchedule.SetLatitude(atof(temp));
      *result = "Change schedule latitude to " + sout.str();
      pthread_mutex_unlock(&mutex);

      // Autoset times if needed
      if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
        mySystem.myConfiguration.mySchedule.AutoSetTimes();

      break;
    }
    case SET_LONGITUDE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_LONGITUDE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySite.SetLongitude(atof(temp));
      *result = "Change site longitude to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_SCHEDULELONGITUDE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SCHEDULELONGITUDE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySchedule.SetLongitude(atof(temp));
      *result = "Change schedule longitude to " + sout.str();
      pthread_mutex_unlock(&mutex);

      // Autoset times if needed
      if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
        mySystem.myConfiguration.mySchedule.AutoSetTimes();

      break;
    }
    case SET_ALTITUDE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_ALTITUDE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.mySite.SetHeight(atof(temp));
      *result = "Change site altitude to " + sout.str();
      pthread_mutex_unlock(&mutex);

      // Autoset times if needed
      if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
        mySystem.myConfiguration.mySchedule.AutoSetTimes();
        
      break;
    }
    case SET_XBIN:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_XBIN";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetBinX(atoi(temp));
      *result = "Change X-bin to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_YBIN:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_YBIN";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetBinY(atoi(temp));
      *result = "Change Y-bin to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_TOP:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_TOP";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetTop(atoi(temp));
      *result = "Change CCD top to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_LEFT:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough argumnets to SET_LEFT";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetLeft(atoi(temp));
      *result = "Change CCD left to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_HEIGHT:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_HEIGHT";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetHeight(atoi(temp));
      *result = "Change CCD height to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_WIDTH:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_WIDTH";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetWidth(atoi(temp));
      *result = "Change CCD height to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_TEMP:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_TEMP";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetCCDTemp(atof(temp));
      *result = "Change CCD temperature setpoint to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_NUMFILTERS:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_TEMP";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myFilterWheel->SetNumFilters(atoi(temp));
      *result = "Changed number of filter positions to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_GAIN:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_GAIN";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetGain(atof(temp));
      *result = "Change CCD gain setpoint to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;      
    }
    case SET_DATA_DIR:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_DATA_DIR";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myDirectories.SetDataPath(sout.str());
      *result = "Change Data directory to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_PNG_DIR:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_PNG_DIR";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myDirectories.SetPNGPath(sout.str());
      *result = "Change PNG directory to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_MPG_DIR:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_MPG_DIR";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myDirectories.SetMoviesPath(sout.str());
      *result = "Change MPG directory to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_QL_DIR:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_QL_DIR";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myDirectories.SetQuickLookPath(sout.str());
      *result = "Change quicklook directory to " + sout.str();
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_DO_PNG:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_DO_PNG";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      if(atoi(temp) == 0)
      {
        mySystem.myConfiguration.myDirectories.SetDoPNG(false);
        *result = "Change doPNG to FALSE";
      }
      else
      {
        mySystem.myConfiguration.myDirectories.SetDoPNG(true);
        *result = "Change doPNG to TRUE";
      }
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_DO_MPG:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_DO_MPG";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      if(atoi(temp) == 0)
      {
        mySystem.myConfiguration.myDirectories.SetDoMovies(false);
        *result = "Change doMovies to FALSE";
      }
      else
      {
        mySystem.myConfiguration.myDirectories.SetDoMovies(true);
        *result = "Change doMovies to TRUE";
      }
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_DO_QL:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_DO_QL";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      if(atoi(temp) == 0)
      {
        mySystem.myConfiguration.myDirectories.SetDoQuickLook(false);
        *result = "Change doQuickLook to FALSE";
      }
      else
      {
        mySystem.myConfiguration.myDirectories.SetDoQuickLook(true);
        *result = "Change doQuickLook to TRUE";
      }
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_CCD_AUTO_TEMP:
    {
    	temp = strtok(NULL, ":");
     	if(temp == NULL)
      {
      	*result = "ERROR: Not enough arguments to SET_CCD_AUTO_TEMP";
       	break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      if(atoi(temp) == 0)
      {
        // If temperature regulation is not set to auto, then we want to turn the
        // fan on all the time
      	mySystem.myConfiguration.myCamera->SetAutoTemp(false);
        *result = "Change AutoTempRegulation to FALSE";
        mySystem.myConfiguration.myCamera->SetTemperatureRegulationOn();
      }
      else
      {
        // If auto temperature regulation, turn the fan off and allow
        // the system to decide when to turn it back on
      	mySystem.myConfiguration.myCamera->SetAutoTemp(true);
      	*result = "Change AutoTempRegulation to TRUE";
       	mySystem.myConfiguration.myCamera->SetTemperatureRegulationOff();
      }
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_FILTERWHEELTYPE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_FILTERHWEELTYPE";
        break;
      }

      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myFilterWheelType = atoi(temp);
      *result = "Change filterwheel type to " + sout.str();
      pthread_mutex_unlock(&mutex);

      // Open the filterwheel
      if(mySystem.myConfiguration.myFilterWheelType == SBIGCFWL)
      {
        // Save the needed variables to temporary holders
        int tempNumFilters, tempShutterPos;

        tempNumFilters = mySystem.myConfiguration.myFilterWheel->myNumFilters;
        tempShutterPos = mySystem.myConfiguration.myFilterWheel->myShutterPos;

        // close the old filterwheel
        mySystem.myConfiguration.myFilterWheel->~FilterWheel();
        delete(mySystem.myConfiguration.myFilterWheel);

        // open the new one
        mySystem.myConfiguration.myFilterWheel = new FilterWheelCFWL();
        mySystem.myConfiguration.myFilterWheel->myNumFilters = tempNumFilters;
        mySystem.myConfiguration.myFilterWheel->myShutterPos = tempShutterPos;
      }
      if(mySystem.myConfiguration.myFilterWheelType == FLIFWL)
      {
        // Save the needed variables to temporary holders
        int tempNumFilters, tempShutterPos;

        tempNumFilters = mySystem.myConfiguration.myFilterWheel->myNumFilters;
        tempShutterPos = mySystem.myConfiguration.myFilterWheel->myShutterPos;

        // close the old filterwheel
        mySystem.myConfiguration.myFilterWheel->~FilterWheel();
        delete(mySystem.myConfiguration.myFilterWheel);

        // open the new one
        mySystem.myConfiguration.myFilterWheel = new FilterWheelFLI();
        mySystem.myConfiguration.myFilterWheel->myNumFilters = tempNumFilters;
        mySystem.myConfiguration.myFilterWheel->myShutterPos = tempShutterPos;
      }
      if(mySystem.myConfiguration.myFilterWheelType == KEOFW)
      {
        // Save the needed variables to temporary holders
        int tempNumFilters, tempShutterPos;
	char tempFWDevice[256];
    	char tempTempControlDevice[256];

        tempNumFilters = mySystem.myConfiguration.myFilterWheel->myNumFilters;
        tempShutterPos = mySystem.myConfiguration.myFilterWheel->myShutterPos;
//	strcpy(tempFWDevice,  mySystem.myConfiguration.myFilterWheel->myFWDevice);
//	strcpy(tempTempControlDevice,  mySystem.myConfiguration.myFilterWheel->myTempControlDevice);

        // close the old filterwheel
        mySystem.myConfiguration.myFilterWheel->~FilterWheel();
        delete(mySystem.myConfiguration.myFilterWheel);

        // open the new one
        mySystem.myConfiguration.myFilterWheel = new FilterWheelKEO();
        mySystem.myConfiguration.myFilterWheel->myNumFilters = tempNumFilters;
        mySystem.myConfiguration.myFilterWheel->myShutterPos = tempShutterPos;
//	strcpy(mySystem.myConfiguration.myFilterWheel->myFWDevice, tempFWDevice);
//	strcpy(mySystem.myConfiguration.myFilterWheel->myTempControlDevice, tempTempControlDevice);
	mySystem.myConfiguration.myFilterWheel->SetUp_FilterWheel();
	mySystem.myConfiguration.myFilterWheel->CloseShutter();
      }
      if(mySystem.myConfiguration.myFilterWheelType == KEOSER)
      {
        // Save the needed variables to temporary holders
        int tempNumFilters, tempShutterPos;
        
        tempNumFilters = mySystem.myConfiguration.myFilterWheel->myNumFilters;
        tempShutterPos = mySystem.myConfiguration.myFilterWheel->myShutterPos;

        // close the old filterwheel
        mySystem.myConfiguration.myFilterWheel->~FilterWheel();
        delete(mySystem.myConfiguration.myFilterWheel);

        // open the new one
        mySystem.myConfiguration.myFilterWheel = new FilterWheelSERIAL();
        mySystem.myConfiguration.myFilterWheel->myNumFilters = tempNumFilters;
        mySystem.myConfiguration.myFilterWheel->myShutterPos = tempShutterPos;
        mySystem.myConfiguration.myFilterWheel->SetUp_FilterWheel();
        mySystem.myConfiguration.myFilterWheel->CloseShutter();
      }

      break;
    }
    case SET_FILTERINFO:
    {
      sout.str() = "";
      // read in the parameters to set the filter info;
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_FILTERINFO";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      int myCurrentFilter;
      myCurrentFilter = atoi(temp);
      *result = "Filterinfo params: Filter: " + sout.str();
      pthread_mutex_unlock(&mutex);

      sout.str() = "";
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_FILTERINFO";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      mySystem.myConfiguration.mySequencer.SetFilterName(myCurrentFilter, str);
      *result += ", Name: " + str;
      pthread_mutex_unlock(&mutex);

      sout.str() = "";
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_FILTERINFO";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      mySystem.myConfiguration.mySequencer.SetEmissionHeight(myCurrentFilter, atof(temp));
      *result += ", Height: " + str;
      pthread_mutex_unlock(&mutex);

      sout.str() = "";
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_FILTERINFO";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      mySystem.myConfiguration.mySequencer.SetExposureTime(myCurrentFilter, atof(temp));
      *result += ", Exp. Time: " + str;
      pthread_mutex_unlock(&mutex);

      sout.str() = "";
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_FILTERINFO";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      if(atoi(temp))
      {
        // This is the shutter position, so set it
        mySystem.myConfiguration.mySequencer.SetShutterPosition(myCurrentFilter);
      }
      *result += ", Shutter: " + str;
      pthread_mutex_unlock(&mutex);
      break;
    }
    case SET_SUBSEQUENCE:
    {
      sout.str() = "";
      // read in the parameters to set the filter info;
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SUBSEQUENCE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      int myCurrentSubSequence;
      myCurrentSubSequence = atoi(temp);    // sequences are 0 indexed
      *result = "SubSequence params: SubSequence: " + sout.str();
      pthread_mutex_unlock(&mutex);

      sout.str() = "";
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SUBSEQUENCE";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      int myNumberInSubSequence;
      myNumberInSubSequence = atoi(temp);
      mySystem.myConfiguration.mySequencer.SetNumberInSubSequence(myCurrentSubSequence, myNumberInSubSequence);
      *result += ", Number in subsequence: " + str;
      pthread_mutex_unlock(&mutex);

      for(int i = 0; i < myNumberInSubSequence; i++)
      {
        sout.str() = "";
        temp = strtok(NULL, ":");
        if(temp == NULL)
        {
          *result = "ERROR: Not enough arguments to SET_SUBSEQUENCE";
          break;
        }
        pthread_mutex_lock(&mutex);
        str = temp;
        mySystem.myConfiguration.mySequencer.SetSubSequence(myCurrentSubSequence, i, atoi(temp));
        *result += ", " + str;
        pthread_mutex_unlock(&mutex);
      }
      break;
    }
    case SET_SUBSEQUENCEORDER:
    {
      sout.str() = "";
      // read in the parameters to set the filter info;
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_SUBSEQUENCEORDER";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      int myNumberInSubSequenceOrder;
      myNumberInSubSequenceOrder = atoi(temp);
      mySystem.myConfiguration.mySequencer.SetNumberInSubSequenceOrder(myNumberInSubSequenceOrder);
      *result += "SubSequenceOrder params: Number in subsequence order: " + str;
      pthread_mutex_unlock(&mutex);

      for(int i = 0; i < myNumberInSubSequenceOrder; i++)
      {
        sout.str() = "";
        temp = strtok(NULL, ":");
        if(temp == NULL)
        {
          *result = "ERROR: Not enough arguments to SET_SUBSEQUENCEORDER";
          break;
        }
        pthread_mutex_lock(&mutex);
        str = temp;
        mySystem.myConfiguration.mySequencer.SetSubSequenceOrder(i, atoi(temp));
        *result += ", " + str;
        pthread_mutex_unlock(&mutex);
      }

      // Now set the sequence
      mySystem.myConfiguration.mySequencer.CreateSequence();
      
      break;
    }
    case SET_DARKSEQUENCE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        * result = "ERROR: Not enough arguments to SET_DARKSEQUENCE";
        break;
      }
      pthread_mutex_lock(&mutex);
      if(atoi(temp) == 0)
      {
        mySystem.myConfiguration.mySequencer.SetDoDark(false);
        *result = "Set darksequencing to FALSE";
      } else
      {
        mySystem.myConfiguration.mySequencer.SetDoDark(true);
        *result = "Set darksequencing to TRUE";
      }
      pthread_mutex_unlock(&mutex);

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_DARKSEQUENCE";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      mySystem.myConfiguration.mySequencer.SetDarkFrequency(atoi(temp));
      *result += ", Frequency to " + str;
      pthread_mutex_unlock(&mutex);

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_DARKSEQUENCE";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      mySystem.myConfiguration.mySequencer.SetDarkExposureTime(atof(temp));
      *result += ", Exposure time to " + str;
      pthread_mutex_unlock(&mutex);
      
      break;
    }
    case SET_PERIODICITY:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        * result = "ERROR: Not enough arguments to SET_PERIODICITY";
        break;
      }
      pthread_mutex_lock(&mutex);
      if(atoi(temp) == 0)
      {
        mySystem.myConfiguration.mySequencer.SetDoPeriodicity(false);
        *result = "Set periodicity to FALSE";
      } else
      {
        mySystem.myConfiguration.mySequencer.SetDoPeriodicity(true);
        *result = "Set periodicity to TRUE";
      }
      pthread_mutex_unlock(&mutex);

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_PERIODICITY";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      mySystem.myConfiguration.mySequencer.SetPeriodicityFrequency(atoi(temp));
      *result += ", Frequency to " + str;
      pthread_mutex_unlock(&mutex);

      break;
    }
    case SET_CAMERATYPE:
    {
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to SET_CAMERATYPE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCameraType = atoi(temp);
      *result = "Change camera type to " + sout.str();
      pthread_mutex_unlock(&mutex);
      
      // Open up the camera
//      if(mySystem.myConfiguration.myCameraType == SBIGUSB)
      if(mySystem.myConfiguration.myCameraType == SBIGUNIV)
      {
        // Save the needed variables to temporary holders
        int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
        double tempTemp;
        float tempGain;
        bool tempRegulation;

        tempXBin = mySystem.myConfiguration.myCamera->GetBinX();
        tempYBin = mySystem.myConfiguration.myCamera->GetBinY();
        tempTop = mySystem.myConfiguration.myCamera->GetTop();
        tempLeft = mySystem.myConfiguration.myCamera->GetLeft();
        tempHeight = mySystem.myConfiguration.myCamera->GetHeight();
        tempWidth = mySystem.myConfiguration.myCamera->GetWidth();
        tempTemp = mySystem.myConfiguration.myCamera->GetCCDSetTemp();
        tempGain = mySystem.myConfiguration.myCamera->GetGain();
        tempRegulation = mySystem.myConfiguration.myCamera->GetAutoTemp();
        
        // close the old camera
        mySystem.myConfiguration.myCamera->~Camera();
        delete(mySystem.myConfiguration.myCamera);

        // open the new one
//        mySystem.myConfiguration.myCamera = new CameraSBIG();
        mySystem.myConfiguration.myCamera = new CameraSBIGUniv();
        mySystem.myConfiguration.myCamera->SetBinX(tempXBin);
        mySystem.myConfiguration.myCamera->SetBinY(tempYBin);
        mySystem.myConfiguration.myCamera->SetTop(tempTop);
        mySystem.myConfiguration.myCamera->SetLeft(tempLeft);
        mySystem.myConfiguration.myCamera->SetHeight(tempHeight);
        mySystem.myConfiguration.myCamera->SetWidth(tempWidth);
        mySystem.myConfiguration.myCamera->SetCCDTemp(tempTemp);
        mySystem.myConfiguration.myCamera->SetGain(tempGain);
        mySystem.myConfiguration.myCamera->SetAutoTemp(tempRegulation);
      }
      
      if(mySystem.myConfiguration.myCameraType == ANDOR)
      {
        // Save the needed variables to temporary holders
	      int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
	      double tempTemp;
	      float tempGain;
	      bool tempRegulation;

	      tempXBin = mySystem.myConfiguration.myCamera->GetBinX();
	      tempYBin = mySystem.myConfiguration.myCamera->GetBinY();
	      tempTop = mySystem.myConfiguration.myCamera->GetTop();
	      tempLeft = mySystem.myConfiguration.myCamera->GetLeft();
	      tempHeight = mySystem.myConfiguration.myCamera->GetHeight();
	      tempWidth = mySystem.myConfiguration.myCamera->GetWidth();
	      tempTemp = mySystem.myConfiguration.myCamera->GetCCDSetTemp();
	      tempGain = mySystem.myConfiguration.myCamera->GetGain();
	      tempRegulation = mySystem.myConfiguration.myCamera->GetAutoTemp();
        
        // close the old camera
	      mySystem.myConfiguration.myCamera->~Camera();
	      delete(mySystem.myConfiguration.myCamera);

        // open the new one
	      mySystem.myConfiguration.myCamera = new CameraAndor();
	      mySystem.myConfiguration.myCamera->SetBinX(tempXBin);
	      mySystem.myConfiguration.myCamera->SetBinY(tempYBin);
	      mySystem.myConfiguration.myCamera->SetTop(tempTop);
	      mySystem.myConfiguration.myCamera->SetLeft(tempLeft);
	      mySystem.myConfiguration.myCamera->SetHeight(tempHeight);
	      mySystem.myConfiguration.myCamera->SetWidth(tempWidth);
	      mySystem.myConfiguration.myCamera->SetCCDTemp(tempTemp);
	      mySystem.myConfiguration.myCamera->SetGain(tempGain);
	      mySystem.myConfiguration.myCamera->SetAutoTemp(tempRegulation);
      }
      
      if(mySystem.myConfiguration.myCameraType == FLIPROLINE)
      {
        // Save the needed variables to temporary holders
	      int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
	      double tempTemp;
	      float tempGain;
	      bool tempRegulation;

	      tempXBin = mySystem.myConfiguration.myCamera->GetBinX();
	      tempYBin = mySystem.myConfiguration.myCamera->GetBinY();
	      tempTop = mySystem.myConfiguration.myCamera->GetTop();
	      tempLeft = mySystem.myConfiguration.myCamera->GetLeft();
	      tempHeight = mySystem.myConfiguration.myCamera->GetHeight();
	      tempWidth = mySystem.myConfiguration.myCamera->GetWidth();
	      tempTemp = mySystem.myConfiguration.myCamera->GetCCDSetTemp();
	      tempGain = mySystem.myConfiguration.myCamera->GetGain();
	      tempRegulation = mySystem.myConfiguration.myCamera->GetAutoTemp();
        
        // close the old camera
	      mySystem.myConfiguration.myCamera->~Camera();
	      delete(mySystem.myConfiguration.myCamera);

        // open the new one
	      mySystem.myConfiguration.myCamera = new CameraFLI();
	      mySystem.myConfiguration.myCamera->SetBinX(tempXBin);
	      mySystem.myConfiguration.myCamera->SetBinY(tempYBin);
	      mySystem.myConfiguration.myCamera->SetTop(tempTop);
	      mySystem.myConfiguration.myCamera->SetLeft(tempLeft);
	      mySystem.myConfiguration.myCamera->SetHeight(tempHeight);
	      mySystem.myConfiguration.myCamera->SetWidth(tempWidth);
	      mySystem.myConfiguration.myCamera->SetCCDTemp(tempTemp);
	      mySystem.myConfiguration.myCamera->SetGain(tempGain);
	      mySystem.myConfiguration.myCamera->SetAutoTemp(tempRegulation);
      }
      
/*      if(mySystem.myConfiguration.myCameraType == PVCAM)
      {
        // Save the needed variables to temporary holders
	      int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
	      double tempTemp;
	      float tempGain;
	      bool tempRegulation;

	      tempXBin = mySystem.myConfiguration.myCamera->GetBinX();
	      tempYBin = mySystem.myConfiguration.myCamera->GetBinY();
	      tempTop = mySystem.myConfiguration.myCamera->GetTop();
	      tempLeft = mySystem.myConfiguration.myCamera->GetLeft();
	      tempHeight = mySystem.myConfiguration.myCamera->GetHeight();
	      tempWidth = mySystem.myConfiguration.myCamera->GetWidth();
	      tempTemp = mySystem.myConfiguration.myCamera->GetCCDSetTemp();
	      tempGain = mySystem.myConfiguration.myCamera->GetGain();
	      tempRegulation = mySystem.myConfiguration.myCamera->GetAutoTemp();
        
        // close the old camera
	      mySystem.myConfiguration.myCamera->~Camera();
	      delete(mySystem.myConfiguration.myCamera);

        // open the new one
	      mySystem.myConfiguration.myCamera = new CameraPVCAM();
	      mySystem.myConfiguration.myCamera->SetBinX(tempXBin);
	      mySystem.myConfiguration.myCamera->SetBinY(tempYBin);
	      mySystem.myConfiguration.myCamera->SetTop(tempTop);
	      mySystem.myConfiguration.myCamera->SetLeft(tempLeft);
	      mySystem.myConfiguration.myCamera->SetHeight(tempHeight);
	      mySystem.myConfiguration.myCamera->SetWidth(tempWidth);
	      mySystem.myConfiguration.myCamera->SetCCDTemp(tempTemp);
	      mySystem.myConfiguration.myCamera->SetGain(tempGain);
	      mySystem.myConfiguration.myCamera->SetAutoTemp(tempRegulation);
      }
*/
/*      if(mySystem.myConfiguration.myCameraType == APOGEE)
      {
        // Save the needed variables to temporary holders
              int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
              double tempTemp;
              float tempGain;
              bool tempRegulation;

              tempXBin = mySystem.myConfiguration.myCamera->GetBinX();
              tempYBin = mySystem.myConfiguration.myCamera->GetBinY();
              tempTop = mySystem.myConfiguration.myCamera->GetTop();
              tempLeft = mySystem.myConfiguration.myCamera->GetLeft();
              tempHeight = mySystem.myConfiguration.myCamera->GetHeight();
              tempWidth = mySystem.myConfiguration.myCamera->GetWidth();
              tempTemp = mySystem.myConfiguration.myCamera->GetCCDSetTemp();
              tempGain = mySystem.myConfiguration.myCamera->GetGain();
              tempRegulation = mySystem.myConfiguration.myCamera->GetAutoTemp();
        
        // close the old camera
              mySystem.myConfiguration.myCamera->~Camera();
              delete(mySystem.myConfiguration.myCamera);

        // open the new one
              mySystem.myConfiguration.myCamera = new CameraAPOGEE();
              mySystem.myConfiguration.myCamera->SetBinX(tempXBin);
              mySystem.myConfiguration.myCamera->SetBinY(tempYBin);
              mySystem.myConfiguration.myCamera->SetTop(tempTop);
              mySystem.myConfiguration.myCamera->SetLeft(tempLeft);
              mySystem.myConfiguration.myCamera->SetHeight(tempHeight);
              mySystem.myConfiguration.myCamera->SetWidth(tempWidth);
              mySystem.myConfiguration.myCamera->SetCCDTemp(tempTemp);
              mySystem.myConfiguration.myCamera->SetGain(tempGain);
              mySystem.myConfiguration.myCamera->SetAutoTemp(tempRegulation);
      }*/
      break;
    }
    case GET_SCHEDULE_MODE:
    {
      // return the schedule mode
      sout << mySystem.myConfiguration.mySchedule.GetScheduleMode();
      *result = sout.str();
      break;
    }
    case GET_SUNANGLE:
    {
      // return the sun angle mask
      sout << mySystem.myConfiguration.mySchedule.GetSunAngle();
      *result = sout.str();
      break;
    }
    case GET_MOONANGLE:
    {
      // return the moon angle mask
      sout << mySystem.myConfiguration.mySchedule.GetMoonAngle();
      *result = sout.str();
      break;
    }
    case GET_DOMOONSET:
    {
      // return the moonset flag
      sout << mySystem.myConfiguration.mySchedule.GetMoonSet();
      *result = sout.str();
      break;
    }
    case GET_STARTTIME:
    {
      // return the starttime
      struct tm StartTime = mySystem.myConfiguration.mySchedule.GetStartTime();
      sout << mktime(&StartTime);
      *result = sout.str();
      break;
    }
    case GET_STOPTIME:
    {
      // return the stoptime
      struct tm StopTime = mySystem.myConfiguration.mySchedule.GetStopTime();
      sout << mktime(&StopTime);
      *result = sout.str();
      break;
    }
    case GET_SITENAME:
    {
      // return the sitename
      sout << mySystem.myConfiguration.mySite.GetName();
      *result = sout.str();
      break;
    }
    case GET_SITEABBR:
    {
      // return the abbreviation
      sout << mySystem.myConfiguration.mySite.GetAbbreviation();
      *result = sout.str();
      break;
    }
    case GET_LATITUDE:
    {
      // returns the latitude
      sout << mySystem.myConfiguration.mySite.GetLatitude();
      *result = sout.str();
      break;
    }
    case GET_SCHEDULELATITUDE:
    {
      // returns the latitude
      sout << mySystem.myConfiguration.mySchedule.GetLatitude();
      *result = sout.str();
      break;
    }
    case GET_LONGITUDE:
    {
      // returns the longitude
      sout << mySystem.myConfiguration.mySite.GetLongitude();
      *result = sout.str();
      break;
    }
    case GET_SCHEDULELONGITUDE:
    {
      // returns the longitude
      sout << mySystem.myConfiguration.mySchedule.GetLongitude();
      *result = sout.str();
      break;
    }
    case GET_ALTITUDE:
    {
      // return the altitude
      sout << mySystem.myConfiguration.mySite.GetHeight();
      *result = sout.str();
      break;
    }
    case GET_XBIN:
    {
      // return the x-binning factor
      sout << mySystem.myConfiguration.myCamera->GetBinX();
      *result = sout.str();
      break;
    }
    case GET_YBIN:
    {
      // return the y-binning factor
      sout << mySystem.myConfiguration.myCamera->GetBinY();
      *result = sout.str();
      break;
    }
    case GET_TOP:
    {
      // return the top of the imaged CCD
      sout << mySystem.myConfiguration.myCamera->GetTop();
      *result = sout.str();
      break;
    }
    case GET_LEFT:
    {
      // return the left of the imaged CCD
      sout << mySystem.myConfiguration.myCamera->GetLeft();
      *result = sout.str();
      break;
    }
    case GET_HEIGHT:
    {
      // return the height of the imaged CCD
      sout << mySystem.myConfiguration.myCamera->GetHeight();
      *result = sout.str();
      break;
    }
    case GET_WIDTH:
    {
      // return the width of the imaged CCD
      sout << mySystem.myConfiguration.myCamera->GetWidth();
      *result = sout.str();
      break;
    }
    case GET_TEMP:
    {
      // return the CCD temperature setpoint
      sout << mySystem.myConfiguration.myCamera->GetCCDSetTemp();
      *result = sout.str();
      break;
    }
    case GET_GAIN:
    {
      // return the CCD gain
      sout << mySystem.myConfiguration.myCamera->GetGain();
      *result = sout.str();
      break;
    }
    case GET_FILTERWHEELTYPE:
    {
      // return the filterwheel type
      sout << mySystem.myConfiguration.myFilterWheelType;
      *result = sout.str();
      break;
    }
    case GET_NUMFILTERS:
    {
      // return the number of filter positions
      cout << "TEMP2: " << mySystem.myConfiguration.myFilterWheel->GetNumFilters() << "::";
      sout << mySystem.myConfiguration.myFilterWheel->GetNumFilters();
      *result = sout.str();
      break;
    }
    case GET_FILTERINFO:
    {
      // returns the filter information for the indicated filter
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to GET_FILTERINFO";
        break;
      }
      pthread_mutex_lock(&mutex);
      int myFilter;
      myFilter = atoi(temp);
      pthread_mutex_unlock(&mutex);

      bool isShutter;
      if(myFilter == mySystem.myConfiguration.mySequencer.GetShutterPosition())
        isShutter = true;
      else
        isShutter = false;

      sout << mySystem.myConfiguration.mySequencer.GetFilterName(myFilter) << ":" << mySystem.myConfiguration.mySequencer.GetEmissionHeight(myFilter) << ":" << mySystem.myConfiguration.mySequencer.GetExposureTime(myFilter) << ":" << isShutter;
      *result = sout.str();
      break;
    }
    case GET_SUBSEQUENCE:
    {
      // returns the subsequence information for the indicated subsequence
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to GET_SUBSEQUENCE";
        break;
      }
      pthread_mutex_lock(&mutex);
      int myRequestedSubSequence;
      myRequestedSubSequence = atoi(temp);    // subsequences are 0 indexed
      pthread_mutex_unlock(&mutex);

      // Find out how many filters are in the subsequence
      int myNumberInRequestedSubSequence;
      myNumberInRequestedSubSequence = mySystem.myConfiguration.mySequencer.GetNumberInSubSequence(myRequestedSubSequence);

      // Begin forming the return string
      sout << myNumberInRequestedSubSequence;

      // Grab each of the filters in the subsequence
      for(int i = 0; i < myNumberInRequestedSubSequence; i++)
        sout << ":" << mySystem.myConfiguration.mySequencer.GetSubSequence(myRequestedSubSequence, i);

      // Form the result string
      *result = sout.str();
      
      break;
    }
    case GET_SUBSEQUENCEORDER:
    {
      // returns the subsequence order information
      // Find out how many subsequences are in the subsequence
      int myNumberInSequenceOrder;
      myNumberInSequenceOrder = mySystem.myConfiguration.mySequencer.GetNumberInSubSequenceOrder();

      // Begin forming the return string
      sout << myNumberInSequenceOrder;

      // Grab each of the filters in the subsequence
      for(int i = 0; i < myNumberInSequenceOrder; i++)
        sout << ":" << mySystem.myConfiguration.mySequencer.GetSubSequenceOrder(i);

      // Form the result string
      *result = sout.str();
      break;
    }
    case GET_DARKSEQUENCE:
    {
      // Return the dark sequencing information
      sout << mySystem.myConfiguration.mySequencer.GetDoDark() << ":" << mySystem.myConfiguration.mySequencer.GetDarkFrequency() << ":" << mySystem.myConfiguration.mySequencer.GetDarkExposureTime();
      *result = sout.str();
      break;
    }
    case GET_PERIODICITY:
    {
      // Return the periodicity information
      sout << mySystem.myConfiguration.mySequencer.GetDoPeriodicity() << ":" << mySystem.myConfiguration.mySequencer.GetPeriodicityFrequency();
      *result = sout.str();
      break;
    }
    case GET_CAMERATYPE:
    {
      // return the camera type
      sout << mySystem.myConfiguration.myCameraType;
      *result = sout.str();
      break;
    }
    case GET_DATA_DIR:
    {
      // return the data directory
      sout << mySystem.myConfiguration.myDirectories.GetDataPath();
      *result = sout.str();
      break;
    }
    case GET_PNG_DIR:
    {
      // return the png directory
      sout << mySystem.myConfiguration.myDirectories.GetPNGPath();
      *result = sout.str();
      break;
    }
    case GET_MPG_DIR:
    {
      // return the mpg directory
      sout << mySystem.myConfiguration.myDirectories.GetMoviesPath();
      *result = sout.str();
      break;
    }
    case GET_QL_DIR:
    {
      // return the quicklook directory
      sout << mySystem.myConfiguration.myDirectories.GetQuickLookPath();
      *result = sout.str();
      break;
    }
    case GET_DO_PNG:
    {
      // return the dopng variable
      sout << mySystem.myConfiguration.myDirectories.GetDoPNG();
      *result = sout.str();
      break;
    }
    case GET_DO_MPG:
    {
      // return the dompg variable
      sout << mySystem.myConfiguration.myDirectories.GetDoMovies();
      *result = sout.str();
      break;
    }
    case GET_DO_QL:
    {
      // return the doql variable
      sout << mySystem.myConfiguration.myDirectories.GetDoQuickLook();
      *result = sout.str();
      break;
    }
    case GET_ACTUAL_TEMP:
    {
      // return the CCD temperature 
      sout << mySystem.myConfiguration.myCamera->GetCCDTemp();
      *result = sout.str();
      break;
    }
    case GET_CCD_AUTO_TEMP:
    {
    	// return the CCD auto temperature regulation parameter
     	sout << mySystem.myConfiguration.myCamera->GetAutoTemp();
     	*result = sout.str();
	break;
    }
    case CAPTURE_IMAGE:
    {
      sout.str() = "";
      // read in the parameters and take a picture;
      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to CAPTURE_IMAGE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;
      mySystem.myConfiguration.myCamera->SetExpTime(atof(temp));
      *result = "Capture params: Exposure time: " + sout.str();
      pthread_mutex_unlock(&mutex);

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to CAPTURE_IMAGE";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      tempBool = atoi(temp);
      if(tempBool == 0)
{
        mySystem.myConfiguration.myCamera->SetDark(false);
      // Open the filterwheel shutter
      mySystem.myConfiguration.myFilterWheel->OpenShutter();
      sleep(1);
      } else {
        mySystem.myConfiguration.myCamera->SetDark(true);
      // Open the filterwheel shutter
      mySystem.myConfiguration.myFilterWheel->CloseShutter();
      sleep(1);
}
      *result = *result + ", dark: " + str;
      pthread_mutex_unlock(&mutex);

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to CAPTURE_IMAGE";
        break;
      }
      pthread_mutex_lock(&mutex);
      sout << temp;

      mySystem.myConfiguration.myFilterWheel->MoveToPosition(atoi(temp));
      // Wait for the filterwheel to get there
      while(mySystem.myConfiguration.myFilterWheel->GetCurrentPosition() == -999)
      {
        // Still moving
        sleep(1);
      }

      *result = "Capture params: Filter position: " + sout.str();
      pthread_mutex_unlock(&mutex);

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to CAPTURE_IMAGE";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = string(temp) + ".tif";
      filename = string(temp);
      *result = *result + ", Filename: " + str;
      pthread_mutex_unlock(&mutex);

      // Add the directory stub for today
      mySystem.myConfiguration.myDirectories.SetTodaysDir();
      filename = mySystem.myConfiguration.myDirectories.GetTodaysDir() + filename;

      temp = strtok(NULL, ":");
      if(temp == NULL)
      {
        *result = "ERROR: Not enough arguments to CAPTURE_IMAGE";
        break;
      }
      pthread_mutex_lock(&mutex);
      str = temp;
      tempBool = atoi(temp);
      *result = *result + ", savePNG: " + str;
      pthread_mutex_unlock(&mutex);

      // Take the image
      mySystem.myConfiguration.myLogger.WriteSystemLogString(*result);
      mySystem.myConfiguration.myCamera->CaptureImage(mySystem.myConfiguration.myImage, mySystem.myConfiguration.myDirectories.GetDataPath() + filename, true);

      // Close the filterwheel shutter
      mySystem.myConfiguration.myFilterWheel->CloseShutter();

      // Save PNG image
      if(tempBool)
        mySystem.myConfiguration.myImage->SavePNG(mySystem.myConfiguration.myDirectories.GetPNGPath() + filename);
   
      break;
    }
    default:
    {
      *result = "Command not recognized.";
      break;
    }
  }
  mySystem.myConfiguration.myLogger.WriteSystemLogString(*result);

  // Save the XML configuration file
  mySystem.myConfiguration.SaveXML(mySystem.myConfiguration.myDirectories.GetQuickLookPath() + mySystem.myXMLFilename);

}

void *control_function(void *threadid)
{
  char tempChar[8];
  double myTemp;
  string str;
  char time_string[64];
  struct tm temp_time;
  bool done = false;
	
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Starting control thread");
  mySystem.myControlRunning = true;

  mySystem.myConfiguration.mySchedule.SetTodaysRiseSetTimes();
 
  if(mySystem.myConfiguration.myCamera->GetAutoTemp())
	mySystem.myConfiguration.myCamera->SetTemperatureRegulationOff();

  while(mySystem.myIsRunning) {
     // Based on the schedule mode, we will either:
     //   1) Sit patiently and do nothing (SCHEDULE_OFF)
     //   2) Run once (SCHEDULE_MANUAL)
     //   3) Run continuously (SCHEDULE_AUTO);
    
     if(mySystem.myConfiguration.mySchedule.GetScheduleMode() != SCHEDULE_OFF)
     {
       mySystem.myConfiguration.mySchedule.SetTodaysRiseSetTimes();
              
       if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
       {
         // Auto sequencing, so get the start/stop time
         mySystem.myConfiguration.mySchedule.AutoSetTimes();
       }

       // If the system is in manual mode, the start/stop times are already set

       // Now we need to wait for the start time
       // We need to keep an eye out for:
       //   1) mySystem.myIsRunning = false
       //   2) mode being switch to off
       
	// JJM mods to test filterwheel homing (removed on 4/23/09)
// 	mySystem.myConfiguration.myFilterWheel->~FilterWheel();
// 	if(mySystem.myConfiguration.myFilterWheelType == SBIGCFWL)
// 		mySystem.myConfiguration.myFilterWheel = new FilterWheelCFWL();
// 	if(mySystem.myConfiguration.myFilterWheelType == FLIFWL)
// 		mySystem.myConfiguration.myFilterWheel = new FilterWheelFLI();

       pthread_mutex_lock(&mutex);
       temp_time = mySystem.myConfiguration.mySchedule.GetStartTime();
       pthread_mutex_unlock(&mutex);
         
       strftime(time_string, 64, "%a %b %d %H:%M:%S %Y", &temp_time);
       str = string("Waiting for ") + time_string + string("...");
       mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
       
       time_t now = time(NULL);
       done = false;
       while(mySystem.myIsRunning &&
             mySystem.myConfiguration.mySchedule.GetScheduleMode() != SCHEDULE_OFF &&
             !done)
       {
         // Update the current time
         now = time(NULL);
         pthread_mutex_lock(&mutex);
         temp_time = mySystem.myConfiguration.mySchedule.GetStartTime();
         pthread_mutex_unlock(&mutex);
	 
	 // Check to see if the time is less than 60 minutes to the
	 // start of the sequence.
	 if((mktime(&temp_time) - now < 3600.0) &&
	    (!mySystem.myConfiguration.myCamera->GetTemperatureRegulation()))
	 {
 		mySystem.myConfiguration.myCamera->SetTemperatureRegulationOn();
	 }

	// Write out the CCD temperature every 5 minutes
	if(now % 300 == 0)
	{
		mySystem.myConfiguration.myCamera->ReadTemperatureStatus();
		myTemp = mySystem.myConfiguration.myCamera->GetCCDTemp();
		sprintf(tempChar, "%7.2f", myTemp);
		str = string("Current CCD temperature: ") + tempChar;
		mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
	}
         
         if(now >= mktime(&temp_time))
         {
           // We are past the start time
           done = true;
         } else
         {
           // Wait for a bit
           sleep(1);
         }
       }

       // check which of the cases for exiting the wait loop were satisfied
       if(done)
       {
         // For some reason, the camera has been losing settings.  Reread the camera settings in
	 // ADDED 20-Sep-2006 by Jonathan J. Makela (jmakela@uiuc.edu)
	 // Removed 29-Jan-2007 to avoid date ambiguity
	 //	 mySystem.myConfiguration.LoadConfiguration(mySystem.myConfFileName);
       
         // We reached the time so need to do sequencing         
         mySystem.myConfiguration.myLogger.WriteSystemLogString("Starting sequence");
         pthread_mutex_lock(&mutex);
         temp_time = mySystem.myConfiguration.mySchedule.GetStopTime();
         pthread_mutex_unlock(&mutex);

	// First open the shutter (this only does something on KEO filter wheels);
	mySystem.myConfiguration.myFilterWheel->OpenShutter();
	
         // Start the sequence (e.g., reset variables)
         mySystem.myConfiguration.mySequencer.StartSequence();
	 
         // Perform the sequencing
	 tm temp_StartTime = mySystem.myConfiguration.mySchedule.GetStartTime();
         while(now < mktime(&temp_time) &&
               mySystem.myConfiguration.mySchedule.GetScheduleMode() != SCHEDULE_OFF  &&
               mySystem.myIsRunning  &&
               now >= mktime(&temp_StartTime))
         {

           // Capture the next image in the sequence
           mySystem.myConfiguration.mySequencer.CaptureNextImage();
	   
           // Update the stop time
           pthread_mutex_lock(&mutex);
           temp_time = mySystem.myConfiguration.mySchedule.GetStopTime();
           pthread_mutex_unlock(&mutex);
           
           // Update the current time
           now = time(NULL);
	   temp_StartTime = mySystem.myConfiguration.mySchedule.GetStartTime();
         }
         mySystem.myConfiguration.myLogger.WriteSystemLogString("Sequence finished");

         // Close the shutter
         mySystem.myConfiguration.myFilterWheel->CloseShutter();

         // Turn the temperature regulation off if we are auto regulating
       	 if(mySystem.myConfiguration.myCamera->GetAutoTemp())
       			mySystem.myConfiguration.myCamera->SetTemperatureRegulationOff();

         // If we are writing PNGs, go back and rewrite the pngs with better contrast

         // If we are creating MPEGs, do so
         
       } else if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_OFF)
       {
         // The schedule mode was switched to off
       } else if(!mySystem.myIsRunning)
       {
         // The system has been turned off, exit the thread
       }
          
       // The current schedule was finished.  If the mode was manual, switch it to off
       if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_MANUAL)
       {
         mySystem.myConfiguration.mySchedule.SetScheduleMode(SCHEDULE_OFF);
       }
    }
  }

  mySystem.myConfiguration.myLogger.WriteSystemLogString("Stopping control thread");
  mySystem.myControlRunning = false;

  pthread_exit(NULL);
}
