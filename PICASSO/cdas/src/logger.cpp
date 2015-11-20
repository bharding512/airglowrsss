/***************************************************************************
                          logger.cpp  -  description
                             -------------------
    begin                : Wed Aug 6 2003
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

#include "logger.h"
#include "system.h"

extern System mySystem;

Logger::Logger(){
  myErrorFileName = "";
  myDoErrorLog = false;
  mySystemLogFileName = "";
  myDoSystemLog = false;
  ofstream myErrorFile();
  ofstream mySystemLogFile();
}
Logger::~Logger(){
  if(myErrorFile)
  {
    WriteErrorString("Closing error file");
    myErrorFile.close();
  }
  if(mySystemLogFile)
  {
    WriteSystemLogString("Closing system log file");
    mySystemLogFile.close();
  }
}
/** Opens the error file */
bool Logger::OpenErrorFile(string FileName){
  myErrorFileName = FileName;

  if(myErrorFile.is_open())
    myErrorFile.close();

  // Open the Error file for output and appending
  myErrorFile.open(myErrorFileName.c_str(), ofstream::out | ofstream::app);
  
  if(!myErrorFile)
  {
    string str;
    str = "Cannot open the error file: " + myErrorFileName;
    WriteSystemLogString(str);
    myDoErrorLog = false;
    return false;
  }

  // We were able to open the error log so let's use it
  myDoErrorLog = true;

  string str;

  str = "Opened Error log file: " + myErrorFileName;
  WriteSystemLogString(str);

  WriteErrorString("Error log opened");
  
  return true;
}
/** Append a string to the error log */
void Logger::WriteErrorString(string str){
  time_t now;
  char time_string[64];

  // Write error to display
  cout << str << endl;

  if(myDoErrorLog)
  {
    // Get the current time
    time(&now);

    // Convert the time to a string (LT)
    strftime(time_string, 64, "%a %b %d %H:%M:%S %Y", localtime(&now));
    
    // Write the timestamp and error string to the Error file if it exists
    if(myErrorFile)
    {
      myErrorFile << "[" << string(time_string) << "]: " << str << endl;
    }
  }
}
/** No descriptions */
bool Logger::OpenSystemLogFile(string FileName){
  mySystemLogFileName = FileName;

  if(mySystemLogFile.is_open())
    mySystemLogFile.close();
  
  // Open the System file for output and appending
  mySystemLogFile.open(mySystemLogFileName.c_str(), ofstream::out | ofstream::app);
  
  if(!mySystemLogFile)
  {
    string str;
    str = "Cannot open the system log file: " + mySystemLogFileName;
    WriteSystemLogString(str);
    myDoSystemLog = false;
    return false;
  }

  // We were able to open the error log so let's use it
  myDoSystemLog = true;

  string str;

  str = "Opened System log file: " + FileName;
  WriteSystemLogString(str);

  return true;
}
/** No descriptions */
void Logger::WriteSystemLogString(string str){
  time_t now;
  char time_string[64];

  // Write string to the display
  cout << str << endl;

  if(myDoSystemLog)
  {
    // Get the current time
    time(&now);

    // Convert the time to a string (LT)
    strftime(time_string, 64, "%a %b %d %H:%M:%S %Y", localtime(&now));

    // Write the timestamp and error string to the system log file if it exists
    if(mySystemLogFile)
    {
      mySystemLogFile << "[" << string(time_string) << "]: " << str << endl;
    }
  }
}
/** copy constructor */
Logger::Logger(const Logger& l){
  myDoErrorLog = l.myDoErrorLog;
  myDoSystemLog = l.myDoSystemLog;
  myErrorFileName = l.myErrorFileName;
  mySystemLogFileName = l.mySystemLogFileName;

  // Now need to open the log files if needed
//  if(myDoErrorLog)
//    OpenErrorFile(myErrorFileName);

//  if(myDoSystemLog)
//    OpenSystemLogFile(mySystemLogFileName);
}
/** Sets the Errorlog flag */
void Logger::SetDoErrorLog(bool flag){
  myDoErrorLog = flag;

  string str("Setting the Error log flag to ");
  
  string state;
  if(flag)
  {
    state = "TRUE";
  }
  else
  {
    state = "FALSE";
  }
  str = str + state;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);  
}
/** Sets the system logger flag */
void Logger::SetDoSystemLog(bool flag){
  myDoSystemLog = flag;

  string str("Setting the System log flag to ");

  string state;
  if(flag)
  {
    state = "TRUE";
  }
  else
  {
    state = "FALSE";
  }
  str = str + state;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
}
/** Write the logger information to the config file */
void Logger::Write(ofstream& os){
  // Put the header label
  os << "#Logger" << endl;

  // Write the filenames
  os << myErrorFileName << endl;
  os << mySystemLogFileName << endl;

  // Write the flags
  os << myDoErrorLog << endl;
  os << myDoSystemLog << endl;
}
/** Read the logger information to the config file */
void Logger::Read(ifstream & is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#Logger" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    WriteErrorString("End of configuration file reached before #Logger found");
    return;
  }

  // Read the paths
  getline(is, myErrorFileName);
  getline(is, mySystemLogFileName);

  // Read the flags
  is >> myDoErrorLog;
  is >> myDoSystemLog;
  
  // Now need to open the log files if needed
  if(myDoSystemLog)
    OpenSystemLogFile(mySystemLogFileName);

  if(myDoErrorLog)
    OpenErrorFile(myErrorFileName);
}
/** Gets the Error log flag */
bool Logger::GetDoErrorLog(){
  return myDoErrorLog;
}
/** Gets the system log flag */
bool Logger::GetDoSystemLog(){
  return myDoSystemLog;
}
/** Returns the error filename */
string Logger::GetErrorFileName(){
  return myErrorFileName;
}
/** returns system log filename */
string Logger::GetSystemLogFileName(){
  return mySystemLogFileName;
}
/** No descriptions */
void Logger::WriteXML(ofstream &os){
  // Put the header label
  os << "<logger>" << endl;

  // Write the filenames
  os << "<errorfilename>" << myErrorFileName << "</errorfilename>" << endl;
  os << "<systemlogfilename>" << mySystemLogFileName << "</systemlogfilename>" << endl;

  // Write the flags
  os << "<doerrorlog>" << myDoErrorLog << "</doerrorlog>" << endl;
  os << "<dosystemlog>" << myDoSystemLog << "</dosystemlog>" << endl;

  os << "</logger>" << endl;
}
