/***************************************************************************
                          logger.h  -  description
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

#ifndef LOGGER_H
#define LOGGER_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

/**The class to handle error and other logs
  *@author Jonathan Makela
  */

class Logger {
public: 
	Logger();
	~Logger();
  /** Opens the error file */
  bool OpenErrorFile(string FileName);
  /** Append a string to the error log */
  void WriteErrorString(string str);
  /** No descriptions */
  void WriteSystemLogString(string str);
  /** No descriptions */
  bool OpenSystemLogFile(string FileName);
  /** copy constructor */
   Logger(const Logger& l);
  /** Sets the Errorlog flag */
  void SetDoErrorLog(bool flag);
  /** Sets the system logger flag */
  void SetDoSystemLog(bool flag);
  /** Read the logger information to the config file */
  void Read(ifstream & is);
  /** Write the logger information to the config file */
  void Write(ofstream& os);
  /** returns system log filename */
  string GetSystemLogFileName();
  /** Returns the error filename */
  string GetErrorFileName();
  /** Gets the system log flag */
  bool GetDoSystemLog();
  /** Gets the Error log flag */
  bool GetDoErrorLog();
  /** No descriptions */
  void WriteXML(ofstream &os);
protected: // Public attributes
  /**  */
  string myErrorFileName;
  /**  */
  ofstream myErrorFile;
  /**  */
  bool myDoErrorLog;
  /**  */
  string mySystemLogFileName;
  /**  */
  ofstream mySystemLogFile;
  /**  */
  bool myDoSystemLog;
};

#endif
