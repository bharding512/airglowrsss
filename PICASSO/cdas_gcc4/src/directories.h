/***************************************************************************
                          directories.h  -  description
                             -------------------
    begin                : Thu Aug 7 2003
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

#ifndef DIRECTORIES_H
#define DIRECTORIES_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

/**
  *@author Jonathan Makela
  */

class Directories {
public: 
	Directories();
	~Directories();
  /** Sets the path where the raw data should be stored */
  void SetDataPath(string path);
  /** Sets the path and flag for doing PNG images */
  void SetPNG(string path, bool flag);
  /** Sets the quick look directory and flag */
  void SetQuickLook(string path, bool flag);
  /** Sets the movies path and flag */
  void SetMovies(string path, bool flag);
  /** Copies the directory information */
  void Copy(Directories dir);
  /** Returns the Quicklook flag */
  bool GetDoQuickLook();
  /** Gets the DoMovies flag */
  bool GetDoMovies();
  /** Returns the PNG flag */
  bool GetDoPNG();
  /** Returns the Quicklook path */
  string GetQuickLookPath();
  /** Returns the movie path */
  string GetMoviesPath();
  /** Returns the PNG path */
  string GetPNGPath();
  /** Returns the current datapath */
  string GetDataPath();
  /** Copy constructor */
   Directories(const Directories& d);
  /** Read the information into the class */
  void Read(ifstream &is);
  /** Write the information in the class to a file */
  void Write(ofstream &os);
  /** No descriptions */
  void SetDoMovies(bool flag);
  /** No descriptions */
  void SetMoviesPath(string path);
  /** No descriptions */
  void SetDoPNG(bool flag);
  /** No descriptions */
  void SetPNGPath(string path);
  /** No descriptions */
  void SetDoQuickLook(bool flag);
  /** No descriptions */
  void SetQuickLookPath(string path);
  /** No descriptions */
  void WriteXML(ofstream &os);
  /** Returns the directory stub for this day.  The directories are in the format of YYYY/DDD where
YYYY is the year and DDD is the day number (0-365/366). */
  string GetTodaysDir();
  /** Sets the directory stub for today's directory */
  void SetTodaysDir();
//  tm* gmtime(&now arg1);
//  void time(&now arg1);
protected: // Public attributes
  /**  */
  string myDataPath;
  /**  */
  bool myDoPNG;
  /**  */
  bool myDoMovies;
  /**  */
  bool myDoQuickLook;
  /**  */
  string myPNGPath;
  /**  */
  string myMoviesPath;
  /**  */
  string myQuickLookPath;
  /**  */
  string myTodaysDir;
};

#endif
