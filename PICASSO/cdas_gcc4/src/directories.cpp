/***************************************************************************
                          directories.cpp  -  description
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

#include "directories.h"
#include "system.h"
#include <sys/types.h>
#include <sys/stat.h>

extern System mySystem;

Directories::Directories(){
  // initialize the variables
  myDataPath = "";
  myPNGPath = "";
  myMoviesPath = "";
  myQuickLookPath = "";
  myDoPNG = false;
  myDoMovies = false;
  myDoQuickLook = false;
}
Directories::~Directories(){
}
/** Sets the path where the raw data should be stored */
void Directories::SetDataPath(string path){
  string str("Changing data path to: ");
  myDataPath = path;
  str = str + path;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
}

/** Sets the path and flag for doing PNG images */
void Directories::SetPNG(string path, bool flag){
  string str("Changing the PNG path to: ");
  myPNGPath = path;
  myDoPNG = flag;
  string state;
  if(flag)
  {
    state = " (TRUE)";
  }
  else
  {
    state = " (FALSE)";
  }
  str = str + path + state;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);  
}
/** Sets the quick look directory and flag */
void Directories::SetQuickLook(string path, bool flag){
  string str("Changing the quicklook path to: ");
  myQuickLookPath = path;
  myDoQuickLook = flag;
  string state;
  if(flag)
  {
    state = " (TRUE)";
  }
  else
  {
    state = " (FALSE)";
  }
  str = str + path + state;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
}
/** Sets the movies path and flag */
void Directories::SetMovies(string path, bool flag){
  string str("Changing the movie path to: ");
  myMoviesPath = path;
  myDoMovies = flag;
  string state;
  if(flag)
  {
    state = " (TRUE)";
  }
  else
  {
    state = " (FALSE)";
  }
  str = str + path + state;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
}
/** Copies the directory information */
void Directories::Copy(Directories dir){
  myDataPath = dir.myDataPath;
  myDoPNG = dir.myDoPNG;
  myDoMovies = dir.myDoMovies;
  myDoQuickLook = dir.myDoQuickLook;
  myPNGPath = dir.myPNGPath;
  myMoviesPath = dir.myMoviesPath;
  myQuickLookPath = dir.myQuickLookPath;
}
/** Returns the current datapath */
string Directories::GetDataPath(){
  return myDataPath;
}
/** Returns the PNG path */
string Directories::GetPNGPath(){
  return myPNGPath;
}
/** Returns the movie path */
string Directories::GetMoviesPath(){
  return myMoviesPath;
}
/** Returns the Quicklook path */
string Directories::GetQuickLookPath(){
  return myQuickLookPath;
}
/** Returns the PNG flag */
bool Directories::GetDoPNG(){
  return myDoPNG;
}
/** Gets the DoMovies flag */
bool Directories::GetDoMovies(){
  return myDoMovies;
}
/** Returns the Quicklook flag */
bool Directories::GetDoQuickLook(){
  return myDoQuickLook;
}
/** Copy constructor */
 Directories::Directories(const Directories& d){
   myDataPath = d.myDataPath;
   myDoPNG = d.myDoPNG;
   myDoMovies = d.myDoMovies;
   myDoQuickLook = d.myDoQuickLook;
   myPNGPath = d.myPNGPath;
   myMoviesPath = d.myMoviesPath;
   myQuickLookPath = d.myQuickLookPath;
}
/** Write the information in the class to a file */
void Directories::Write(ofstream &os){
  // Put the header label
  os << "#Directories" << endl;

  // Write the paths
  os << myDataPath << endl;
  os << myPNGPath << endl;
  os << myMoviesPath << endl;
  os << myQuickLookPath << endl;

  // Write the flags
  os << myDoPNG << endl;
  os << myDoMovies << endl;
  os << myDoQuickLook << endl;
}
/** Read the information into the class */
void Directories::Read(ifstream &is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#Directories" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #Directories found");
    return;
  }
  
  // Read the paths
  getline(is, myDataPath);
  getline(is, myPNGPath);
  getline(is, myMoviesPath);
  getline(is, myQuickLookPath);

  // Read the flags
  is >> myDoPNG;
  is >> myDoMovies;
  is >> myDoQuickLook;
}
/** No descriptions */
void Directories::SetPNGPath(string path){
  SetPNG(path, myDoPNG);
}
/** No descriptions */
void Directories::SetDoPNG(bool flag){
  SetPNG(myPNGPath, flag);
}
/** No descriptions */
void Directories::SetMoviesPath(string path){
  SetMovies(path, myDoMovies);
}
/** No descriptions */
void Directories::SetDoMovies(bool flag){
  SetMovies(myMoviesPath, flag);
}
/** No descriptions */
void Directories::SetQuickLookPath(string path){
  SetQuickLook(path, myDoQuickLook);
}
/** No descriptions */
void Directories::SetDoQuickLook(bool flag){
  SetQuickLook(myQuickLookPath, flag);
}
/** No descriptions */
void Directories::WriteXML(ofstream &os){
  // Put the header label
  os << "<directories>" << endl;

  // Write the paths
  os << "<datapath>" << myDataPath << "</datapath>" << endl;
  os << "<pngpath>" << myPNGPath << "</pngpath>" << endl;
  os << "<moviespath>" << myMoviesPath << "</moviespath>" << endl;
  os << "<quicklookpath>" << myQuickLookPath << "</quicklookpath>" << endl;

  // Write the flags
  os << "<dopng>" << myDoPNG << "</dopng>" << endl;
  os << "<domovies> " << myDoMovies << "</domovies>" << endl;
  os << "<doquicklook>" << myDoQuickLook << "</doquicklook>" << endl;

  os << "</directories>" << endl;
}
/** Returns the directory stub for this day.  The directories are in the format of YYYY/DDD where
YYYY is the year and DDD is the day number (0-365/366). */
string Directories::GetTodaysDir(){
  return myTodaysDir;
}
/** Sets the directory stub for today's directory */
void Directories::SetTodaysDir(){
  time_t now;
  struct tm* temp_tm;

  // Get the current time
  time(&now);
  temp_tm = gmtime(&now);

  // Grab the year and day of year
  int year = temp_tm->tm_year + 1900;
  int doy = temp_tm->tm_yday + 1;

  // Create the YYYY/DDD
  char myChar[12];
  sprintf(myChar, "%.4d/%.3d/", year, doy);

  // Save it to the myTodaysDir variable
  string myDir(myChar);
  myTodaysDir = myDir;

  string temp;

  // This is the default mode for creating our directories (RWXR-XR-X)
  mode_t myMode;
  myMode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  
  // Go through and create these directories
  if(myDataPath != "")
  {
    mkdir(myDataPath.c_str(),myMode);
    sprintf(myChar,"%.4d", year);
    temp = myDataPath + string(myChar);
    mkdir(temp.c_str(),myMode);
    sprintf(myChar,"/%.3d", doy);
    temp = temp + string(myChar);
    mkdir(temp.c_str(),myMode);
  }

  if(myMoviesPath != "")
  {
    mkdir(myMoviesPath.c_str(),myMode);
    sprintf(myChar,"%.4d", year);
    temp = myMoviesPath + string(myChar);
    mkdir(temp.c_str(),myMode);
    sprintf(myChar,"/%.3d", doy);
    temp = temp + string(myChar);
    mkdir(temp.c_str(),myMode);
  }

  if(myPNGPath != "")
  {
    mkdir(myPNGPath.c_str(),myMode);
    sprintf(myChar,"%.4d", year);
    temp = myPNGPath + string(myChar);
    mkdir(temp.c_str(),myMode);
    sprintf(myChar,"/%.3d", doy);
    temp = temp + string(myChar);
    mkdir(temp.c_str(),myMode);
  }
}
