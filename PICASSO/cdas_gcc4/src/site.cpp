/***************************************************************************
                          site.cpp  -  description
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

#include "site.h"
#include "system.h"

extern System mySystem;

Site::Site(){
  myHeight = myLatitude = myLongitude = 0.0;
  myName = "";
  strcpy(myAbbreviation, "XXX");
}
Site::Site(float lat, float lon, float ht, string name, char abb[4]){
  // Constructor to construct a full class
  myHeight = ht;  // in km

  // Make sure the lat/lon make sense
  while(lon < 0) {
    mySystem.myConfiguration.myLogger.WriteErrorString("Longitude out of range (Site::Site)");
    lon = lon + 360.;
  }
  while(lon > 360) {
    mySystem.myConfiguration.myLogger.WriteErrorString("Longitude out of range (Site::Site)");
    lon = lon - 360.;
  }
  if(lat < -90 || lat > 90){
    mySystem.myConfiguration.myLogger.WriteErrorString("Latitude out of range (Site::Site)");
    lat = -999.;
  }

  myLatitude = lat;
  myLongitude = lon;

  char str[64];
  sprintf(str, "Set Lat: %8.3f, Long: %8.3f, Height: %8.3f", myLatitude, myLongitude, myHeight);
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);

  myName = name;
  strcpy(myAbbreviation, abb);

  sprintf(str, "Set Name: %s, Abbreviation: %s", myName.c_str(), myAbbreviation);
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
}
Site::~Site(){
}
/** Sets the location of the site */
void Site::SetLocation(float lat, float lon, float ht){
  myHeight = ht;  // in km

  // Make sure the lat/lon make sense
  while(lon < 0) {
    mySystem.myConfiguration.myLogger.WriteErrorString("Longitude out of range (Site::SetLocation)");
    lon = lon + 360.;
  }
  while(lon > 360) {
    mySystem.myConfiguration.myLogger.WriteErrorString("Longitude out of range (Site::SetLocation)");
    lon = lon - 360.;
  }
  if(lat < -90 || lat > 90){
    mySystem.myConfiguration.myLogger.WriteErrorString("Latitude out of range (Site::SetLocation)");
    lat = -999.;
  }

  myLatitude = lat;
  myLongitude = lon;

  char str[64];
  sprintf(str, "Set Lat: %8.3f, Long: %8.3f, Height: %8.3f", myLatitude, myLongitude, myHeight);
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);

  if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
    mySystem.myConfiguration.mySchedule.AutoSetTimes();
}
/** Sets the name and abbreviation for the site */
void Site::SetNames(string name, char abb[4]){
  char str[64];

  myName = name;
  strcpy(myAbbreviation, abb);

  sprintf(str, "Set Name: %s, Abbreviation: %s", myName.c_str(), myAbbreviation);
  mySystem.myConfiguration.myLogger.WriteSystemLogString(str);
}
/** Returns the site's latitude */
float Site::GetLatitude(){
  return myLatitude;
}
/** Returns the site's longitude */
float Site::GetLongitude(){
  return myLongitude;
}
/** Returns the site's height */
float Site::GetHeight(){
  return myHeight;
}
/** Return's the site's name */
string Site::GetName(){
  return myName;
}
/** Returns the site's abbreviation. */
char* Site::GetAbbreviation(){
  return myAbbreviation;
}
/** Copy constructor */
Site::Site(const Site& s){
  strcpy(myAbbreviation, s.myAbbreviation);
  myHeight = s.myHeight;
  myLatitude = s.myLatitude;
  myLongitude = s.myLongitude;
  myName = s.myName;
}
/** Read in the class information */
void Site::Read(ifstream& is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#Site" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #Site found");
    return;
  }

  string str;
  
  // Read the site name
  getline(is, myName);
  getline(is, str);
  strcpy(myAbbreviation, str.c_str());

  // Read the data
  is >> myLatitude;
  is >> myLongitude;
  is >> myHeight;
}
/** Write out the class structure */
void Site::Write(ofstream& os){
  // Put the header label
  os << "#Site" << endl;

  // Write the site name
  os << myName << endl;
  os << myAbbreviation << endl;

  // Write the Data
  os << myLatitude << endl;
  os << myLongitude << endl;
  os << myHeight << endl;
}
/** No descriptions */
void Site::SetSiteName(string name){
  myName = name;
}
/** No descriptions */
void Site::SetAbbreviation(char abbr[4]){
  strcpy(myAbbreviation, abbr);
}
/** No descriptions */
void Site::SetLatitude(float lat){
  if(lat < -90 || lat > 90){
    mySystem.myConfiguration.myLogger.WriteErrorString("Latitude out of range (Site::SetLocation)");
    lat = -999.;
  }

  myLatitude = lat;

//  if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
//    mySystem.myConfiguration.mySchedule.AutoSetTimes();
}
/** No descriptions */
void Site::SetLongitude(float lon){
  // Make sure the lat/lon make sense
  while(lon < 0) {
    mySystem.myConfiguration.myLogger.WriteErrorString("Longitude out of range (Site::SetLocation)");
    lon = lon + 360.;
  }
  while(lon > 360) {
    mySystem.myConfiguration.myLogger.WriteErrorString("Longitude out of range (Site::SetLocation)");
    lon = lon - 360.;
  }

  myLongitude = lon;

//  if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
//    mySystem.myConfiguration.mySchedule.AutoSetTimes();
}
/** No descriptions */
void Site::SetHeight(float ht){
  myHeight = ht;

//  if(mySystem.myConfiguration.mySchedule.GetScheduleMode() == SCHEDULE_AUTO)
//    mySystem.myConfiguration.mySchedule.AutoSetTimes();
}
/** No descriptions */
void Site::WriteXML(ofstream &os){
  // Put the header label
  os << "<site>" << endl;

  // Write the site name
  os << "<sitename>" << myName << "</sitename>" << endl;
  os << "<siteabbreviation> " << myAbbreviation << "</siteabbreviation>" << endl;

  // Write the Data
  os << "<latitude>" << myLatitude << "</latitude>" << endl;
  os << "<longitude>" << myLongitude << "</longitude>" << endl;
  os << "<height>" << myHeight << "</height>" << endl;

  os << "</site>" << endl;
}
