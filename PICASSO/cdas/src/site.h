/***************************************************************************
                          site.h  -  description
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

#ifndef SITE_H
#define SITE_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

/**Class for the site information.
  *@author Jonathan Makela
  */

class Site {
public: 
	Site();
  Site(float lat, float lon, float ht, string name, char abb[4]);
	~Site();
  /** Sets the location of the site */
  void SetLocation(float lat, float lon, float ht);
  /** No descriptions */
  void SetNames(string name, char abb[4]);
  /** Return's the site's name */
  string GetName();
  /** Returns the site's height */
  float GetHeight();
  /** Returns the site's longitude */
  float GetLongitude();
  /** No descriptions */
  float GetLatitude();
  /** Returns the site's abbreviation. */
  char* GetAbbreviation();
  /** Copy constructor */
   Site(const Site& s);
  /** Write out the class structure */
  void Write(ofstream& os);
  /** Read in the class information */
  void Read(ifstream& is);
  /** No descriptions */
  void SetSiteName(string name);
  /** No descriptions */
  void SetHeight(float ht);
  /** No descriptions */
  void SetLongitude(float lon);
  /** No descriptions */
  void SetLatitude(float lat);
  /** No descriptions */
  void SetAbbreviation(char abbr[4]);
  /** No descriptions */
  void WriteXML(ofstream &os);
protected: // Public attributes
  /**  */
  float myHeight;
  /**  */
  float myLatitude;
  /**  */
  float myLongitude;
  /**  */
  string myName;
  /**  */
  char myAbbreviation[4];
};

#endif
