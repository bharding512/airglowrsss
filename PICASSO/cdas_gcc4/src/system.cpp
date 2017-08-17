/***************************************************************************
                          system.cpp  -  description
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

#include "system.h"

System::System(){
  myConfFileName = "CDAS.conf";
  myIsRunning = false;
  myServerRunning = false;
  myControlRunning = false;
  myIsExposing = false;
  myXMLFilename = "CDAS.xml";
}
System::~System(){
}
