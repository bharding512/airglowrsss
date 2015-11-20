/***************************************************************************
                          client.cpp  -  description
                             -------------------
    begin                : Fri Aug 15 2003
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

#include "client.h"

#include "SocketException.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include "ServerCommands.h"
#include "schedule.h"

Client::Client(){
  myServerAddress = "";
  myServerPort = 0;
}
Client::~Client(){
}
/** Tries to connect to the server */
void Client::Connect(){
  try
  {
    cout << myServerAddress << " " << myServerPort << endl;
    ClientSocket myClient(myServerAddress, myServerPort);
  }
  catch(SocketException& e)
  {
    cout << "Exception was caught: " << e.description() << endl;
  }

  cout << "Connected" << endl;
}
/** No descriptions */
bool Client::SendCommand(string cmd){

  cout << cmd << endl;

  ClientSocket myClient(myServerAddress, myServerPort);
  
  try
  {
    myClient << cmd;
    myClient >> myLastReply;
  }
  catch(SocketException &e)
  {
    cout << e.description() << endl;
    return false;
  }

  return true;
}
