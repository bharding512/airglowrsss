/***************************************************************************
                          client.h  -  description
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

#ifndef CLIENT_H
#define CLIENT_H

#include "ClientSocket.h"

/**Handles the communication with the remote server
  *@author Jonathan Makela
  */

using namespace std;

class Client {
public: 
	Client();
	~Client();
  /** Tries to connect to the server */
  void Connect();
  /** No descriptions */
  bool SendCommand(string cmd);
protected: // Protected attributes
  /**  */
  ClientSocket myClient;
public: // Public attributes
  /**  */
  string myLastReply;
  /**  */
  int myServerPort;
  /**  */
  string myServerAddress;
};

#endif
