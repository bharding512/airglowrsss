/***************************************************************************
                          user.h  -  description
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

#ifndef USER_H
#define USER_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

/**Class to hold the user's information.
  *@author Jonathan Makela
  */

class User {
public: 
	User();
  User(string name, string email, bool daily, bool error);
	~User();
  /** Set's the user's email address */
  void SetUserEmail(string email);
  /** Sets the user's name */
  void SetUserName(string name);
  /** Sets the error email flag */
  void SetEmailOnError(bool flag);
  /** Sets the daily email flag. */
  void SetEmailDaily(bool flag);
  /** Gets the email on error flag */
  bool GetEmailOnError();
  /** Gets the email daily flag */
  bool GetEmailDaily();
  /** Gets the user's name */
  string GetUserName();
  /** Gets the user's email address */
  string GetUserEmail();
  /** Copies the user class */
  void Copy(User user);
  /** Copy constructor */
   User(const User& u);
  /** Write the class structure */
  void Write(ofstream& os);
  /** Read in the class structure */
  void Read(ifstream& is);
  /** No descriptions */
  void WriteXML(ofstream &os);
protected: // Public attributes
  /**  */
  string myUserName;
  string myUserEmail;
  bool myEmailOnError;
  bool myEmailDaily;
};

#endif
