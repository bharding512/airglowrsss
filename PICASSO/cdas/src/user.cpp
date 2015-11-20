/***************************************************************************
                          user.cpp  -  description
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

#include "user.h"
#include "system.h"

extern System mySystem;

User::User(){
  myEmailDaily = myEmailOnError = false;
  myUserEmail = "";
  myUserName = "";
}
User::User(string name, string email, bool daily, bool error){
  myEmailDaily = daily;
  myEmailOnError = error;
  myUserName = name;
  myUserEmail = email;
}  
User::~User(){
}
/** Set's the user's email address */
void User::SetUserEmail(string email){
  myUserEmail = email;
}
/** Sets the user's name */
void User::SetUserName(string name){
  myUserName = name;
}
/** Sets the daily email flag. */
void User::SetEmailDaily(bool flag){
  myEmailDaily = flag;
}
/** Sets the error email flag */
void User::SetEmailOnError(bool flag){
  myEmailOnError = flag;
}
/** Gets the user's email address */
string User::GetUserEmail(){
  return myUserEmail;
}
/** Gets the user's name */
string User::GetUserName(){
  return myUserName;
}
/** Gets the email daily flag */
bool User::GetEmailDaily(){
  return myEmailDaily;
}
/** Gets the email on error flag */
bool User::GetEmailOnError(){
  return myEmailOnError;
}
/** Copies the user class */
void User::Copy(User user){
  myEmailDaily = user.GetEmailDaily();
  myEmailOnError = user.GetEmailOnError();
  myUserEmail = user.GetUserEmail();
  myUserName = user.GetUserName();
}
/** Copy constructor */
 User::User(const User& u){
   myEmailDaily = u.myEmailDaily;
   myEmailOnError = u.myEmailOnError;
   myUserEmail = u.myUserEmail;
   myUserName = u.myUserName;
}
/** Read in the class structure */
void User::Read(ifstream& is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#User" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #User found");
    return;
  }

  // Read in the strings
  getline(is, myUserEmail);
  getline(is, myUserName);

  // Read in the flags
  is >> myEmailDaily;
  is >> myEmailOnError;
}
/** Write the class structure */
void User::Write(ofstream& os){
  // Put the header label
  os << "#User" << endl;

  // Write the string
  os << myUserEmail << endl;
  os << myUserName << endl;
  
  // Write the flags
  os << myEmailDaily << endl;
  os << myEmailOnError << endl;
}
/** No descriptions */
void User::WriteXML(ofstream &os){
  // Put the header label
  os << "<user>" << endl;

  // Write the string
  os << "<useremail>" << myUserEmail << "</useremail>" << endl;
  os << "<username>" << myUserName << "</username>" << endl;

  // Write the flags
  os << "<useremaildaily>" << myEmailDaily << "</useremaildaily>" << endl;
  os << "<useremailonerror>" << myEmailOnError << "</useremailonerror>" << endl;

  os << "</user>" << endl;
}
