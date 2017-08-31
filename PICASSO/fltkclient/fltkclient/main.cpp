/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Fri Aug 15 13:07:41 EDT 2003
    copyright            : (C) |YEAR| by Jo2003 Makela
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>

#include <FL/Fl.H>
#include "fltkclientUI.h"

int main(int argc, char *argv[])
{
  fltkClientUI *myGUI = new fltkClientUI;

  Fl::visual(FL_SINGLE|FL_INDEX);

  myGUI->show(argc, argv);

  return Fl::run();

  return EXIT_SUCCESS;
}
