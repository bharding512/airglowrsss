/***************************************************************************
                          schedule.h  -  description
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

#ifndef SCHEDULE_H
#define SCHEDULE_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "SunAndMoon.h"

using namespace std;

#define SCHEDULE_OFF    0   // The system has no schedule
#define SCHEDULE_AUTO   1   // The system should create its own schedule
#define SCHEDULE_MANUAL 2   // The user is specifying the schedule

#define ANGLE_MOON 0        // Set the moon angle
#define ANGLE_SUN 1         // Set the sun angle
#define ANGLE_DOMOON 2       // Set the myDoMoon flag

/**The class to handle the start and stop time for the system.
  *@author Jonathan Makela
  */

class Schedule {
public: 
	Schedule();
	~Schedule();
  /** Set the moon angle */
  void SetMoonAngle(double angle);
  /** Set the sun angle above horizon */
  void SetSunAngle(double angle);
  /** Copy constructor */
   Schedule(const Schedule& s);
  /** Read Schedule class from config file */
  void Read(ifstream& is);
  /** Write the schedule class to config file */
  void Write(ofstream& os);
  /** Returns the moonset flag */
  bool GetMoonSet();
  /** Set the moonset flag */
  void SetMoonSet(bool flag);
  /** Returns the start time */
  tm GetStartTime();
  /** Returns the sun angle variable */
  double GetSunAngle();
  /** Return the moonangle variable */
  double GetMoonAngle();
  /** Sets the stop time */
  void SetStopTime(tm stop);
  /** Gets the stop time */
  tm GetStopTime();
  /** Sets the start time */
  void SetStartTime(tm start);
  /** Returns whether the current time is between rise and set or not */
  int isUp(double time, double rise, double set);
  /** Sets the start and stop time based upon the current scheduling configuration and the sun and moon rise/set times. */
  bool AutoSetTimes();
  /** Sets the schedule mode */
  void SetScheduleMode(int mode);
  /** returns the schedulemode flag */
  int GetScheduleMode();
  /** Calculates the rise and set times of the sun and moon for the given day and location. */
  void GetRiseSetTimes(int month, int day, int year, double lat, double lon, int timezone, double *sunrise, double *sunset, double *moonrise, double *moonset);

protected: // Protected attributes
  /**  */
  tm myStartTime;
  /**  */
  tm myStopTime;
  /**  */
  double myMoonAngle;
  /**  */
  double mySunAngle;
  /**  */
  int myScheduleMode;
  /**  */
  bool myMoonSet;
  /**  */
  CSunAndMoon mySunAndMoon;
protected: // Protected methods
  /** Converts from degrees to radians */
  double DegRad(float Degrees);
};

#endif
