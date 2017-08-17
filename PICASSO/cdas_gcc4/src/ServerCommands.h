/***************************************************************************
                          ServerCommands.h  -  description
                             -------------------
    begin                : Mon Aug 11 2003
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

 // Defines the commands that can be sent to the server

 #define SERVER_STOP    0   // Stops the system from running  (O:)
 #define SCHEDULE_MODE  1   // Sets the schedule mode (1:MODE)
                            // MODE in schedule.h
                            // If SCHEDULE_MANUAL, needs to include start and stop time
                            // (time_t) e.g., 1:SCHEDULE_MANUAL:start:stop
 #define ANGLES_SET     2   // Sets the horizon information (2:ANGLE)
                            // ANGLE in schedule.h
                            // If ANGLES_MOON or ANGLES_SUN needs to be followed by
                            // the angle in degrees
 #define SET_SITENAME   3   // Sets the site name
 #define SET_SITEABBR   4   // Sets the site abbreviation
 #define SET_LATITUDE   5   // Sets the latitude
 #define SET_LONGITUDE  6   // Sets the longitude
 #define SET_ALTITUDE   7   // Sets the altitude
 #define SET_XBIN       8   // Sets the xbin factor
 #define SET_YBIN       9   // Sets the ybin factor
 #define SET_TOP        10  // Sets the top corner of the imaged CCD
 #define SET_LEFT       11  // Sets the left corner of the imaged CCD
 #define SET_HEIGHT     12  // Sets the height of the imaged CCD
 #define SET_WIDTH      13  // Sets the width of the imaged CCD
 #define SET_TEMP       14  // Sets the temperature setpoint for the CCD
 #define SET_CAMERATYPE 15  // Sets the camera type
 #define CAPTURE_IMAGE  16  // Captures an image
 #define SET_DATA_DIR   17  // Sets the data path for saving *tif images
 #define SET_PNG_DIR    18  // Sets the path for saving *png images
 #define SET_MPG_DIR    19  // Sets the path for saving *mpg movies
 #define SET_QL_DIR     20  // Sets the path for saving quicklook images
 #define SET_DO_PNG     21  // Sets the variable on whether or not to save *png images
 #define SET_DO_MPG     22  // Sets the variable on whether or not to save *mpg movies
 #define SET_DO_QL      23  // Sets the variable on whether or not to save quicklook images
 #define SET_FILTERWHEELTYPE 24 // Sets the filterwheel type
 #define SET_FILTERINFO 25  // Sets the variables for each individual filter
 #define SET_SUBSEQUENCE 26 // Sets the sub-sequence for the indicated sequence
 #define SET_SUBSEQUENCEORDER 27  // Sets the sub-sequence order for the master sequence
 #define SET_DARKSEQUENCE 28      // Sets the variable on whether to take dark images in the sequence
 #define SET_PERIODICITY  29      // Sets the variable on whether to wait in between images in the sequence
 #define SET_GAIN         30      // Sets the CCD gain
 #define SET_CCD_AUTO_TEMP 31		  // Sets the CCD auto temperature regulation parameter
 #define SET_SCHEDULELATITUDE 32	// Sets the latitude used in the schedule calculations
 #define SET_SCHEDULELONGITUDE 33	// Sets the longitude used in the schedule calculations
 #define SET_NUMFILTERS 34  // Sets the number of filter positions in the filter wheel

 // Commands the will return a value
 #define GET_SCHEDULE_MODE  1001  // Gets the schedule mode
 #define GET_SUNANGLE       1002  // Gets the sun angle
 #define GET_MOONANGLE      1003  // Get the moon angle
 #define GET_DOMOONSET      1004  // Gets the moonset flag
 #define GET_STARTTIME      1005  // Gets the current start time
 #define GET_STOPTIME       1006  // Gets the current stop time
 #define GET_SITENAME       1007  // Gets the site name
 #define GET_SITEABBR       1008  // Gets the site abbreviation
 #define GET_LATITUDE       1009  // Gets the latitude
 #define GET_LONGITUDE      1010  // Gets the longitude
 #define GET_ALTITUDE       1011  // Gets the altitude
 #define GET_XBIN           1012  // Gets the xbin factor
 #define GET_YBIN           1013  // Gets the ybin factor
 #define GET_TOP            1014  // Gets the top corner of the CCD
 #define GET_LEFT           1015  // Gets the left corner of the CCD
 #define GET_HEIGHT         1016  // Gets the height of the imaged CCD
 #define GET_WIDTH          1017  // Gets the width of the imaged CCD
 #define GET_TEMP           1018  // Gets the temperature setpoint for the CCD
 #define GET_CAMERATYPE     1019  // Gets the camera type
 #define GET_DATA_DIR       1020  // Gets the data directory
 #define GET_PNG_DIR        1021  // Gets the png directory
 #define GET_MPG_DIR        1022  // Gets the mpg directory
 #define GET_QL_DIR         1023  // Gets the quicklook directory
 #define GET_DO_PNG         1024  // Gets the variable on whether or not to save png image
 #define GET_DO_MPG         1025  // Gets the variable on whether or not to save mpg movies
 #define GET_DO_QL          1026  // Gets the variable on whether or not to save quicklook images
 #define GET_FILTERWHEELTYPE 1027 // Gets the filterwheel type
 #define GET_FILTERINFO     1028  // Gets the variables pertaining to the filters
 #define GET_SUBSEQUENCE    1029  // Gets the sub-sequence indicated
 #define GET_SUBSEQUENCEORDER 1030  // Gets the sub-sequence order
 #define GET_DARKSEQUENCE   1031  // Gets the variable on whether or not to do dark images in sequence
 #define GET_PERIODICITY    1032  // Gets the variable on whether or not to wait in between images in the sequence
 #define GET_GAIN           1033  // Gets the CCD gain
 #define GET_CCD_AUTO_TEMP  1034  // Gets the CCD auto temperature regulation parameter
 #define GET_ACTUAL_TEMP    1035  // Read actual CCD temperature
 #define GET_SCHEDULELATITUDE 1036	// Gets the latitude used in the scheduler
 #define GET_SCHEDULELONGITUDE 1037	// Gets the longitude used in the scheduler
 #define GET_NUMFILTERS 1038      // Get the number of filter positions in the filter wheel
