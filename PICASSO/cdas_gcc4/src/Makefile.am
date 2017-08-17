bin_PROGRAMS = cdas
cdas_SOURCES = ClientSocket.cpp ServerSocket.cpp Socket.cpp SunAndMoon.cpp \
	calibration.cpp camera.cpp camerafli.cpp camerasbiguniv.cpp \
	configuration.cpp directories.cpp filterwheel.cpp filterwheelcfwl.cpp filterwheelfli.cpp \
	filterwheelkeo.cpp image.cpp logger.cpp main.cpp schedule.cpp sequencer.cpp site.cpp \
	cameraandor.cpp camerapvcam.cpp user.cpp system.cpp filterwheelserial.cpp

cdas_LDADD = -lpthread -lpng -lfli -lsbigudrv -lm -landor -ltcl8.5 -lpvcam -lraw1394

# set the include path found by configure
INCLUDES= $(all_includes)

# the library search path.
cdas_LDFLAGS = $(all_libraries)

# Define the APOGEE_PPI
##USERDEFINES = -DAPOGEE_PPI -DLINUX
