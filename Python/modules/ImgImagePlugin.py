# See http://www.pythonware.com/library/pil/handbook/decoder.htm

try:
    import Image, ImageFile
except ImportError:
    from PIL import Image, ImageFile
import string
import struct
import datetime

class ImgImageFile(ImageFile.ImageFile):

    format = "IMG"
    format_description = "Sherwood img format for FPI images"

    def _open(self):

        # check header
	FileIdentifier = self.fp.read(4)
        if FileIdentifier != "A3OI":
            raise SyntaxError, "not an IMG file"

	# read in the header structure
	vn = struct.unpack('i',self.fp.read(4))
    	pgmVersion = struct.unpack('i',self.fp.read(4))

	# Total size of the header, followed by offsets/sizes for:
    	#   AndorVersion
    	#   DetectorInformation
    	#   CameraPameters
	#   Conditions
    	#   ImageInfo
    	hdrSize = struct.unpack('i',self.fp.read(4))
    	avOffset = struct.unpack('h',self.fp.read(2))
    	avSize = struct.unpack('h',self.fp.read(2))
    	diOffset = struct.unpack('h',self.fp.read(2))
    	diSize = struct.unpack('h',self.fp.read(2))
    	cpOffset = struct.unpack('h',self.fp.read(2))
    	cpSize = struct.unpack('h',self.fp.read(2))
    	cnOffset = struct.unpack('h',self.fp.read(2))
    	cnSize = struct.unpack('h',self.fp.read(2))
    	iiOffset = struct.unpack('h',self.fp.read(2))
    	iiSize = struct.unpack('h',self.fp.read(2))

	# Read in information on AndorVersion
    	self.fp.seek(avOffset[0])
    	vn = struct.unpack('i',self.fp.read(4))
    	PCBVer = struct.unpack('I',self.fp.read(4))
    	FlexVer = struct.unpack('I',self.fp.read(4))
    	DummyVer = struct.unpack('IIII',self.fp.read(16))
    	EpromVer = struct.unpack('I',self.fp.read(4))
    	CofFileVer = struct.unpack('I',self.fp.read(4))
    	VxdRev = struct.unpack('I',self.fp.read(4))
    	CxdVer = struct.unpack('I',self.fp.read(4))
    	DIIRev = struct.unpack('I',self.fp.read(4))
    	DIIVer = struct.unpack('I',self.fp.read(4))

	# Read in the Detector Information
    	self.fp.seek(diOffset[0]);
    	vn = struct.unpack('i',self.fp.read(4));
    	rows = struct.unpack('i',self.fp.read(4));
    	cols = struct.unpack('i',self.fp.read(4));

	# Read in the CameraParameters
    	self.fp.seek(cpOffset[0]);
    	vn = struct.unpack('i',self.fp.read(4));
    	ReadMode = struct.unpack('i',self.fp.read(4));
    	AcqMode = struct.unpack('i',self.fp.read(4));
    	TrigMode = struct.unpack('i',self.fp.read(4));
    	ShiftSpeedH = struct.unpack('i',self.fp.read(4));
    	ShiftSpeedV = struct.unpack('i',self.fp.read(4));
    	ShutterMode = struct.unpack('i',self.fp.read(4));
    	ShutterType = struct.unpack('i',self.fp.read(4));
    	ShutterOpenT = struct.unpack('i',self.fp.read(4));
    	ShutterCloseT = struct.unpack('i',self.fp.read(4));
    	FinalTemp = struct.unpack('i',self.fp.read(4));
    	AccumCycleT = struct.unpack('i',self.fp.read(4));
    	DDGGateStep = struct.unpack('i',self.fp.read(4));
    	ExposureTimeRq = struct.unpack('f',self.fp.read(4));
    	ExposureTime = struct.unpack('f',self.fp.read(4));
    	FKHeight = struct.unpack('i',self.fp.read(4));
    	FKNumber = struct.unpack('i',self.fp.read(4));
    	FKExposureT = struct.unpack('f',self.fp.read(4));
    	FKBinMode = struct.unpack('i',self.fp.read(4));
    	KineticCycleT = struct.unpack('i',self.fp.read(4));
    	FKVShiftSpeed = struct.unpack('i',self.fp.read(4));
    	XBinning = struct.unpack('i',self.fp.read(4)); 
    	YBinning = struct.unpack('i',self.fp.read(4));
    	FullImageCol0 = struct.unpack('i',self.fp.read(4));
    	FullImageCol1 = struct.unpack('i',self.fp.read(4));
    	FullImageRow0 = struct.unpack('i',self.fp.read(4));
    	FullImageRow1 = struct.unpack('i',self.fp.read(4));
    	MultiTrackNum = struct.unpack('i',self.fp.read(4));
    	MultiTrackHeight = struct.unpack('i',self.fp.read(4));
    	MultiTrackOffset = struct.unpack('i',self.fp.read(4));
    	MultiTrackBottom = struct.unpack('i',self.fp.read(4));
    	MultiTrackGap = struct.unpack('i',self.fp.read(4));
    	AccumScanNum = struct.unpack('i',self.fp.read(4));
    	KineticScanNum = struct.unpack('i',self.fp.read(4));
    	RandomTrackNum = struct.unpack('i',self.fp.read(4));
    	SingleTrackCenter = struct.unpack('i',self.fp.read(4));
    	SingleTrackHeight = struct.unpack('i',self.fp.read(4));

	# Read in Conditions
    	self.fp.seek(cnOffset[0])
    	vn = struct.unpack('i',self.fp.read(4))
    	ObsVn = struct.unpack('i',self.fp.read(4))
    	digitalInputs = struct.unpack('d',self.fp.read(8))
    	azAngle = struct.unpack('d',self.fp.read(8))
    	zeAngle = struct.unpack('d',self.fp.read(8))
    	OutsideTemperature = struct.unpack('d',self.fp.read(8))
    	Pressure = struct.unpack('d',self.fp.read(8))
        aux = struct.unpack('dddddddddddd',self.fp.read(8*12))

	# Read in Image Information
    	self.fp.seek(iiOffset[0])
    	vn = struct.unpack('i',self.fp.read(4))
    	tmp = struct.unpack('hhhhhhhh',self.fp.read(16))
    	LocalTime = datetime.datetime(tmp[0],tmp[1],tmp[3],tmp[4],tmp[5],tmp[6])
    	tmp = struct.unpack('hhhhhhhh',self.fp.read(16))
    	getCopyTime = datetime.datetime(tmp[0],tmp[1],tmp[3],tmp[4],tmp[5],tmp[6])
    	CCDTemperature = struct.unpack('i',self.fp.read(4));
        imageTemperatureStatus = struct.unpack('i',self.fp.read(4))
        sz = struct.unpack('i',self.fp.read(4))

        # Save pertinent header info
        self.info = {'ExposureTime':ExposureTime[0],'XBinning':XBinning[0],
                     'YBinning':YBinning[0],'azAngle':azAngle[0],'zeAngle':zeAngle[0],
                     'OutsideTemperature':OutsideTemperature[0],'Pressure':Pressure[0],
                     'LocalTime':LocalTime,'CCDTemperature':CCDTemperature[0]}

        # size in pixels (width, height)
        self._size = rows[0]/YBinning[0], cols[0]/XBinning[0]

	# mode setting
	self.mode = "I;16L"

        # data descriptor
        self.tile = [
            ("raw", (0, 0) + self.size, 1024, (self.mode, 0, 1))
        ]

Image.register_open("IMG", ImgImageFile)

Image.register_extension("IMG", ".img")
