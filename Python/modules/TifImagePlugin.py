# See http://www.pythonware.com/library/pil/handbook/decoder.htm

try:
    import Image, ImageFile
except ImportError:
    from PIL import Image, ImageFile
import string
import struct
import datetime

class TifImageFile(ImageFile.ImageFile):

    format = "TIF"
    format_description = "Garcia Image Format"

    def _open(self):

        # Skip the TIF header stuff
        self.fp.seek(200)

        # check header
        FileIdentifier = self.fp.read(3)
        if FileIdentifier != "CDF":
            raise SyntaxError("not a CDF TIF file")

        # Read in the number of rows and columns
        rows = struct.unpack('h', self.fp.read(2))
        cols = struct.unpack('h', self.fp.read(2))

        # Read in the autoscaled min and max pixel values
        pmin = struct.unpack('h', self.fp.read(2))
        pmax = struct.unpack('h', self.fp.read(2))

        # Read in the local and universal time structures
        tmp = struct.unpack('18h', self.fp.read(2*18))
        LocalTime = datetime.datetime(tmp[10]+1900,tmp[8]+1,tmp[6],tmp[4],tmp[2],tmp[0])
        tmp = struct.unpack('18h', self.fp.read(2*18))
        UniversalTime = datetime.datetime(tmp[10]+1900,tmp[8]+1,tmp[6],tmp[4],tmp[2],tmp[0])

        # Read the filter name
        tmp = struct.unpack('12s', self.fp.read(12))
        Filter = tmp[0].split('\x00')[0]

        # Read the emission height
        EmissionHeight = struct.unpack('f', self.fp.read(4))

        # Read the exposure time
        ExposureTime = struct.unpack('f', self.fp.read(4))

        # Read in CCD parameters
        gain = struct.unpack('B', self.fp.read(1))
        XBinning = struct.unpack('B', self.fp.read(1))
        YBinning = struct.unpack('B', self.fp.read(1))
        CCDTemperature = struct.unpack('f', self.fp.read(4))

        # Read in the FW temperature
        FWTemperature = struct.unpack('f', self.fp.read(4))

        # Read in XY spacing
        XSpacing = struct.unpack('f', self.fp.read(4))
        YSpacing = struct.unpack('f', self.fp.read(4))

        # Read in the spatial calibration coefficients
        a0 = struct.unpack('d', self.fp.read(8))
        a1 = struct.unpack('d', self.fp.read(8))
        a2 = struct.unpack('d', self.fp.read(8))
        b0 = struct.unpack('d', self.fp.read(8))
        b1 = struct.unpack('d', self.fp.read(8))
        b2 = struct.unpack('d', self.fp.read(8))
        Calibration = (a0[0], a1[0], a2[0], b0[0], b1[0], b2[0])

        # Read in site location
        lat = struct.unpack('f', self.fp.read(4))
        lon = struct.unpack('f', self.fp.read(4))
        alt = struct.unpack('f', self.fp.read(4))
        tmp = struct.unpack('30s', self.fp.read(30))
        Observatory = (tmp[0].split('\x00')[0], lat[0], lon[0], alt[0])

        # Read in comment
        tmp = struct.unpack('100s', self.fp.read(100))
        Comment = tmp[0].split('\x00')[0]

        # Read in additional information
        isAllSky = struct.unpack('B', self.fp.read(1))
        CenterAz = struct.unpack('d', self.fp.read(8))
        CenterEl = struct.unpack('d', self.fp.read(8))
        ProjectionAltitude =  struct.unpack('f', self.fp.read(4))
        ProjectionLon = struct.unpack('f', self.fp.read(4))
        tmp = struct.unpack('63B', self.fp.read(63))

        # Save pertinent header info
        self.info = {'LocalTime':LocalTime, 'UniversalTime': UniversalTime,
                     'pmin': pmin[0], 'pmax': pmax[0], 'Filter': Filter,
                     'EmissionHeight': EmissionHeight[0], 'ExposureTime': ExposureTime[0],
                     'gain': gain[0], 'XBinning': XBinning[0], 'YBinning': YBinning[0],
                     'CCDTemperature': CCDTemperature[0], 'FWTemperature': FWTemperature[0],
                     'Calibration': Calibration, 'Observatory': Observatory,
                     'Comment': Comment, 'isAllSky': isAllSky[0], 'CenterAz': CenterAz[0],
                     'CenterEl': CenterEl[0], 'ProjectionAltitude': ProjectionAltitude[0],
                     'ProjectionLon': ProjectionLon[0]}

        # size in pixels (width, height)
        self.size = rows[0], cols[0]

        # mode setting
        self.mode = "I;16L"

        # data descriptor
        self.tile = [
            ("raw", (0, 0) + self.size, 600, (self.mode, 0, 1))
        ]

Image.register_open("TIF", TifImageFile)

Image.register_extension("TIF", ".tif")
