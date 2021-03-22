import sys
import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''
CREDIT FOR MAJORITY OF spefile CODE GOES TO https://github.com/timrae/pylase/blob/master/winspec.py
I MERELY ADAPTED AND MODIFIED AS NEEDED FOR A MUCH MORE SIMPLIFIED SOLUTION IN PYTHON 3
'''

def read_spe(spefile):

    # open SPE file as binary input
    spe = open(spefile, "rb")

    # Header length is a fixed number
    nBytesInHeader = 4100

    # Read the entire header
    header = spe.read(nBytesInHeader)

    # Date of the observation
    # (format is MM/DD/YY e.g. 02/25/21)
    date = struct.unpack_from("9s", header, offset=20)[0].decode('utf-8','replace')

    # Exposure time (float)
    exp_sec = struct.unpack_from("f", header, offset=10)[0]

    # Data type (0=float, 1=long integer, 2=integer, 3=unsigned int)
    data_type = struct.unpack_from("h", header, offset=108)[0]

    # CCD Chip Temperature (Degrees C)
    # detectorTemperature = struct.unpack_from("f", header, offset=36)[0]

    accumulations = struct.unpack_from("l", header, offset=668)[0]

    #total run time (in seconds)
    total_collection_time = exp_sec * accumulations

    # print ("exp_date     = ", date)
    # print ("exp_sec      = ", exp_sec)
    # print ("accums       = ", accumulations)
    # print ("run time (s) = ", total_collection_time)
    # print ("data_type    = ", data_type)
    # print ("detectorTemperature [C] = ", detectorTemperature)

    # Number of pixels on x-axis and y-axis
    nx = struct.unpack_from("H", header, offset=42)[0]
    ny = struct.unpack_from("H", header, offset=656)[0]

    # Number of image frames in this SPE file
    nframes = struct.unpack_from("l", header, offset=1446)[0]

    # print ("nx, ny, nframes = ", nx, ", ", ny, ", ", nframes)

    if data_type == 0:
        # float (4 bytes)
        dataTypeStr = "f"  #untested
        bytesPerPixel = 4
        dtype = "float32"
    elif data_type == 1:
        # long (4 bytes)
        dataTypeStr = "l"  #untested
        bytesPerPixel = 4
        dtype = "int32"
    elif data_type == 2:
        # short (2 bytes)
        dataTypeStr = "h"  #untested
        bytesPerPixel = 2
        dtype = "int32"
    elif data_type == 3:  
        # unsigned short (2 bytes)
        dataTypeStr = "H"  # 16 bits in python on intel mac
        bytesPerPixel = 2
        dtype = "int32"  # for numpy.array().
        # other options include:
        # IntN, UintN, where N = 8,16,32 or 64
        # and Float32, Float64, Complex64, Complex128
        # but need to verify that pyfits._ImageBaseHDU.ImgCode cna handle it
        # right now, ImgCode must be float32, float64, int16, int32, int64 or uint8
    else:
        print ("unknown data type")
        print ("returning...")
        sys.exit()

    npixels = nx*ny
    npixStr = str(npixels)
    fmtStr  = npixStr+dataTypeStr
    # print ("fmtStr = ", fmtStr)

    # How many bytes per image?
    nbytesPerFrame = npixels*bytesPerPixel
    # print ("nbytesPerFrame = ", nbytesPerFrame)

    """Start of X Calibration Structure (although I added things to it that I thought were relevant,
        like the center wavelength..."""
    xcalib = {}

    #float SpecCenterWlNm # 72 Center Wavelength in Nm
    xcalib['SpecCenterWlNm'] = struct.unpack_from("f", header, offset=72)[0]

    #SHORT SpecGlueFlag 76 T/F File is Glued
    xcalib['SpecGlueFlag'] = bool( struct.unpack_from("h", header, offset=76)[0] )

    #float SpecGlueStartWlNm 78 Starting Wavelength in Nm
    xcalib['SpecGlueStartWlNm'] = struct.unpack_from("f", header, offset=78)[0]

    #float SpecGlueEndWlNm 82 Starting Wavelength in Nm
    xcalib['SpecGlueEndWlNm'] = struct.unpack_from("f", header, offset=82)[0]

    #float SpecGlueMinOvrlpNm 86 Minimum Overlap in Nm
    xcalib['SpecGlueMinOvrlpNm'] = struct.unpack_from("f", header, offset=86)[0]

    #float SpecGlueFinalResNm 90 Final Resolution in Nm
    xcalib['SpecGlueFinalResNm'] = struct.unpack_from("f", header, offset=90)[0]

    #  short   BackGrndApplied              150  1 if background subtraction done
    xcalib['BackgroundApplied'] = struct.unpack_from("h", header, offset=150)[0]
    BackgroundApplied=False
    if xcalib['BackgroundApplied']==1: BackgroundApplied=True

    #  float   SpecGrooves                  650  Spectrograph Grating Grooves
    xcalib['SpecGrooves'] = struct.unpack_from("f", header, offset=650)[0]

    #  short   flatFieldApplied             706  1 if flat field was applied.
    xcalib['flatFieldApplied'] = struct.unpack_from("h", header, offset=706)[0]
    flatFieldApplied=False
    if xcalib['flatFieldApplied']==1: flatFieldApplied=True

    #double offset # 3000 offset for absolute data scaling */
    xcalib['offset'] = struct.unpack_from("d", header, offset=3000)[0]

    #double factor # 3008 factor for absolute data scaling */
    xcalib['factor'] = struct.unpack_from("d", header, offset=3008)[0]

    #char current_unit # 3016 selected scaling unit */
    xcalib['current_unit'] = ord(struct.unpack_from("c", header, offset=3016)[0])

    #char calib_valid # 3098 flag if calibration is valid */
    xcalib['calib_valid'] = ord(struct.unpack_from("c", header, offset=3098)[0])

    #char input_unit # 3099 current input units for */
    xcalib['input_unit'] = ord(struct.unpack_from("c", header, offset=3099)[0])
    """/* "calib_value" */"""

    #char polynom_unit # 3100 linear UNIT and used */
    xcalib['polynom_unit'] = ord(struct.unpack_from("c", header, offset=3100)[0])
    """/* in the "polynom_coeff" */"""

    #char polynom_order # 3101 ORDER of calibration POLYNOM */
    xcalib['polynom_order'] = ord(struct.unpack_from("c", header, offset=3101)[0])

    #char calib_count # 3102 valid calibration data pairs */
    xcalib['calib_count'] = ord(struct.unpack_from("c", header, offset=3102)[0])

    #double pixel_position[10];/* 3103 pixel pos. of calibration data */
    xcalib['pixel_position'] = struct.unpack_from("10d", header, offset=3103)

    #double calib_value[10] # 3183 calibration VALUE at above pos */
    xcalib['calib_value'] = struct.unpack_from("10d", header, offset=3183)

    #double polynom_coeff[6] # 3263 polynom COEFFICIENTS */
    xcalib['polynom_coeff'] = struct.unpack_from("6d", header, offset=3263)

    #double laser_position # 3311 laser wavenumber for relativ WN */
    xcalib['laser_position'] = struct.unpack_from("d", header, offset=3311)[0]

    #setup data dictionary
    spedict = {'data':[], 
                'EXPOSURE':exp_sec,
                'SPEFNAME':spefile,
                'OBSDATE':date,
                'XCALIB':xcalib,
                'ACCUMULATIONS':accumulations,
                'FLATFIELD':flatFieldApplied,
                'BACKGROUND':BackgroundApplied
                }

    for ii in range(nframes):
        iistr = str(ii)
        data = spe.read(nbytesPerFrame)
        # read pixel values into a 1-D numpy array. the "=" forces it to use
        # standard python datatype size (4bytes for 'l') rather than native
        # (which on 64bit is 8bytes for 'l', for example).
        # See http://docs.python.org/library/struct.html
        dataArr = np.array(struct.unpack_from("="+fmtStr, data, offset=0),
                            dtype=dtype)

        # Resize array to nx by ny pixels
        # notice order... (y,x)
        dataArr.resize((ny, nx))
        # print (dataArr.shape)

        # Push this image frame data onto the end of the list of images
        # but first cast the datatype to float (if it's not already)
        # this isn't necessary, but shouldn't hurt and could save me
        # from doing integer math when i really meant floating-point...
        spedict['data'].append(dataArr.astype(float))

    xpixel_calib_coeff_array = np.flip(np.asarray(spedict['XCALIB']['polynom_coeff']))
    xpixel_array = range(1,1+int(spedict['XCALIB']['pixel_position'][1]))
    # print(xpixel_calib_coeff_array)
    # print(xpixel_array)

    # Finds wavelength from pixels using calibration coefficients
    wavelengthData = np.polyval(xpixel_calib_coeff_array,xpixel_array)
    # print("Wavelength data is... ", wavelengthData)
    # print("Does this roughly match?: ", spedict['XCALIB']['calib_value'][:2])

    # Only pulls first spectrum. (we generally only collect one spectrum and accumulate anyway)
    # Could easily change to take average of multiple spectra.
    intensityData = spedict['data'][0][0]

    #test to see that x and y data is correct
    # plt.plot(wavelengthData,intensityData)
    # plt.show()

    return wavelengthData, intensityData, spedict

path = '/home/jkusz/github/king_jason/Project-2/CFS1-2.SPE'
xdata, ydata, spedict = read_spe(path)
# print("xdata: ", xdata, "\nydata: ", ydata)
# print(spedict)