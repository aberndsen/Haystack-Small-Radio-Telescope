#!/usr/bin/env python
"""
A small routine to plot SRT data files. 
The program should be smart enough to sort multiple files by date,
 and plot the sorted data.

One could also import this file to access the SRTdata module,
which will parse most of the data for you.
Eg. import srt_dataparser as SRT
d = SRT.SRTdata('inputfile.rad')

Then use d's or SRT's useful functions for some simple analysis.
d.data (ntsteps x nchans) and 
d.data_info (ntsteps x (time az alt daz dalt freq mode))

"""
import datetime
import os
import time
from optparse import OptionParser
import pylab as plt
import numpy as np
import re
import ephem

#create a small font for the legend
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D 


def main(files,opts):
    """
    Given our list of data-file objects "files" generate the plots 
    specified in opts

    """
    if len(files) > 1:
        print "Using sorted list of files (name,date,duration):"
    for f in files:
        print "\t (%s, %s, %s)" % \
            (f.fname,f.tsteps[0],f.tsteps[-1]-f.tsteps[0])

#make a plot of channel vs time?
    if opts.chanplot:
        plot_channels(files,channels=opts.chanlist,ignoreYMD=opts.ignoreYMD)

    if opts.spectrum:
        plot_spectrum(files)

    if opts.plottracks:
        plot_tracks(files)

    if opts.plottrackoffsets:
        plot_trackoffsets(files, threshold=opts.threshold)

    if opts.plottrackoffsets3d:
        plot_trackoffsets3d(files, opts.fracofyear)

    if opts.plotnpointoffsets:
        plot_npointoffsets(files)

    if opts.plotnpointoffsets3d:
        plot_npointoffsets3d(files, opts.fracofyear)

############################################
# SRTdata 
    

class SRTdata(object):
    """
    a class to store each SRT data file. 

    The object reads in the data, storing:
    timesteps in self.tsteps
    frequency intensity in self.data (ntsteps,nchannels)

    other observing information in 
    self.data_info : 
          (a numpy.recarray holding date, az, alt, az_off, alt_off, lofreq, 
                                    chan_bw, mode, nchan)
    self.data : a numpy array (ntsteps x nchans)
    self.obs : an ephem.Observer for the lat/lon of the data (epoch 1950) 
    self.body0 : an ephem.FixedBody for the start of the observation
    self.body1 : en ephemFixedBody for the end of the observation
    

    """

    
#convert the time string into a datetime object, and the other cols
    func = lambda s: datetime.datetime.strptime(s,"%Y:%j:%H:%M:%S")
    converters = {0: func}
    srtype = [('date',datetime.datetime),('az',float),('alt',float),\
                  ('az_off',float),('alt_off',float),\
                  ('lofreq',float),('chan_bw',float),\
                  ('mode',int),('nchan',int)]

    
    def __init__(self,fname,elev=105.):
        """
        given a .rad file, create a SRTdata object for convenient access to the data. 
        The 'elev' argument defaults to Vancouver's configuration, so you may want to change it.
        It will have a small effect on the computed Az/Alts

        """
        self.fname = fname #name of SRT data file
        self.obs = ephem.Observer() #create an Observer to help with coord conversions
        self.obs.elevation = 105.
        self.obs.epoch = ephem.B1950 #'1950'
        
        self.srtypes = SRTdata.srtype


# ignore above, since the data file includes the number of channels anyways.
        self.get_nchan()


# store the telescope information in 'data_info'
# note: genfromtxt (using dtype) returns a 1d array, where
# each column is accessed as data['name of column']
        self.data_info = np.genfromtxt(fname,
                                       dtype=self.srtypes,
                                       skip_header=1,
                                       comments='*',
                                       usecols=range(len(self.srtypes)),
                                       converters=SRTdata.converters)

# store the channel data in 'data'
        usecols = [i + len(self.srtypes) for i in range(self.nchan)]
        self.data = np.genfromtxt(fname,
                                  dtype=None,
                                  comments='*',
                                  skip_header=1,
                                  usecols=usecols,
                                  filling_values = 0.,
                                  invalid_raise = False
                                  )
        ndata = self.data.shape[0]
        self.data_info = self.data_info[0:ndata]

        self.tsteps = self.data_info['date']

# determine some basic information about the observation
# including 'telescope mode', 'telescope lat/lon'
# later is used for converting (Az,Alt) to (Ra,Dec)
        self.get_info()



    def get_info(self):
        """
        read the first data line of the SRT data file to determine
        the observing mode, latitude, and longitude

        This routine should be run last, since we use the data arrays

        """
        self.has_npoint_scan = False

        f = open(self.fname)

# determine the lat/lon from first line of data file
        r1 = f.readline()
        lat, lon = re.findall("\d+.\d*",r1)
        if 'LONGW' in r1:
            lon = '-%s' % lon
        self.obs.lat = lat
        self.obs.lon = lon
        npoints = []
# determine date and mode of observation from second line 
# and the data rows just after npointing
        line_num = 0
        for line in  f:
            if re.findall('NPOINT',line):
                self.has_npoint_scan = True
#                print "Observation %s has an npoint scan" % self.fname
                lsplt = line.split()
                npoints.append((line_num,float(lsplt[5]),float(lsplt[6])))

            if line[0] != '*': 
                self.mode = int(line.split()[7])
                self.obs.date = datetime.datetime.strptime(line.split()[0],\
                                                               "%Y:%j:%H:%M:%S")
            else:
# the data and data_info arrays skip over '*', so don't count these rows
                line_num += -1

            line_num += 1

# prettify the npoint data,
# recording the datetime, az, alt, az_off, alt_off
        self.npoints = []
        for idx,az_off,alt_off in npoints:
            dat = self.data_info[idx-1]
            self.npoints.append((dat[0],dat[1],dat[2],az_off,alt_off))

        f.close()

# create ephem.FixedBody's for the start and end of the observation
# recall, angles in degrees
        azs = (self.data_info['az'] + self.data_info['az_off'])
        alts = (self.data_info['alt'] + self.data_info['alt_off'])

# takes str(deg), or rad
        ra0,dec0 = self.obs.radec_of(str(azs[0]),str(alts[0])) 
        self.body0 = ephem.FixedBody()
        self.body0.name = 'startObs'
        self.body0._ra = str(ra0)
        self.body0._dec = str(dec0)
        self.body0._epoch = ephem.B1950 #'1950'

        ra1,dec1 = self.obs.radec_of(azs[-1],alts[-1])
        self.body1 = ephem.FixedBody()
        self.body1.name = 'endObs'
        self.body1._ra = ra1
        self.body1._dec = dec1
        self.body1._epoch = ephem.B1950 #'1950'

    def get_nchan(self):
        """
        read the first data line of the SRT data file to determine
        the number of frequency channels. We use this to help construct the 'dtype'
        needed to read in the data
        
        """
        f = open(self.fname)
        for line in f.readlines():
            if line[0] == '*':
                continue
            else:
                self.nchan = int(line.split()[8])
                break
        f.close()

    def plot_track(self,sun=True,xlim=(),ylim=()):
        """
        Call this modules 'plottracks' routine, 
        for this data set only
        
        args: accepts same arguments as 'plottracks'
        
        """
        plot_tracks([self],sun=sun,xlim=xlim,ylim=ylim)

######### end class SRTdata ###############



def plot_spectrum(files):
    """
    plot the spectrum

    """
    print "Plotting spectrum. Close window to continue...\n"
    
    nplots = len(files)
    ncols = 2
    nrows = nplots/ncols
    if nrows*ncols < nplots:
        nrows += 1
    idx = 0
    for f in files:
        ax = plt.subplot(nrows,ncols,idx)#,subplots_adjust(wspace=0,hspace=0))
        img = ax.imshow(f.data.transpose(), 
                          aspect='auto',
                          cmap='gray',
                          origin='lower')

        plt.title(f.fname)
# make a legend
        plt.colorbar(img, cmap='gray')
        if idx == 0:
            plt.xlabel('timestep')
            plt.ylabel('frequency channel')
        idx += 1
    plt.show()
    
def plot_channels(files,channels=None,chanbn=1,ignoreYMD=False):
    """
    given a list of SRTdata objects/files,
    plot the frequency vs channel

    args:
    channels: a list channels to plot
              default: 'None' plots the average (Intensity)

    chanbn: ** Not yet implemented **
            average the channels into this many bins
            default: 1 does no binning
    ignoreYMD: True/False ignore the year/month/day timestamps?
             useful if comparing tracks tacking on sequential days
             (like the day before and day-of an eclipse)
             default: False

    """
    print "Plotting channel intensities. Close window to continue...\n"
    if not channels:
        for f in files:
            dataavg = f.data.mean(axis=1)
            if ignoreYMD:
                tsteps = [replaceYMD(t) for t in f.tsteps]
            else:
                tsteps = f.tsteps
            plt.plot(tsteps,dataavg,label=f.fname)
            
        plt.xlabel('timestamp (UTC)')
        plt.ylabel('temperature [k]')
        plt.legend(prop=fontP)
        plt.title('Intensity vs time')
    
    else:
        for f in files:
#first plot the frequency average
            plt.plot(f.tsteps,f.data.mean(axis=1),label='average intensity')
            for chan in [int(i) for i in channels.split(',')]:
                if ignoreYMD:
                    tsteps = [replaceYMD(t) for t in f.tsteps]
                else:
                    tsteps = f.tsteps
                plt.plot(tsteps,f.data[:,chan],label='channel %s' % chan)
            
    plt.legend(prop=fontP)
    plt.show()
    

def plot_tracks(files=[],sun=True,xlim=(),ylim=()):
    """
    given a list of SRTdata objects,
    plot their az/alt tracks
    
    args:
    files: list of SRTdata objects
    sun: True/False, plot the sun's track?
    xlim: (xmin, xmax) [deg]
    ylim: (ymin, ymax) [deg]

    """

    if sun: suno = ephem.Sun(epoch=ephem.B1950)

# keep track of separation from the sun (if sun=True)    
    maxsep = 0
    if files:
        for f in files:

            azs = (f.data_info['az'] + f.data_info['az_off'])
            alts = (f.data_info['alt'] + f.data_info['alt_off'])
            #plot start and end, keeping same color for each file
            ax = plt.plot(azs[0],alts[0],'*',markersize=15)[0]
            clr = ax.get_color()
            plt.plot(azs[-1],alts[-1],'%sx'%clr,markersize=12)
            if len(files) > 5:
                plt.plot(azs,alts,'%s+'%clr,label=f.fname)
            else:
                plt.plot(azs,alts,clr,label=f.fname)
            if sun:
                tsteps = f.tsteps
                sazs = []
                salts = [] 
                tmax = 0
                for ti,tv in enumerate(tsteps):
                    f.obs.date = tv
                    suno.compute(f.obs)
# create object at latest telescope Alt/Az
# so we can calculate the separation from sun
                    tmpbod = ephem.FixedBody()
                    ra,dec = f.obs.radec_of(str(azs[ti]),str(alts[ti])) 
                    tmpbod._ra = str(ra)
                    tmpbod._dec = str(dec)
                    tmpbod._epoch = ephem.B1950 #'1950'
                    tmpbod.compute(f.obs)
                    if ephem.separation(tmpbod,suno) > maxsep:
                        maxsep = ephem.separation(tmpbod,suno)#max(maxsep,ephem.separation(tmpbod,suno))
                        tmax = tv
                    sazs.append(suno.az*180./np.pi)
                    salts.append(suno.alt*180./np.pi)
                    
# skip some points
                skp = int(np.log(len(tsteps)))*4
                ax = plt.plot(sazs[::skp],salts[::skp],'+',label='sun during: %s' % f.fname)[0]
                clr = ax.get_color()
                #special marks for start and end:
                plt.plot(sazs[0],salts[0],'%s*' % clr,markersize=15)
                plt.plot(sazs[-1],salts[-1],'%sx' % clr,markersize=12)

    plt.xlabel('Azimuth [deg]')
    plt.ylabel('Altitude [deg]')
    if xlim:
         plt.xlim(xlim)
    if ylim:
        plt.ylim((0,95))
                                 
    if sun:
        plt.title('Tracks [B1950 coords]\n (max separation of %s at %s)'\
                      % (maxsep,tmax))
    else:
        plt.title('Tracks [B1950 coords]')
    plt.legend(bbox_to_anchor=(1.23, 1.0),prop=fontP)
    plt.show()


def plot_trackoffsets(files, threshold=360.):
    """
    given a list of SRTdata objects, plot the alt-az pointing
    offset from the sun.

    Args:
    files : list of srt data files/objects
    threshold : only plot points within this many degrees.

    """

    sun = ephem.Sun(epoch=ephem.B1950)
    tazs = []
    talts = []
    trk_az_off = []
    trk_alt_off = []
    seps = []
    for f in files:
    
        azs = (f.data_info['az'] + f.data_info['az_off'])
        alts = (f.data_info['alt'] + f.data_info['alt_off'])
        for az in azs:
            tazs.append(az)
        for alt in alts:
            talts.append(alt)

        tsteps = f.tsteps
        tmax = 0
#find (alt,az) of sun, where we were pointing, and their separation
        for ti,tv in enumerate(tsteps):
            f.obs.date = tv
            sun.compute(f.obs)
# create object at latest telescope Alt/Az
# so we can calculate the separation from sun
            srtbod = ephem.FixedBody(epoch=ephem.B1950)
            ra,dec = f.obs.radec_of(str(azs[ti]),str(alts[ti]))  #current ra/dec in B1950
            srtbod._ra = str(ra)
            srtbod._dec = str(dec)
            srtbod.compute(f.obs)
#            trk_az_off.append(srtbod.az - sun.az) #rad
#            trk_alt_off.append(srtbod.alt - sun.alt) #rad

# line up the objects to compute offset in each direction
            trk_az_off.append(ephem.separation((srtbod.az, srtbod.alt), (sun.az, srtbod.alt)))
            trk_alt_off.append(ephem.separation((srtbod.az, srtbod.alt), (srtbod.az, sun.alt)))
            seps.append(float(ephem.separation((srtbod.az, srtbod.alt), (sun.az, sun.alt))))  #ra sep.


    idcs = np.where(np.array(seps) < threshold*np.pi/180.)[0]
#convert rad --> deg
    trk_az_off = np.array(trk_az_off)*180./np.pi
    trk_alt_off = np.array(trk_alt_off)*180./np.pi
    tazs = np.array(tazs)
    talts = np.array(talts)
#plot the az and alt offsets
    plt.subplot(2,1,1)
    plt.plot(tazs[idcs],trk_az_off[idcs],'b+')
    plt.xlabel('telescope azimuth')
    plt.ylabel('aziumuthal separation [deg]')
    plt.title('telescope.az - sun.az')
    plt.subplot(2,1,2)
    plt.plot(talts[idcs],trk_alt_off[idcs],'b+')
    plt.xlabel('telescope altitude')
    plt.ylabel('altitude separation [deg]')
    plt.title('telescope.alt - sun.alt')

    plt.show()

def plot_trackoffsets3d(files, foy=False, threshold=360.):
    """
    given a list of SRTdata objects, plot the alt-az pointing
    offset from the sun using (az/alt or fraction of year as 3rd axis)/.

    Args:
    files = list of srt data files/objects
    foy = use fraction of year as third axis (True), or az/alt (default)
    threshold = only plot points within this many degrees

    """

#prep the data
    sun = ephem.Sun(epoch=ephem.B1950)
    tazs = []
    talts = []
    trk_az_off = []
    trk_alt_off = []
    foys = []
    seps = []
    for f in files:
    
        azs = (f.data_info['az'] + f.data_info['az_off'])
        alts = (f.data_info['alt'] + f.data_info['alt_off'])
        for az in azs:
            tazs.append(az)
        for alt in alts:
            talts.append(alt)

        tsteps = f.tsteps
        tmax = 0
#find (alt,az) of sun, where we were pointing, and their separation
        for ti,tv in enumerate(tsteps):
            f.obs.date = tv
            sun.compute(f.obs)
# create object at latest telescope Alt/Az
# so we can calculate the separation from sun
            srtbod = ephem.FixedBody()
            ra,dec = f.obs.radec_of(str(azs[ti]),str(alts[ti])) 
            srtbod._ra = str(ra)
            srtbod._dec = str(dec)
            srtbod._epoch = ephem.B1950 #'1950'
            srtbod.compute(f.obs)
#            trk_az_off.append(srtbod.az - sun.az) #rad
#            trk_alt_off.append(srtbod.alt - sun.alt) #rad

# line up the objects to compute offset in each direction
            trk_az_off.append(ephem.separation((srtbod.az, srtbod.alt), (sun.az, srtbod.alt)))
            trk_alt_off.append(ephem.separation((srtbod.az, srtbod.alt), (srtbod.az, sun.alt)))
            seps.append(ephem.separation(srtbod, sun))
            foys.append(year_fraction(tv))
            
    idcs = np.where(np.array(seps) < threshold*np.pi/180.)[0]
#convert rad --> deg
    trk_az_off = np.array(trk_az_off)*180./np.pi
    trk_alt_off = np.array(trk_alt_off)*180./np.pi
    tazs = np.array(tazs)
    talts = np.array(talts)
#plot the az and alt offsets
    fig = plt.figure(1)
    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    if not foy:
        ax1.scatter(tazs[idcs], talts[idcs], trk_az_off[idcs], '+')
        ax1.set_ylabel('altitude [deg]')
    else:
        ax1.scatter(tazs[idcs], foys[idcs], trk_az_off[idcs], '+')
        ax1.set_ylabel('fraction of year')
    ax1.set_xlabel('telescope azimuth')
    ax1.set_zlabel('aziumuthal separation [deg]')
    ax1.set_title('telescope.az - sun.az')
    print "Close window to continue..."
    plt.show()
    plt.clf()

    fig = plt.figure(1)
    ax2 = fig.add_subplot(1, 1, 1, projection='3d')
    if not foy:
        ax2.scatter(tazs[idcs], talts[idcs], trk_alt_off[idcs], '+')
        ax2.set_ylabel('altitude [deg]')
    else:
        ax2.scatter(tazs[idcs], foys[idcs], trk_alt_off[idcs], '+')
        ax2.set_ylabel('fraction of year')
    ax2.set_xlabel('telescope azimuth')
    ax2.set_zlabel('altitude separation [deg]')
    ax2.set_title('telescope.alt - sun.alt')
    plt.show()

def plot_npointoffsets(files):
    """
    given a list of SRTdata objects, plot the alt-az pointing
    offset from the sun as determined by the NPOINT scans

    """
    
    azs = []
    alts = []
    dazs = []
    dalts = []

    for f in files:
        npoints = f.npoints
        for t,az,alt,daz,dalt in npoints:
            azs.append(az)
            alts.append(alt)
            dazs.append(daz)
            dalts.append(dalt)

    plt.subplot(221)
    ax1 = plt.plot(azs,dazs,'+')
    plt.xlabel('azimuthal [deg]')
    plt.ylabel('azimuth offset')

    plt.subplot(222)
    ax3 = plt.plot(alts,dazs,'+')
    plt.xlabel('altitude [deg]')
    plt.ylabel('azimuthal offset')

    plt.subplot(223)
    ax2 = plt.plot(azs,dalts,'+')
    plt.xlabel('azimuth [deg]')
    plt.ylabel('altitude offset')

    plt.subplot(224)
    ax4 = plt.plot(alts,dalts,'+')
    plt.xlabel('altitude [deg]')
    plt.ylabel('azimuthal offset')
    
    plt.show()

def plot_npointoffsets3d(files, foy=False):
    """
    given a list of SRTdata objects, plot the alt-az pointing
    offset reported by the NPOINT scans of the sun (in 3d).

    Args:
    files : list of srt data files
    foy : use fraction of year as the 3rd dimensions [default = False]

    """
    
    azs = []
    alts = []
    dazs = []
    dalts = []
    foys = []

    for f in files:
        npoints = f.npoints
        for t, az, alt, daz, dalt in npoints:
            azs.append(az)
            alts.append(alt)
            dazs.append(daz)
            dalts.append(dalt)
            foys.append(year_fraction(t))

    fig = plt.figure(1)
    ax1 = fig.add_subplot(1,2,1, projection='3d')
#    ax1.plot(azs,dazs,'+')
    if not foy:
        ax1.scatter(azs, alts, dazs,'+')
        ax1.set_ylabel('altitude [deg]')
    else:
        ax1.scatter(azs, foys, dazs,'+')
        ax1.set_ylabel('fraction of year')
    ax1.set_title('az offset')
    ax1.set_xlabel('azimuth [deg]')
    ax1.set_zlabel('azimuthal offset')
        

    ax3 = fig.add_subplot(1,2,2, projection='3d')
#    ax3.plot(alts,dazs,'+')
    if not foy:
        ax3.scatter(alts, azs, dalts, '+')
        ax3.set_ylabel('azimuth [deg]')
    else:
        ax3.scatter(alts, foys, dalts, '+')
        ax3.set_ylabel('fraction of year')
    ax3.set_title('alt offset')
    ax3.set_xlabel('altitude [deg]')
    ax3.set_zlabel('altitude offset')

#    ax2 = fig.add_subplot(2,2,3, projection='3d')
#    ax2.plot(azs,dalts,'+')
#    ax2.set_xlabel('azimuth [deg]')
#    ax2.set_ylabel('altitude offset')
#    ax2.set_zlabel('altitude offset')

#    ax4 = fig.add_subplot(2,2,4, projection='3d')
#    ax4.plot(alts,dalts,'+')
#    ax4.set_xlabel('altitude [deg]')
#    ax4.set_ylabel('azimuthal offset')
#    ax4.set_zlabel('')
    
    plt.show()


#######################

def replaceYMD(date):
    """
    replace the year, month and date of an observation
    keeping the HH:MM:ss.

    Useful for comparing observations taking on sequential days, 
    such as the day before and day of an eclipse

    We set the Y/M/D params to the B1950-ish epoch for fun
    """
    return date.replace(year=1950,month=1,day=1)
    
def juldate2ephem(num):
    """Convert Julian date to ephem date, measured from noon, Dec. 31, 1899."""
    return ephem.date(num - 2415020.)

def ephem2juldate(num):
    """Convert ephem date (measured from noon, Dec. 31, 1899) to Julian date."""
    return float(num + 2415020.)

def ephem2mjd(num):
    """Convert ephem date to MJD (measured from 00:00 Nov 17, 1858)."""
    return num + 15019.5  #= - 2415020.+2400000.5

def mjd2ephem(num):
    """Convert MJD to ephem date"""
    return ephem.date(num - 15019.5)

def mjd2juldate(num):
    """Convert MJD date to Julian date."""
    return float(num + 2400000.5)

def juldate2mjd(num):
    """Convert Juliand date to MJD"""
    return float(num - 2400000.5)

def year_fraction(date):
    """return the fraction of the year for this datetime object"""
    def sinceEpoch(date): #seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return fraction
        

if __name__ == '__main__':
    
    usage = "usage: python srt_plotdata.py [options] *.rad"
    parser = OptionParser(usage=usage)
    parser.add_option("-s","--spectrum",
                      action='store_true',
                      dest="spectrum",
                      default=False,
                      help="Plot the spectrum")
    parser.add_option("-c","--intensity",
                      action='store_true',
                      dest="chanplot",
                      default=False,
                      help="Plot intensity vs time")
    parser.add_option("--chanlist",dest="chanlist",
                      default=None,
                      help="Plot the individual channel intensities."\
                          "Eg. --chanlist=1,3,55  (zero-indexed)")
    parser.add_option("--ignoreYMD",dest="ignoreYMD",
                      action='store_true',
                      default=False,
                      help="ignore the year/month/day timestamps? in intensity plots"\
                          "useful if comparing data on sequential days"\
                          "(like the day before and day-of an eclipse)"
                      )
    parser.add_option("-t","--plottracks",
                      action='store_true',
                      dest='plottracks',
                      default=False,
                      help="Plot AZ ALT tracks")
    parser.add_option("-d","--plottracksoffsets",
                      action='store_true',
                      dest='plottrackoffsets',
                      default=False,
                      help="Plot AZ ALT tracking offsets from sun")
    parser.add_option("-e","--plottracksoffsets3d",
                      action='store_true',
                      dest='plottrackoffsets3d',
                      default=False,
                      help="Plot AZ ALT tracking offsets from sun (in 3d)")
    parser.add_option("-n","--plotnpointoffsets",
                      action='store_true',
                      dest='plotnpointoffsets',
                      default=False,
                      help="Plot AZ ALT tracking offsets from sun, as determined by npoint scans")
    parser.add_option("-m","--plotnpointoffsets3d",
                      action='store_true',
                      dest='plotnpointoffsets3d',
                      default=False,
                      help="Plot AZ ALT tracking offsets from sun, as determined by npoint scans (in 3d)")
    parser.add_option("--fracofyear",
                      action='store_true',
                      dest='fracofyear',
                      default=False,
                      help="Plot the offsets, using fraction of year as 3rd dimension. (used in 3d plots only)")
    parser.add_option("--threshold",
                      dest='threshold',
                      default=360.,
                      type='float',
                      help="For all plots, only plot points within this distance [DEG] (helps cut npoint scanning).")

    (opts,args) = parser.parse_args()
    files = []
#do preliminary sort based on filename (later sorted by data timestamp)
    for f in sorted(args):
        if os.path.exists(f):
            files.append(SRTdata(f))
        else:
            print "Couldn't find %s. Skipping" % f

    print "THRS", opts.threshold, type(opts.threshold)
    if files:
        #sort based on the first timestamp
        files = sorted(files,key=lambda f:f.tsteps[0])
        main(files,opts)
    else:
        print "No data files found. Exiting..."
