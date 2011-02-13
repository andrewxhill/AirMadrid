#!/usr/bin/env python26
from pylab import *
import math
import datetime
import os
import numpy as np
import random
import subprocess
#from generate_tiles import *

VALUE_FILE = 'airdata2011.csv'
COORDINATE_FILE = 'estaciones.csv'
TILE_DIR = 'tiles/'
         
def GetCoordinates(fname):
    f = open(fname,'r')
    f.readline() #remove header
    coords = {}
    for line in f.readlines():
        line = line.strip().split(",")
        """
        gp = GoogleProjection()
        ll = fromLLtoPixel([line[2],line[3]],9)
        print ll
        """
        try:
            coords[int(line[0])] = {'lat': float(line[2]),'lon':float(line[3]),'name':line[1]}
        except:
            pass
    return coords

class Stations():
    maxVal = -1000000
    minVal = 10000000000
    coords = {}
    cmap = None
    
    def __init__(self,coords,fname):
        self.coords = coords
        self.datafile = fname
        self.hrData = self._getDailyData(self.datafile)
        self._imgLimits()
        
    def _imgLimits(self):    
        self.xmax,self.xmin = -100000000,100000000
        self.ymax,self.ymin = -100000000,100000000
        for a,i in self.coords.items():
            if float(i['lon']) > self.xmax:
                self.xmax = i['lon']
            if float(i['lon']) < self.xmin:
                self.xmin = i['lon']
            if float(i['lat']) > self.ymax:
                self.ymax = i['lat']
            if float(i['lat']) < self.ymin:
                self.ymin = i['lat']
        self.zeros = []
        
        
        #sets the final tile dimensions
        
        self.xdim = math.fabs(self.xmax - self.xmin)
        self.ydim = math.fabs(self.ymax - self.ymin)
        """
        edge = self.xdim if self.xdim > self.ydim else self.ydim
        self.xdim = edge
        self.ydim = edge
        """
        
        xfrac = (self.xmax - self.xmin)/7
        yfrac = (self.ymax - self.ymin)/7
        self.xfrac = xfrac
        self.yfrac = yfrac
        
        self.zeros.append({'lat':self.ymin + yfrac,'lon':self.xmin -xfrac})
        self.zeros.append({'lat':self.ymin - yfrac,'lon':self.xmin +xfrac})
        
        self.zeros.append({'lat':self.ymax - yfrac,'lon':self.xmin -xfrac})
        self.zeros.append({'lat':self.ymax + yfrac,'lon':self.xmin +xfrac})
        
        self.zeros.append({'lat':self.ymax - yfrac,'lon':self.xmax +xfrac})
        self.zeros.append({'lat':self.ymax + yfrac,'lon':self.xmax -xfrac})
        
        self.zeros.append({'lat':self.ymin - yfrac,'lon':self.xmax -xfrac})
        self.zeros.append({'lat':self.ymin + yfrac,'lon':self.xmax +xfrac})
        
        #self.ydim = int(1000 * (yfrac * 20)/(20*(xfrac+yfrac)))
        dim = self.xdim if self.xdim > self.ydim else self.ydim
        dim = dim/2
        xmid = (self.xmax+self.xmin)/2
        ymid = (self.ymax+self.ymin)/2
        self.xmax = xmid + dim
        self.xmin = xmid - dim
        self.ymax = ymid + dim
        self.ymin = ymid - dim
        self.zeros.append({'lat':self.ymin - yfrac,'lon':self.xmax +xfrac})
        self.zeros.append({'lat':self.ymin - yfrac,'lon':self.xmin -xfrac})
        self.zeros.append({'lat':self.ymax + yfrac,'lon':self.xmin -xfrac})
        self.zeros.append({'lat':self.ymax + yfrac,'lon':self.xmax +xfrac})
        
        """
        self.zeros.append({'lat':self.ymin - 2*yfrac,'lon':(self.xmax + self.xmin)/2})
        self.zeros.append({'lat':self.ymax + 2*yfrac,'lon':(self.xmax + self.xmin)/2})
        
        self.zeros.append({'lat':(self.ymax + self.ymin)/2,'lon':self.xmax - 2*xfrac})
        self.zeros.append({'lat':(self.ymin - self.ymin)/2,'lon':self.xmax - 2*xfrac})
        """
    def _getDailyData(self,fname):
        f = open(fname,'r')
        f.readline() #remove header
        hrData = {} #format = hr-day-yr
        for line in f.readlines():
            line = line.strip()
            line = line.split(',')
            station = int(line[0])
            if int(station) in self.coords.keys():
                dt = line[1].strip()
                dt = dt.split(' ')
                dtt = "%s %s %s" % (dt[2], dt[1], dt[4])
                yr = int(dt[4])
                day = int(datetime.datetime.strptime(dtt,"%d %b %Y").timetuple().tm_yday)
                hr = int(dt[3].split(':')[0])
                m = line[2]
        
                v = float(line[3])
                if v > self.maxVal:
                    self.maxVal = v
                elif v < self.minVal:
                    self.minVal = v
            
                try:
                    hrData[yr]
                except:
                    hrData[yr] = {}
                try:
                    hrData[yr][day]
                except:
                    hrData[yr][day] = {}
                try:
                    hrData[yr][day][hr]
                except:
                    hrData[yr][day][hr] = {}
                hrData[yr][day][hr][station] = v
        print self.minVal,self.maxVal
        return hrData
        
    def _getCMAP(self):
        
        x, y, z = [],[],[]
        for a in self.zeros:
            x.append(a['lon'])
            y.append(a['lat'])
            z.append(random.randint(10,100))
            
        z[0] = self.minVal
        z[-1] = self.maxVal
        
        xmin, xmax = min(x),max(x)
        ymin, ymax = min(y),max(y)
        x,y,z = np.array(x),np.array(y),np.array(z)
        
        nx = 500
        ny = 500
        
        xi = np.linspace(xmin, xmax, nx)
        yi = np.linspace(ymin, ymax, ny)
        xi, yi = np.meshgrid(xi, yi)
        
        zi = mlab.griddata(x,y,z,xi,yi)
        
        fig = plt.figure()
        pcm = plt.pcolormesh(xi,yi,zi,cmap='binary')
        self.cmap = pcm.get_cmap()
        plt.colorbar()
        fig.savefig('colorbar.png')
        
        
    def _renderLayer(self,yr,dy,hr,lyrName):
        x1, y1, z1 = [],[],[]
        for s,v in self.hrData[yr][dy][hr].items():
            x1.append(self.coords[s]['lon'])
            y1.append(self.coords[s]['lat'])
            z1.append(float(v))
        for a in self.zeros:
            x1.append(a['lon'])
            y1.append(a['lat'])
            z1.append(self.minVal)
            
        x,y,z = np.array(x1),np.array(y1),np.array(z1)
        
        nx,ny = 1000,1000
        
        xmin, xmax = self.xmin-self.xfrac, self.xmax+self.xfrac
        ymin, ymax = self.ymin-self.yfrac, self.ymax+self.yfrac
        
        xi = np.linspace(xmin, xmax, nx)
        yi = np.linspace(ymin, ymax, ny)
        xi, yi = np.meshgrid(xi, yi)
        
        zi = mlab.griddata(x,y,z,xi,yi)
        
        fig = plt.figure()
        fig.set_size_inches(10,10)
        
        ax = fig.add_subplot(111,frame_on=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.pcolormesh(xi,yi,zi,cmap=self.cmap)
        
        #plt.scatter(x,y,c=z)
        #plt.axis([xmin, xmax, ymin, ymax])
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(lyrName + '.png', bbox_inches=extent,dpi=100)
        
        #continue here
        
    def _renderTIF(self,lyrName):
        nulfp= open("error.log", "a+")
        p = subprocess.Popen(
            ["gdal_translate",
             "-a_srs",
             "EPSG:4326",
             "-a_ullr",
             "%s" % self.xmin,
             "%s" % self.ymax,
             "%s" % self.xmax,
             "%s" % self.ymin,
             "%s.png" % lyrName,
             "%s_4326.tif" % lyrName
            ], stderr=nulfp)
        p.wait()
        p = subprocess.Popen(
            ["gdalwarp",
             "-of",
             "GTiff",
             "-t_srs",
             "epsg:900913",
             "%s_4326.tif" % lyrName,
             "%s_900913.tif" % lyrName
            ], stderr=nulfp)
        p.wait()
        
    def createMaps(self,test=False):
        if self.cmap == None:
            self._getCMAP()
            
        if test:
            yr, dy, hr = 2011, 1, 8
            lyrDir = TILE_DIR + "%s/%s/" % (yr,dy)
            if not os.path.exists(lyrDir):
                os.makedirs(lyrDir)
            for hr in self.hrData[yr][dy].keys():
                lyrName = lyrDir + "%s" % hr
                self._renderLayer(yr,dy,hr,lyrName)
                self._renderTIF(lyrName)
                
        else:
            for yr in self.hrData.keys():
                for dy in self.hrData[yr].keys():
                    lyrDir = TILE_DIR + "%s/%s/" % (yr,dy)
                    if not os.path.exists(lyrDir):
                        os.makedirs(lyrDir)
                    for hr in self.hrData[yr][dy].keys():
                        lyrName = lyrDir + "%s" % hr
        
if __name__ == "__main__":  
    coords = GetCoordinates(COORDINATE_FILE)
    stations = Stations(coords,VALUE_FILE)
    stations.createMaps(test=True)