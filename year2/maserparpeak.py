#!/Users/Lucas/anaconda2/bin/python

from astropy.time import Time as aTime
from astropy.coordinates import SkyCoord
from astropy import units as u
import datetime, difflib, numpy as np
from numpy import cos, sin, tan, pi
import pdb
import warnings
warnings.filterwarnings("ignore")

###############################################################################



class Maser:
    kind = 'maser' 
    def __init__(self, name, ra, dec, vel, flux, alias):
        self.name  = name  
        self.ra    = ra
        self.dec   = dec
        self.vel   = float(vel)
        self.cflux  = float(flux)
        self.alias = alias

def get_maser(template):
    masers = np.array([
        ['G232.620+0.996','07:32:09.79','-16:58:12.4', 22.9, 11.5,'s1' ],
        ['G287.371+0.644','10:48:04.44','-58:27:01.0', -1.9, 21.9,'s3' ],
        ['G309.921+0.479','13:50:41.78','-61:35:10.2',-57.9, 57.6,'s4' ],
        ['G323.740-0.263','15:31:45.45','-56:30:50.1',-50.4,346.6,'s5' ],
        ['G327.402+0.445','15:49:19.50','-53:45:13.9',-82.9, 37.7,'s6' ],
        ['G328.254-0.532','15:57:59.75','-53:58:00.4',-36.8, 20.9,'s7' ],
        ['G328.808+0.633','15:55:48.45','-52:43:06.6',-44.4, 30.6,'s8' ],
        ['G339.622-0.121','16:46:05.99','-45:36:43.3',-33.2, 23.3,'s9' ],
        ['G339.884-1.259','16:52:04.67','-46:08:34.2',-35.6,424.1,'s10'],
        ['G345.505+0.348','17:04:22.91','-40:44:21.7',-14.1, 24.5,'s11'],
        ['G291.274-0.709','11:11:53.35','-61:18:23.7',-30.7, 10.7,'s14'],
        ['G299.772-0.005','12:23:48.97','-62:42:25.3', -6.7, 12.3,'s15'],
        ['G318.948-0.196','15:00:55.40','-58:58:52.1',-36.3, 12.5,'s16'],
        ['G326.475+0.703','15:43:16.64','-54:07:14.6',-38.4, 13.5,'s17'],
        ['G328.237-0.547','15:57:58.28','-53:59:22.7',-44.7, 41.9,'s18'],
        ['G329.029-0.205','16:00:31.80','-53:12:49.6',-36.1, 11.1,'s19'],
        ['G332.295+2.280','16:05:41.72','-49:11:30.3',-23.7, 10.5,'s20'],
        ['G337.920-0.456','16:41:06.05','-47:07:02.5',-38.6, 12.7,'s21'],
        ['G345.010+1.792','16:56:47.58','-40:14:25.8',-17.0, 14.2,'s22'],
        ['G348.550-0.979','17:19:20.41','-39:03:51.6',-10.4, 10.5,'s23'],
        ['G352.630-1.067','17:31:13.91','-35:44:08.7', -3.3, 17.6,'s24']])
    try:
        match = difflib.get_close_matches(template, masers[:,0])[0]
        index = masers[:,0]==match
        return Maser(*masers[index,:][0])
    except IndexError:
        try: 
            match = difflib.get_close_matches(template, masers[:,-1])[0]
            index = masers[:,-1]==match
            return Maser(*masers[index,:][0])
        except IndexError:
            print('Cannot identify maser, defaulting to G232.62')
            return Maser(*masers[0,:])

def parpropmot(p,t,t_r,ra,dec):
    X  = cos(2*pi*(t-0.22))
    Y  = sin(2*pi*(t-0.22))*cos(23.4/180.*pi)
    Z  = sin(2*pi*(t-0.22))*sin(23.4/180.*pi)
    Te = 1.0 + 0.0167*sin(2*pi*(t-0.257))
    xm = p[0]*Te*(Y*cos(ra) - X*sin(ra)) + p[1]*(t - t_r)
    ym = p[0]*Te*(Z*cos(dec) - X*cos(ra)*sin(dec) - Y*sin(ra)*sin(dec)) + p[2]*(t - t_r)
    return np.array([xm, ym])

def printpeakmag(src,tnow):
    trend={1:'\033[1;34m ^\033[0;0m',-1:'\033[1;31m v\033[0;0m'}
    c = SkyCoord(get_maser(src).ra,get_maser(src).dec,unit=(u.hourangle,u.deg))
    name = get_maser(src).name
    current  = parpropmot([1.0,0,0],tnow,0.0         ,c.ra.rad,c.dec.rad) 
    nextweek = parpropmot([1.0,0,0],tnow+7/365.25,0.0,c.ra.rad,c.dec.rad)

    if abs(current)[0]>0.6: end1 = '*'
    else:                   end1 = ' '

    if abs(current)[1]>0.6: end2 = '*'
    else:                   end2 = ' '
    ratrend = trend[np.sign(nextweek-current)[0]]+end1
    dctrend = trend[np.sign(nextweek-current)[1]]+end2

    text1 = 'SOURCE {0:3s}|{1:s} \033[0;0m RA {5:+3.2f}{3:3s} DEC {6:+3.2f}{4:3s}'.format(src,name,tnow*365.25,ratrend,dctrend,*current)
    print(text1)

###############################################################################

def main():
	dt = 0.0*u.week
	tnow = (aTime.now()+dt).decimalyear%1
        mzrs = ['s1','s3','s4','s5','s6','s7','s8','s9','s10','s11','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24']
        #mzrs = ['s3','s4','s7','s11']
	print('DAY {:<4.0f}'.format(tnow*365.25))
	for src in mzrs:
	    dt = 0
	    tnow = (aTime.now()+dt*u.day).decimalyear%1
	    printpeakmag(src,tnow)
	print('\n')

###############################################################################

if __name__=='__main__':
	main()

###############################################################################
