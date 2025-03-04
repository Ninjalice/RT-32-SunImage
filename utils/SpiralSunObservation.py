import numpy as np
from astropy.time import Time
from astropy.coordinates import  AltAz, get_sun 
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from scipy.interpolate import Rbf
import os

from models import Weather, Antenna
from utils.mapFunctions import RT32_SUN_PARA, bintable_to_pandas, getFinalProcessedData, process_all_heliocentric_coordinates, processData, time_to_seconds


class SpiralSunObservation:
    def __init__(self , weather: Weather , antenna_instance : Antenna,year,month,day,hour_start,minute_start):
        self.weather = weather
        self.antenna = antenna_instance

        self.year = year
        self.month = month
        self.day = day
        self.hour_start = hour_start
        self.minute_start = minute_start

        self.start_time = Time(f"{year}-{month:02d}-{day:02d} {hour_start:02d}:{minute_start:02d}:00")
        self.t_start_seconds = time_to_seconds(self.start_time.datetime)

        self.sun_location = get_sun(self.start_time)
        self.az_start, self.el_start = self.sun_location.transform_to(AltAz(obstime=self.start_time, location=self.antenna.location, pressure=self.weather.pressure , temperature=self.weather.temperature, relative_humidity=self.weather.relative_humidity ,obswl=self.weather.obswl )).az, \
                                self.sun_location.transform_to(AltAz(obstime=self.start_time, location=self.antenna.location, pressure=self.weather.pressure , temperature=self.weather.temperature, relative_humidity=self.weather.relative_humidity ,obswl=self.weather.obswl)).alt

        self.file_out = f"{(year - 2000):02d}{month:02d}{day:02d}_{hour_start:02d}{minute_start:02d}"
        self.file_name1 = f"sun_scan_{self.file_out}.ptf"       

        self.num_scan = 5
        self.step1 = 6.  # arcmin
        self.step2 = 6.
        self.step3 = 6.
        self.step4 = 10.
        self.step5 = 12.
        self.sky = 16. * 4  # arcmin

        self.t_cal = 60.  # seconds
        self.t1 = 40.
        self.t2 = 100.
        self.t3 = 120.
        self.t4 = 120.
        self.t5 = 120.
        self.t_slew = 20.
        self.t_scan = self.t_cal + self.t1 + self.t2 + self.t3 + self.t4 + self.t5 + self.t_slew + self.t_cal + self.t_slew  # duration of scan, seconds

        self.t_spirals = np.array([self.t1, self.t2, self.t3, self.t4, self.t5])

        self.times = np.array([self.t_cal, self.t1, self.t2, self.t3, self.t4, self.t5, self.t_slew, self.t_cal, self.t_slew])
        self.labels = np.array(['t_cal', 't1', 't2', 't3', 't4', 't5', 't_slew', 't_cal_2', 't_slew_2'])
        self.times_scan_cumsum = np.cumsum(self.times)

        self.times_sec_dic = dict(zip(self.labels, self.times_scan_cumsum))

        self.utime_start = self.start_time.jd  # start time in JD
        self.time_session = self.num_scan * self.t_scan / 3600. / 24  # session duration in hours
        self.utime_mean = self.utime_start + self.time_session / 2.  # mean UT time of session
        self.utime_end = self.utime_start + self.time_session  # end UT time of session
        

    def calculatePositions(self):
        # Calculate coordinates
        x = np.zeros(700)
        y = np.zeros(700)

        # Calibration Sun center
        
        for i in range(int(self.t_cal)):
            x[i] = 0.
            y[i] = 0.

        # 1 turn
        for i in range(int(self.t1)):
            i0 = int(self.t_cal) 
            fi = i * 360. / self.t1
            r = i * self.step1 / self.t1
            x[int(i + i0)] = r * np.cos(np.deg2rad(fi))
            y[int(i + i0)] = r * np.sin(np.deg2rad(fi))

        # 2 turn
        for i in range(int(self.t2)):
            i0 = int(self.t_cal) + int(self.t1)
            fi = i * 360. / self.t2
            r = self.step1 + i * (self.step2 / self.t2)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # 3 turn
        for i in range(int(self.t3)):
            i0 = int(self.t_cal) + int(self.t1) + int(self.t2)
            fi = i * 360. / self.t3
            r = self.step1 + self.step2 + i * (self.step3 / self.t3)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # 4 turn
        for i in range(int(self.t4)):
            i0 = int(self.t_cal) + int(self.t1) + int(self.t2) + int(self.t3)
            fi = i * 360. / self.t4
            r = self.step1 + self.step3 + self.step3 + (i * self.step4 / self.t4)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # 5 turn
        for i in range(int(self.t5)):
            i0 = int(self.t_cal) + int(self.t1) + int(self.t2) + int(self.t3) + int(self.t4)
            fi = i * 360. / self.t5
            r = self.step1 + self.step2 + self.step3 + self.step4 + i * (self.step5 / self.t5)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # slew to clibration, Sky
        for i in range(int(self.t_slew)):
            i0 = int(self.t_cal) + int(self.t1) + int(self.t2) + int(self.t3) + int(self.t4) + int(self.t5)
            y[i + i0] = 0.
            x0 = self.step1 + self.step2 + self.step3 + self.step4 + self.step5
            x[i + i0] = x0 + i * (self.sky - x0) / self.t_slew

        # calibration sky
        for i in range(int(self.t_cal)):
            i0 = int(self.t_cal) + int(self.t1) + int(self.t2) + int(self.t3) + int(self.t4) + int(self.t5) + int(self.t_slew)
            x[i + i0] = self.sky
            y[i + i0] = 0.

        # Seventh loop: slew to Sun center
        for i in range(int(self.t_slew)):
            i0 = int(self.t_cal) + int(self.t1) + int(self.t2) + int(self.t3) + int(self.t4) + int(self.t5) + int(self.t_slew) + int(self.t_cal)
            y[i + i0] = 0.
            x[i + i0] = (int(self.t_slew) - i - 1) * self.sky / self.t_slew

            xx = np.zeros(self.num_scan * int(self.t_scan))
            yy = np.zeros(self.num_scan * int(self.t_scan))

        for j in range(self.num_scan):
            ff = j * (360.0 / self.num_scan)
            for i in range(int(self.t_scan)):
                ii = j * int(self.t_scan) + i
                xx[ii] = x[i] * np.cos(np.deg2rad(ff)) - y[i] * np.sin(np.deg2rad(ff))
                yy[ii] = x[i] * np.sin(np.deg2rad(ff)) + y[i] * np.cos(np.deg2rad(ff))
                
        
        utc = self.utime_start + np.arange(self.num_scan * self.t_scan) / 3600. / 24   
        q = RT32_SUN_PARA(utc , self.antenna.location)

        xx1 = xx * np.cos(np.deg2rad(q)) - yy * np.sin(np.deg2rad(q))
        yy1 = xx * np.sin(np.deg2rad(q)) + yy * np.cos(np.deg2rad(q))

        
        az_sun = self.sun_location.transform_to(AltAz(obstime=Time(utc, format='jd'), location=self.antenna.location, pressure=self.weather.pressure , temperature=self.weather.temperature, relative_humidity=self.weather.relative_humidity ,obswl=self.weather.obswl)).az.deg
        el_sun = self.sun_location.transform_to(AltAz(obstime=Time(utc, format='jd'), location=self.antenna.location, pressure=self.weather.pressure , temperature=self.weather.temperature, relative_humidity=self.weather.relative_humidity ,obswl=self.weather.obswl)).alt.deg

        az_anten = az_sun + xx1 / np.cos(np.deg2rad(el_sun)) / 60.
        el_anten = el_sun + yy1 / 60.    
       

        return az_anten, el_anten , az_sun , el_sun , xx1 , yy1, utc 

    def generateFile(self, path  , az_anten , el_anten , utc):
        with open(path + self.file_name1, 'w') as file1:
            file1.write("# Table of horizontal coordinates for Sun scanning\n")
            file1.write("# Celestial Object ID: SUN\n")
            file1.write("# Spiral scanning, number of scans: " + str(self.num_scan) + "\n")
            file1.write("# Scan time schedule:\n")
            file1.write("# Sun center calibration      60 sec\n")
            file1.write("# 1 spiral turn              40 sec  step 0->6 arcmin\n")
            file1.write("# 2 spiral                   100 sec  step 6->12 arcmin\n")
            file1.write("# 3 spiral                   120 sec step 12->18 arcmin\n")
            file1.write("# 4 spiral                   120 sec step 18-28 arcmin\n")
            file1.write("# 5 spiral                   120 sec step 28-40 arcmin\n")
            file1.write("# Slew to sky 4 Rsun        20 sec\n")
            file1.write("# Sky calibration           60 sec\n")
            file1.write("# Slew to Sun center        20 sec\n")

            
            file1.write(f"# Az@Start Time      : {self.az_start:.5f} deg.\n")
            file1.write(f"# El@Start Time      : {self.el_start:.5f} deg.\n")
            file1.write(f"# Start Time (UTC+0) : {self.start_time.isot}\n")        
            
            file1.write(f"# End Time, (UTC), hours: {Time(self.utime_end, format='jd').isot}\n")
            file1.write(f"# Mean observation time (UTC): {Time(self.utime_mean, format='jd').isot}\n")
            file1.write(f"# Az offset          : {self.antenna.az_offset0} deg.\n")
            file1.write(f"# El offset          : {self.antenna.el_offset0} deg.\n")
            file1.write("#!!! Offset of the antenna pointing system must be set 0 !!!\n")
            file1.write("\n[Interpolation]\n")
            file1.write("Newton\n")
            file1.write("\n[Load Mode]\n")
            file1.write("New\n")
            file1.write("\n[Start Time]\n")
            file1.write(f"{self.start_time.iso}\n")
            file1.write("\n[Table Data]\n")
            file1.write("#  Time   Az source [degr.]     El source[degr.]\n")        

        
            for i in range(len(utc)):
                file1.write(f"{Time(utc[i], format='jd').isot} \t\t {az_anten[i]:.5f} \t\t {el_anten[i]:.5f}\n")

        print('-------------------------------------------------------------')
        print('Saved: ', self.file_name1, "  ",  len(utc), "  points")
        return True

    def createImages(self,fit_file_path):

        #CONSTANTS
        hdu_number = 1  # Number of the extension containing the binary table
        path = ''

        #DATA CALCULATION:
        az_anten, el_anten , az_sun , el_sun , xx1 , yy1, utc = self.calculatePositions()
        self.generateFile(path, az_anten , el_anten , utc)  

        print('-------------------------------------------------------------')
        print('Start processing the FITS file:', fit_file_path)
        # Converts the binary table to a Pandas DataFrame
        data_df = bintable_to_pandas(fit_file_path, hdu_number)        
  
        print('Processing the data...')
        band_data_dfs = processData(data_df)

        sunPositionDf = pd.DataFrame({'UTC': utc,'SunX': xx1, 'SunY': yy1  })

        print('Getting final processed data...')
        processed_dfs = getFinalProcessedData(self , sunPositionDf,band_data_dfs)

        print('Processing heliocentric coordinates...')
        band_processed_helio_dfs = process_all_heliocentric_coordinates(processed_dfs, self )

        # Define the target directory
        directory = f'{self.year}-{self.month}-{self.day}T{self.hour_start}_{self.minute_start}_00'
        
        # Create the directory if it doesn't exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        bands = ['4.07GHZ', '6.42GHZ', '8.40GHZ', '9.80GHZ', '11.90GHZ']

        print('Creating solar maps...')
        for band in bands:
            SunX = band_processed_helio_dfs[band]['tx_helio_anten']
            SunY = band_processed_helio_dfs[band]['ty_helio_anten']
            STOKE_I = band_processed_helio_dfs[band][f'STOKE_I_{band}'].values
            STOKE_V = band_processed_helio_dfs[band][f'STOKE_V_{band}'].values

            # Define the grid covering the helioprojective coordinate space
            tx_min, tx_max = -1200, 1200
            ty_min, ty_max = -1200, 1200
            grid_step = 10  # Adjust as needed

            # Create a grid
            tx, ty = np.meshgrid(np.arange(tx_min, tx_max, grid_step),
                                np.arange(ty_min, ty_max, grid_step))

            # Interpolate power values for each point on the grid using Rbf
            rbf = Rbf(SunX, SunY, STOKE_I, function='linear')
            interp_power_STOKE_I = rbf(tx, ty)

            # Create a grid
            tx, ty = np.meshgrid(np.arange(tx_min, tx_max, grid_step),
                                np.arange(ty_min, ty_max, grid_step))

            # Interpolate power values for each point on the grid using Rbf
            rbf = Rbf(SunX, SunY, STOKE_V, function='linear')
            interp_power_STOKE_V = rbf(tx, ty)    


            # Define metadata for the solar map
            metadata = {
                'date-obs': f'{self.year}-{self.month}-{self.day}T{self.hour_start}:{self.minute_start}:00',  # Adjust this to the correct observation date
                'crval1': 0,
                'crval2': 0,
                'cdelt1': grid_step,
                'cdelt2': grid_step,
                'cunit1': 'arcsec',
                'cunit2': 'arcsec',
                'ctype1': 'HPLN-TAN',
                'ctype2': 'HPLT-TAN',
                'crpix1': (tx_max - tx_min) / (2 * grid_step),
                'crpix2': (ty_max - ty_min) / (2 * grid_step),
                'waveunit': 'm',
                'wavelnth': 0.0262897 * u.m,
                'obsrvtry': 'Ventspils International Radio Astronomy Center',
                'detector': 'LNSP4',
                'dsun_ref': 149597870691,
                'dsun_obs': 151846026489 ,
                'rsun': 1573.89688496,
                'rsun_ref': 696000000 ,                           
                'hglt_obs': 0 * u.deg,
                'hgln_obs': 0 * u.deg,              
                
            }

            # Create a map using the interpolated power values and metadata STOKE I
            interpolated_map = sunpy.map.Map((interp_power_STOKE_I, metadata))

            # Plot the interpolated map using a heatmap with the 'hot' colormap
            plt.ioff()
            plt.rcParams['text.color'] = 'white'
            plt.rcParams['axes.labelcolor'] = 'white'
            plt.rcParams['xtick.color'] = 'white'
            plt.rcParams['ytick.color'] = 'white'   
            fig = plt.figure(figsize=(8, 8))
            fig.patch.set_facecolor('black')            
            ax = fig.add_subplot(projection=interpolated_map)           
            interpolated_map.plot(axes=ax,cmap='gist_heat')
            interpolated_map.draw_limb(axes=ax)
            interpolated_map.draw_grid(axes=ax)
            plt.colorbar(label='Power')
            plt.title(f'STOKE I | {band} {self.year}-{self.month:02d}-{self.day:02d}T{self.hour_start:02d}:{self.minute_start:02d}')
            plt.xlabel('X (arcsec)')
            plt.ylabel('Y (arcsec)')
            circle = plt.Circle((0, 0), 960, color='white', fill=False, linewidth=2)
            ax.add_artist(circle)
            plt.grid(True)   
            name = f'{directory}/LNSP4-{directory}-STOKE_I-{band}.jpeg'
            plt.savefig(name, format='jpeg', dpi=300)     
            plt.close()

            # Create a map using the interpolated power values and metadata STOKE V
            interpolated_map = sunpy.map.Map((interp_power_STOKE_V, metadata))

            # Plot the interpolated map using a heatmap with the 'hot' colormap
            plt.ioff()
            plt.rcParams['text.color'] = 'white'
            plt.rcParams['axes.labelcolor'] = 'white'
            plt.rcParams['xtick.color'] = 'white'
            plt.rcParams['ytick.color'] = 'white'   
            fig = plt.figure(figsize=(8, 8))
            fig.patch.set_facecolor('black')            
            ax = fig.add_subplot(projection=interpolated_map)           
            interpolated_map.plot(axes=ax,cmap='nipy_spectral')
            interpolated_map.draw_limb(axes=ax)
            interpolated_map.draw_grid(axes=ax)         
            plt.colorbar(label='Power')
            plt.title(f'STOKE V | {band} {self.year}-{self.month:02d}-{self.day:02d}T{self.hour_start:02d}:{self.minute_start:02d}')
            plt.xlabel('X (arcsec)')
            plt.ylabel('Y (arcsec)')
            plt.grid(True)
            name = f'{directory}/LNSP4-{directory}-STOKE_V-{band}.jpeg'
            plt.savefig(name, format='jpeg', dpi=300)   
            plt.close()