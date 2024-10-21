import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun ,SkyCoord
from astropy.table import Table
import pandas as pd
import astropy.units as u
from sunpy.coordinates import frames, sun
import matplotlib.pyplot as plt
import sunpy.map
from scipy.interpolate import Rbf
import os

class Weather:
    def __init__(self , temperature, pressure, relative_humidity, obswl):
        self.temperature  = temperature
        self.pressure  = pressure
        self.relative_humidity  = relative_humidity
        self.obswl  = obswl               
    

class Antenna:
    def __init__(self, az_offset0=0.08, el_offset0=-0.17, az_offset2=0.0, el_offset2=0.0, description="Default Antenna"):
        self.az_offset0 = az_offset0
        self.el_offset0 = el_offset0
        self.az_offset2 = az_offset2
        self.el_offset2 = el_offset2
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.location = None
        self.description = description
        
        
    def set_location(self, latitude, longitude, elevation):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.location = EarthLocation(lat=self.latitude, lon=self.longitude, height=self.elevation)

class RT32(Antenna):
    def __init__( description="RT-32 Antenna"):
        super().__init__(description=description)


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
        

def RT32_SUN_PARA(utc , location):
    # Define constants
    lat_rt32 = 57.5535171694  # RT32 geographic latitude in degrees
    lon_rt32 = 21.8545525000  # RT32 geographic longitude in degrees

    # Ensure utc is an array
    if np.isscalar(utc):
        utime = np.array([utc])
    else:
        utime = np.array(utc)

    q = np.zeros(len(utc))

    for i in range(len(utime)):
        # Define observer's location
        observer_location = location

        # Create Time object
        observation_time = Time(utc[i], format='jd', location=observer_location)

        # Calculate solar alt-az coordinates
        sun_altaz = get_sun(observation_time).transform_to(AltAz(obstime=observation_time, location=observer_location))
        
        # Calculate local sidereal time
        sidereal_time = observation_time.sidereal_time('mean')

        # Calculate hour angle of the Sun
        ha = sidereal_time - sun_altaz.az
   
        # Calculate tangent of parallax angle
        q_tan = np.sin(ha) / (np.tan(lat_rt32) * np.cos(sun_altaz.alt) - 
                              np.sin(sun_altaz.alt) * np.cos(ha))

        
        # Parallax angle in degrees
        q[i] = -np.rad2deg(np.arctan(q_tan)).value # Positive values are to the left, negative to the right
        
    
    # Return result
    if len(q) == 1:
        return q[0]
    else:
        return q


def time_to_seconds(time):
    hour = time.hour
    minute = time.minute
    second = time.second
    microsecond = time.microsecond
    decimal_hour = hour + minute / 60 + second / 3600 + microsecond / 3600000000
    return decimal_hour * 3600

def seconds_to_time(year, month, day, seconds):
    hours = int(seconds // 3600)
    remaining_seconds = seconds % 3600
    minutes = int(remaining_seconds // 60)
    remaining_seconds %= 60
    remaining_seconds = int(remaining_seconds)

    time = Time(f"{year}-{month:02d}-{day:02d} {hours}:{minutes}:{remaining_seconds}")
    return time.isot

def bintable_to_pandas(file_path, hdu_number):
    try:
        # Leer la tabla binaria del archivo FITS
        table = Table.read(file_path, hdu=hdu_number)

        # Convertir la tabla en un DataFrame de Pandas
        df = table.to_pandas()
       
        return df
    except Exception as e:
        print("Error:", e)
        return None
    
def bintable_to_pandas_OLD(file_path, hdu_number):
    try:
        # Leer la tabla binaria del archivo FITS
        table = Table.read(file_path, hdu=hdu_number)        
        # Crear una lista para almacenar las filas descomprimidas
        rows = []
        # Iterar sobre cada columna de la tabla
        for colname in table.colnames:
            # Verificar si la columna contiene un array
            if isinstance(table[colname][0], np.ndarray):
                # Iterar sobre cada elemento del array
                for i in range(len(table[colname][0])):
                    # Crear un diccionario para la nueva fila
                    new_row_dict = {}
                    # Iterar sobre cada columna nuevamente para llenar la nueva fila
                    for colname_inner in table.colnames:
                        if isinstance(table[colname_inner][0], np.ndarray):
                            new_row_dict[colname_inner] = table[colname_inner][0][i]
                        else:
                            new_row_dict[colname_inner] = table[colname_inner][0]
                    rows.append(new_row_dict)
            else:
                # Si la columna no contiene un array, agregar la fila original
                row_dict = {colname: table[colname][0] for colname in table.colnames}
                rows.append(row_dict)
            # Convertir la lista de diccionarios en una nueva tabla
            table_descompressed = Table(rows)        
            # Convertir la tabla filtrada en un DataFrame de Pandas
            df = table_descompressed.to_pandas()     

            return df 
        
    except Exception as e:
        print("Error:", e)
        return None
     
def bintable_to_pandas_OLD_BROKEN(file_path, hdu_number):
    try:
        # Leer la tabla binaria del archivo FITS
        table = Table.read(file_path, hdu=hdu_number)     
        # Crear una lista para almacenar las filas descomprimidas
        rows = []

        # Iterar sobre cada fila de la tabla
        for row_index, row in enumerate(table):
            # Crear un diccionario para la fila original
            row_dict = {}
            # Iterar sobre cada columna
            for colname in table.colnames:
                # Verificar si la columna contiene un array y si es la primera fila
                if isinstance(row[colname], np.ndarray) and row_index == 0:
                    # Si es la primera fila y es un array, agregar cada elemento como una fila adicional
                    for i, value in enumerate(row[colname]):
                        new_row_dict = row_dict.copy()  # Copiar el diccionario de la fila original
                        new_row_dict[colname] = value  # Actualizar el valor de la columna con el elemento del array
                        rows.append(new_row_dict)  # Agregar la fila a la lista de filas
                else:
                    # Si no es un array o no es la primera fila, agregar el valor a la fila original
                    row_dict[colname] = row[colname]
            # Agregar la fila original a la lista de filas
            rows.append(row_dict)

        # Convertir la lista de diccionarios en una nueva tabla
        table_descompressed = Table(rows)

        # Convertir la tabla filtrada en un DataFrame de Pandas
        df = table_descompressed.to_pandas()       

        return df
    except Exception as e:
        print("Error:", e)
        return None

def getFinalProcessedData(observation, sunPositionDf, data_dfs):

    # Convertir Julian_Time a objetos Time de astropy
    time = Time(sunPositionDf['UTC'], format='jd')

    # Calcular las horas decimales con precisión de milisegundos para cada valor
    decimal_seconds = [time_to_seconds(datetime_obj) for datetime_obj in time.datetime]

    # Agregar las horas decimales al DataFrame
    sunPositionDf['UTC'] = np.round(np.array(decimal_seconds).astype(float), 3)

    columns = sunPositionDf.columns
    interpolated_values = np.empty((sunPositionDf.shape[0] * 1000, sunPositionDf.shape[1]))
    for i, col in enumerate(columns):
        for j in range(len(sunPositionDf)-1):
            interpolated_values[j*1000:(j+1)*1000, i] = np.linspace(sunPositionDf.iloc[j][col], sunPositionDf.iloc[j+1][col], 1001)[0:-1]

    interpolated_df = pd.DataFrame(interpolated_values, columns=columns)
    pd.set_option('display.float_format', '{:.10f}'.format)

    # Diccionario para almacenar los DataFrames procesados por banda
    band_processed_dfs = {}

    for band, data_df in data_dfs.items():
        data_df[f'UTC_{band}'] = data_df[f'UTC_{band}'].astype('float64').round(3)
        
        interpolated_df['UTC'] = interpolated_df['UTC'].round(3)               

        
        merged_df = pd.merge(interpolated_df, data_df, left_on='UTC', right_on=f'UTC_{band}')

        merged_df['SunX'] = merged_df['SunX'] + observation.antenna.az_offset0
        merged_df['SunY'] = merged_df['SunY'] + observation.antenna.el_offset0

      
        t_scan_cumsum = np.linspace(0, 4, 5) * observation.t_scan + observation.t_start_seconds

        def filter_dataframe(df, start, end):
            return df[(df['UTC'] >= start) & (df['UTC'] <= end)]

        filtered_dfs = {}
        cal_df_centre = []
        cal_df_sky = []

        columns_to_remove = ["UTC", "SunX", "SunY", f'UTC_{band}']       

        for i, start_seconds in enumerate(t_scan_cumsum):
            start_t_cal_1 = start_seconds
            end_t_cal_1 = start_seconds + observation.times_sec_dic['t_cal']

            start_t_slew_1 = start_seconds + observation.times_sec_dic['t5']
            end_t_slew_1 = start_seconds + observation.times_sec_dic['t_slew']

            start_t_cal_2 = start_seconds + observation.times_sec_dic['t_slew']
            end_t_cal_2 = start_seconds + observation.times_sec_dic['t_cal_2']

            start_t_slew_2 = start_seconds + observation.times_sec_dic['t_cal_2']
            end_t_slew_2 = start_seconds + observation.times_sec_dic['t_slew_2']

            filters = [
                (f'filter_it_{i+1}_t_cal_1', start_t_cal_1, end_t_cal_1),
                (f'filter_it_{i+1}_t_slew_2', start_t_slew_1, end_t_slew_1),
                (f'filter_it_{i+1}_t_cal_2', start_t_cal_2, end_t_cal_2),
                (f'filter_it_{i+1}_t_slew_3', start_t_slew_2, end_t_slew_2)
            ]

            for filter_name, start, end in filters:
                if filter_name.endswith("t_cal_1"):
                    cal_df_centre.append(filter_dataframe(merged_df, start, end).drop(columns=columns_to_remove).mean().values)
                elif filter_name.endswith("t_cal_2"):
                    cal_df_sky.append(filter_dataframe(merged_df, start, end).drop(columns=columns_to_remove).mean().values)
                filtered_dfs[filter_name] = filter_dataframe(merged_df, start, end)

        filtered_indices = set().union(*[filtered_df.index for filtered_df in filtered_dfs.values()])
        rest_of_df = merged_df[~merged_df.index.isin(filtered_indices)]

        cal_df_centre = np.array(cal_df_centre)
        cal_df_sky = np.array(cal_df_sky)

        sun_centre_means = np.mean(cal_df_centre, axis=0) 
        sky_means = np.mean(cal_df_sky, axis=0)      

        rest_of_df[f'RCP_{band}'] = (rest_of_df[f'RCP_{band}'] - sky_means[0]) / (sun_centre_means[0] - sky_means[0])
        rest_of_df[f'LCP_{band}'] = (rest_of_df[f'LCP_{band}'] - sky_means[0]) / (sun_centre_means[0] - sky_means[0])

        rest_of_df[f'STOKE_I_{band}'] = (rest_of_df[f'RCP_{band}'] + rest_of_df[f'LCP_{band}']) / 2
        rest_of_df[f'STOKE_V_{band}'] = (rest_of_df[f'RCP_{band}'] - rest_of_df[f'LCP_{band}']) / 2

        rest_of_df["isoT_time"] = rest_of_df.apply(lambda row: seconds_to_time(observation.year, observation.month, observation.day, row["UTC"]), axis=1)
        band_processed_dfs[band] = rest_of_df

    return band_processed_dfs


def rename_columns(data_df):
    # Renombrar las columnas eliminando los dos primeros números y manteniendo los prefijos LCP y RCP
    renamed_columns = {
        col: f"{col.split()[0]} {col.split()[2]}"
        if col.startswith(('LCP', 'RCP'))
        else col
        for col in data_df.columns
    }
    data_df = data_df.rename(columns=renamed_columns)

     # Tabla de mapeo de números a bandas
    band_mapping = {
        '01': '4.07GHZ',
        '04': '6.42GHZ',
        '07': '8.40GHZ',
        '09': '9.80GHZ',
        '11': '11.90GHZ'
    }

    # Renombrar columnas UTC
    renamed_columns = {
        col: f'UTC {col.split()[1]} {band_mapping[col.split()[2]]}' 
        if col.startswith('UTC') and col.split()[2] in list(band_mapping.keys())
        else col
        for col in data_df.columns
    }
    data_df = data_df.rename(columns=renamed_columns)

    return data_df

def processData(data_df):   


    data_df = rename_columns(data_df)
    
    # Lista de bandas a procesar (sin los dos primeros números)
    bands = ['4.07GHZ', '6.42GHZ', '8.40GHZ', '9.80GHZ', '11.90GHZ']

    # Diccionario para almacenar los DataFrames de cada banda
    band_dfs = {}

    for band in bands:
        # Extracción y redondeo de UTC
        UTC_RCP = np.round(data_df[f'UTC RCP {band}'].dropna() * 3600, 3)
        UTC_LCP = np.round(data_df[f'UTC LCP {band}'].dropna() * 3600, 3)
        RCP = data_df[f'RCP {band}'].dropna()
        LCP = data_df[f'LCP {band}'].dropna()

        # STOKE_I = (RCP.values + LCP.values) / 2 
        # STOKE_V = ( RCP.values - LCP.values) / 2 

        # print(STOKE_I_11_4_11_90GHZ.size)
        # print(STOKE_V_11_4_11_90GHZ.size)
        # print(UTC_RCP_11.size)      

        # Crear DataFrame para la banda actual
        band_df = pd.DataFrame({
            f'UTC_{band}': UTC_RCP.values,          
            f'RCP_{band}': RCP,
            f'LCP_{band}': LCP
        })

        # Almacenar en el diccionario
        band_dfs[band] = band_df

    return band_dfs



def process_all_heliocentric_coordinates(band_processed_dfs, observation):
    def process_heliocentric_coordinates(band_df, observation):       

        coordsXHelio = []
        coordsYHelio = []

        # Convert the array of times into a list of Time objects
        times = [Time(t) for t in band_df['isoT_time']]

        # Calculate the sun's positions for each time in the list
        az_sun = observation.sun_location.transform_to(AltAz(obstime=times, location=observation.antenna.location, pressure=observation.weather.pressure , temperature=observation.weather.temperature, relative_humidity=observation.weather.relative_humidity ,obswl=observation.weather.obswl)).az.deg  + observation.antenna.az_offset0
        el_sun = observation.sun_location.transform_to(AltAz(obstime=times, location=observation.antenna.location, pressure=observation.weather.pressure , temperature=observation.weather.temperature, relative_humidity=observation.weather.relative_humidity ,obswl=observation.weather.obswl)).alt.deg + observation.antenna.el_offset0

        az_anten = az_sun + band_df['SunX'] / np.cos(np.deg2rad(el_sun)) / 60.
        el_anten = el_sun + band_df['SunY'] / 60.    

        band_df['az_anten'] = az_anten
        band_df['el_anten'] = el_anten

        # Iterar sobre los diferentes momentos de tiempo
        for index, row in band_df.iterrows():
            # Convert the AltAz coordinates from the DataFrame from degrees to radians
            el_deg = row['el_anten'] * u.deg
            az_deg = row['az_anten'] * u.deg

            # Set the observation time
            obstime = row['isoT_time']

            # Convertir a coordenadas heliocéntricas
            frame_altaz = AltAz(obstime=Time(obstime), location=observation.antenna.location, pressure=observation.weather.pressure , temperature=observation.weather.temperature, relative_humidity=observation.weather.relative_humidity ,obswl=observation.weather.obswl)
            sun_helio = SkyCoord(alt=el_deg, az=az_deg, observer='earth', distance=sun.earth_distance(obstime), frame=frame_altaz).transform_to(frames.Helioprojective)

            
            # Append the transformed coordinates to the list
            coordsXHelio.append(sun_helio.Tx.value)
            coordsYHelio.append(sun_helio.Ty.value)

        band_df['tx_helio_anten'] = coordsXHelio
        band_df['ty_helio_anten'] = coordsYHelio

        return band_df
    
    # Aplicar la función de procesamiento a cada DataFrame en el diccionario
    for band, df in band_processed_dfs.items():
        band_processed_dfs[band] = process_heliocentric_coordinates(df, observation)
    
    return band_processed_dfs