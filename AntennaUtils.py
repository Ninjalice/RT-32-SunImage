import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
from astropy.table import Table
from AntennaUtils import *


class Antenna:
    def __init__(self, az_offset0=+0.01, el_offset0=+0.18, az_offset2=0.0, el_offset2=0.0, description="Default Antenna"):
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
    def __init__(self , antenna_instance : Antenna,year,month,day,hour_start,minute_start):

        self.antenna = antenna_instance

        self.year = year
        self.month = month
        self.day = day
        self.hour_start = hour_start
        self.minute_start = minute_start

        self.start_time = Time(f"{year}-{month:02d}-{day:02d} {hour_start:02d}:{minute_start:02d}:00")
        self.t_start_seconds = time_to_seconds(self.start_time.datetime)

        self.sun_location = get_sun(self.start_time)
        self.az_start, self.el_start = self.sun_location.transform_to(AltAz(obstime=self.start_time, location=self.antenna.location)).az, \
                                self.sun_location.transform_to(AltAz(obstime=self.start_time, location=self.antenna.location)).alt

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

        
        az_sun = self.sun_location.transform_to(AltAz(obstime=Time(utc, format='jd'), location=self.antenna.location)).az.deg
        el_sun = self.sun_location.transform_to(AltAz(obstime=Time(utc, format='jd'), location=self.antenna.location)).alt.deg

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
            file1.write(f"# Start Time (UTC+0) : {self.start_time.iso}\n")        
            
            file1.write(f"# End Time, (UTC), hours: {Time(self.utime_end, format='jd').iso}\n")
            file1.write(f"# Mean observation time (UTC): {Time(self.utime_mean, format='jd').iso}\n")
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
                file1.write(f"{Time(utc[i], format='jd').iso}.00 \t\t {az_anten[i]:.5f} \t\t {el_anten[i]:.5f}\n")

        print('-------------------------------------------------------------')
        print('Saved: ', self.file_name1, "  ",  len(utc), "  points")
        return True


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


def bintable_to_pandas(file_path, hdu_number):
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
    
