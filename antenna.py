import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
from astropy.units import Quantity

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



# RT32 location (Ventspils, Latvia)

latitude = 57.5535171694 
longitude = 21.8545525000
elevation = 20  # Elevation in meters (optional)


# Define constants
path = ''
# year = 2024
# month = 2
# day = 15
# hour_start = 11
# minute_start = 25

az_offset0 = +0.01
el_offset0 = +0.18

az_offset2 = 0.
el_offset2 = 0.

num_scan = 5
step1 = 6.  # arcmin
step2 = 6.
step3 = 6.
step4 = 10.
step5 = 12.
sky = 16. * 4  # arcmin

t_cal = 60.  # seconds
t1 = 40.
t2 = 100.
t3 = 120.
t4 = 120.
t5 = 120.
t_slew = 20.
t_scan = t_cal + t1 + t2 + t3 + t4 + t5 + t_slew + t_cal + t_slew  # duration of scan, seconds

t_spirals = np.array([t1,t2,t3,t4,t5])

times = np.array([t_cal, t1, t2, t3, t4, t5, t_slew, t_cal, t_slew])
labels = np.array(['t_cal', 't1', 't2', 't3', 't4', 't5','t_slew','t_cal_2','t_slew_2'])
times_scan_cumsum = np.cumsum(times)

times_sec_dic = dict(zip(labels, times_scan_cumsum))


def create_antenna_file(year,month,day,hour_start,minute_start, pressure, temperature, relative_humidity,obswl):
    
    # Create EarthLocation object with the given location
    location = EarthLocation(lat=latitude, lon=longitude, height=elevation)
    
    pressure_quantity = Quantity(pressure, unit='hPa')
    temperature_quantity = Quantity(temperature, unit='deg_C')
    humidity_quantity = Quantity(relative_humidity)
    obswl_quantity = Quantity(obswl, unit='um')  


    # Calculate start time
    start_time = Time(f"{year}-{month:02d}-{day:02d} {hour_start:02d}:{minute_start:02d}:00")
    t_start_seconds = time_to_seconds(start_time.datetime)

    sun_location = get_sun(start_time)
    az_start, el_start = sun_location.transform_to(AltAz(obstime=start_time, location=location , pressure=pressure_quantity ,temperature=temperature_quantity , relative_humidity=humidity_quantity , obswl=obswl_quantity )).az, \
                            sun_location.transform_to(AltAz(obstime=start_time, location=location , pressure=pressure_quantity ,temperature=temperature_quantity , relative_humidity=humidity_quantity , obswl=obswl_quantity )).alt

    file_out = f"{(year - 2000):02d}{month:02d}{day:02d}_{hour_start:02d}{minute_start:02d}"
    file_name1 = f"sun_scan_{file_out}.ptf"

    with open(path + file_name1, 'w') as file1:
        file1.write("# Table of horizontal coordinates for Sun scanning\n")
        file1.write("# Celestial Object ID: SUN\n")
        file1.write("# Spiral scanning, number of scans: " + str(num_scan) + "\n")
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

        
        file1.write(f"# Az@Start Time      : {az_start:.5f} deg.\n")
        file1.write(f"# El@Start Time      : {el_start:.5f} deg.\n")
        file1.write(f"# Start Time (UTC+0) : {start_time.iso}\n")

        utime_start = start_time.jd  # start time in JD
        time_session = num_scan * t_scan / 3600. / 24  # session duration in hours
        utime_mean = utime_start + time_session / 2.  # mean UT time of session
        utime_end = utime_start + time_session  # end UT time of session
        
        file1.write(f"# End Time, (UTC), hours: {Time(utime_end, format='jd').iso}\n")
        file1.write(f"# Mean observation time (UTC): {Time(utime_mean, format='jd').iso}\n")
        file1.write(f"# Az offset          : {az_offset0} deg.\n")
        file1.write(f"# El offset          : {el_offset0} deg.\n")
        file1.write("#!!! Offset of the antenna pointing system must be set 0 !!!\n")
        file1.write("\n[Interpolation]\n")
        file1.write("Newton\n")
        file1.write("\n[Load Mode]\n")
        file1.write("New\n")
        file1.write("\n[Start Time]\n")
        file1.write(f"{start_time.iso}\n")
        file1.write("\n[Table Data]\n")
        file1.write("#  Time   Az source [degr.]     El source[degr.]\n")

        # Calculate coordinates
        x = np.zeros(700)
        y = np.zeros(700)

        # Calibration Sun center

        for i in range(int(t_cal)):
            x[i] = 0.
            y[i] = 0.

        # 1 turn
        for i in range(int(t1)):
            i0 = int(t_cal) 
            fi = i * 360. / t1
            r = i * step1 / t1
            x[int(i + i0)] = r * np.cos(np.deg2rad(fi))
            y[int(i + i0)] = r * np.sin(np.deg2rad(fi))

        # 2 turn
        for i in range(int(t2)):
            i0 = int(t_cal) + int(t1)
            fi = i * 360. / t2
            r = step1 + i * (step2 / t2)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # 3 turn
        for i in range(int(t3)):
            i0 = int(t_cal) + int(t1) + int(t2)
            fi = i * 360. / t3
            r = step1 + step2 + i * (step3 / t3)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # 4 turn
        for i in range(int(t4)):
            i0 = int(t_cal) + int(t1) + int(t2) + int(t3)
            fi = i * 360. / t4
            r = step1 + step3 + step3 + (i * step4 / t4)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # 5 turn
        for i in range(int(t5)):
            i0 = int(t_cal) + int(t1) + int(t2) + int(t3) + int(t4)
            fi = i * 360. / t5
            r = step1 + step2 + step3 + step4 + i * (step5 / t5)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # slew to clibration, Sky
        for i in range(int(t_slew)):
            i0 = int(t_cal) + int(t1) + int(t2) + int(t3) + int(t4) + int(t5)
            y[i + i0] = 0.
            x0 = step1 + step2 + step3 + step4 + step5
            x[i + i0] = x0 + i * (sky - x0) / t_slew

        # calibration sky
        for i in range(int(t_cal)):
            i0 = int(t_cal) + int(t1) + int(t2) + int(t3) + int(t4) + int(t5) + int(t_slew)
            x[i + i0] = sky
            y[i + i0] = 0.

        # Seventh loop: slew to Sun center
        for i in range(int(t_slew)):
            i0 = int(t_cal) + int(t1) + int(t2) + int(t3) + int(t4) + int(t5) + int(t_slew) + int(t_cal)
            y[i + i0] = 0.
            x[i + i0] = (int(t_slew) - i - 1) * sky / t_slew

            xx = np.zeros(num_scan * int(t_scan))
            yy = np.zeros(num_scan * int(t_scan))

        for j in range(num_scan):
            ff = j * (360.0 / num_scan)
            for i in range(int(t_scan)):
                ii = j * int(t_scan) + i
                xx[ii] = x[i] * np.cos(np.deg2rad(ff)) - y[i] * np.sin(np.deg2rad(ff))
                yy[ii] = x[i] * np.sin(np.deg2rad(ff)) + y[i] * np.cos(np.deg2rad(ff))
                
        
        utc = utime_start + np.arange(num_scan * t_scan) / 3600. / 24   
        q = RT32_SUN_PARA(utc , location)

        xx1 = xx * np.cos(np.deg2rad(q)) - yy * np.sin(np.deg2rad(q))
        yy1 = xx * np.sin(np.deg2rad(q)) + yy * np.cos(np.deg2rad(q))


        az_sun = sun_location.transform_to(AltAz(obstime=Time(utc, format='jd'), location=location ,pressure=pressure_quantity ,temperature=temperature_quantity , relative_humidity=humidity_quantity , obswl=obswl_quantity)).az.deg
        el_sun = sun_location.transform_to(AltAz(obstime=Time(utc, format='jd'), location=location, pressure=pressure_quantity ,temperature=temperature_quantity , relative_humidity=humidity_quantity , obswl=obswl_quantity)).alt.deg

        az_anten = az_sun + xx1 / np.cos(np.deg2rad(el_sun)) / 60.
        el_anten = el_sun + yy1 / 60.

    
        for i in range(len(utc)):
            file1.write(f"{Time(utc[i], format='jd').iso}.00 \t\t {az_anten[i]:.5f} \t\t {el_anten[i]:.5f}\n")

    print('-------------------------------------------------------------')
    print('Saved: ', file_name1, "  ",  len(utc), "  points")
    return file_name1
