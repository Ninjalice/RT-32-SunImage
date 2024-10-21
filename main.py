import pandas as pd
import astropy.units as u
from AntennaUtils import *  
import warnings

# Imprimimos estadísticas resumidas del DataFrame final
pd.set_option('display.float_format', '{:.10f}'.format)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)

# RT32 location (Ventspils, Latvia)
rt32_antenna = RT32()
rt32_antenna.set_location(latitude=57.5535171694, longitude=21.8545525000, elevation=20)

# Define constants
year = 2024
month = 4
day = 25
hour_start = 8
minute_start = 55

temperature = u.Quantity(20.0, unit=u.deg_C)
pressure = u.Quantity(1013.25, unit=u.hPa)
relative_humidity = u.Quantity(60.0, unit=u.percent)
obswl =u.Quantity(50000, unit=u.nm) 

weather = Weather(temperature, pressure, relative_humidity, obswl)

observation = SpiralSunObservation(weather,rt32_antenna , year , month , day , hour_start , minute_start)

fit_file_path = "lnsp4_5ch_241020_173655_241020_173725.fit"

observation.createImages(fit_file_path)