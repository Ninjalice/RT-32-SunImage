from astropy.coordinates import EarthLocation

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
