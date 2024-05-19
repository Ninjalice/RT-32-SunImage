import openmeteo_requests
import requests_cache
from retry_requests import retry


def weatherForecast(option):
    cache_session = requests_cache.CachedSession('.cache', expire_after = 3600)
    retry_session = retry(cache_session, retries = 5, backoff_factor = 0.2)
    openmeteo = openmeteo_requests.Client(session = retry_session)
    url = "https://api.open-meteo.com/v1/forecast"
    if(option == 1):
        params = {
        "latitude": 57.55668258666992,
        "longitude": 21.849212646484375,
        "current": ["temperature_2m", "relative_humidity_2m", "precipitation", "rain", "showers", "snowfall", "cloud_cover", "surface_pressure",                                   "wind_speed_10m", "wind_direction_10m", "wind_gusts_10m"],
        "wind_speed_unit": "ms"
        }
        responses = openmeteo.weather_api(url, params=params)
        response = responses[0]
        current = response.Current()
        current_temperature_2m = current.Variables(0).Value()
        current_relative_humidity_2m = current.Variables(1).Value()
        current_precipitation = current.Variables(2).Value()
        current_rain = current.Variables(3).Value()
        current_showers = current.Variables(4).Value()
        current_snowfall = current.Variables(5).Value()
        current_cloud_cover = current.Variables(6).Value()
        current_surface_pressure = current.Variables(7).Value()
        current_wind_speed_10m = current.Variables(8).Value()
        current_wind_direction_10m = current.Variables(9).Value()
        current_wind_gusts_10m = current.Variables(10).Value()
        print(f"Current time {current.Time()}")
        print(f"Current temperature_2m {current_temperature_2m}")
        print(f"Current relative_humidity_2m {current_relative_humidity_2m}")
        print(f"Current precipitation {current_precipitation}")
        print(f"Current rain {current_rain}")
        print(f"Current showers {current_showers}")
        print(f"Current snowfall {current_snowfall}")
        print(f"Current cloud_cover {current_cloud_cover}")
        print(f"Current surface_pressure {current_surface_pressure}")
        print(f"Current wind_speed_10m {current_wind_speed_10m}")
        print(f"Current wind_direction_10m {current_wind_direction_10m}")
        print(f"Current wind_gusts_10m {current_wind_gusts_10m}")
        
    elif(option == 2):
        params = {
            "latitude": 57.55668258666992,
            "longitude": 21.849212646484375,
            "daily": ["temperature_2m_max", "temperature_2m_min", "daylight_duration", "precipitation_sum", "rain_sum", "snowfall_sum", "precipitation_hours",                       "precipitation_probability_max", "wind_speed_10m_max", "wind_gusts_10m_max", "wind_direction_10m_dominant"],
            "timezone": "auto",
            "forecast_days": 1
        }
        responses = openmeteo.weather_api(url, params=params)
        response = responses[0]
        daily = response.Daily()
        daily_temperature_2m_max = daily.Variables(0).ValuesAsNumpy()
        daily_temperature_2m_min = daily.Variables(1).ValuesAsNumpy()
        daily_daylight_duration = daily.Variables(2).ValuesAsNumpy()
        daily_precipitation_sum = daily.Variables(3).ValuesAsNumpy()
        daily_rain_sum = daily.Variables(4).ValuesAsNumpy()
        daily_snowfall_sum = daily.Variables(5).ValuesAsNumpy()
        daily_precipitation_hours = daily.Variables(6).ValuesAsNumpy()
        daily_precipitation_probability_max = daily.Variables(7).ValuesAsNumpy()
        daily_wind_speed_10m_max = daily.Variables(8).ValuesAsNumpy()
        daily_wind_gusts_10m_max = daily.Variables(9).ValuesAsNumpy()
        daily_wind_direction_10m_dominant = daily.Variables(10).ValuesAsNumpy()
        print(f"Max temperature {daily_temperature_2m_max}")
        print(f"Min temperature {daily_temperature_2m_min}")
        print(f"Daylight duration {daily_daylight_duration}")
        print(f"Precipitation sum {daily_precipitation_sum}")
        print(f"Rain sum {daily_rain_sum}")
        print(f"Snowfall sum {daily_snowfall_sum}")
        print(f"Precipitation hours {daily_precipitation_hours}")
        print(f"Precipitation max probability {daily_precipitation_probability_max}")
        print(f"Max wind speed {daily_wind_speed_10m_max}")
        print(f"Max wind gusts {daily_wind_gusts_10m_max}")
        print(f"Dominant wind direction{daily_wind_direction_10m_dominant}")
        
    else:
        params = {
            "latitude": 57.55668258666992,
            "longitude": 21.849212646484375,
            "daily": ["temperature_2m_max", "temperature_2m_min", "daylight_duration", "precipitation_sum", "rain_sum", "snowfall_sum", "precipitation_hours",                       "precipitation_probability_max", "wind_speed_10m_max", "wind_gusts_10m_max", "wind_direction_10m_dominant"],
            "timezone": "auto",
            "forecast_days": 3
        }
        responses = openmeteo.weather_api(url, params=params)
        response = responses[0]
        daily = response.Daily()
        daily_temperature_2m_max = daily.Variables(0).ValuesAsNumpy()
        daily_temperature_2m_min = daily.Variables(1).ValuesAsNumpy()
        daily_daylight_duration = daily.Variables(2).ValuesAsNumpy()
        daily_precipitation_sum = daily.Variables(3).ValuesAsNumpy()
        daily_rain_sum = daily.Variables(4).ValuesAsNumpy()
        daily_snowfall_sum = daily.Variables(5).ValuesAsNumpy()
        daily_precipitation_hours = daily.Variables(6).ValuesAsNumpy()
        daily_precipitation_probability_max = daily.Variables(7).ValuesAsNumpy()
        daily_wind_speed_10m_max = daily.Variables(8).ValuesAsNumpy()
        daily_wind_gusts_10m_max = daily.Variables(9).ValuesAsNumpy()
        daily_wind_direction_10m_dominant = daily.Variables(10).ValuesAsNumpy()

        print(f"Max temperature {daily_temperature_2m_max}")
        print(f"Min temperature {daily_temperature_2m_min}")
        print(f"Daylight duration {daily_daylight_duration}")
        print(f"Precipitation sum {daily_precipitation_sum}")
        print(f"Rain sum {daily_rain_sum}")
        print(f"Snowfall sum {daily_snowfall_sum}")
        print(f"Precipitation hours {daily_precipitation_hours}")
        print(f"Precipitation max probability {daily_precipitation_probability_max}")
        print(f"Max wind speed {daily_wind_speed_10m_max}")
        print(f"Max wind gusts {daily_wind_gusts_10m_max}")
        print(f"Dominant wind direction{daily_wind_direction_10m_dominant}")



def getCurrentWeather():
    cache_session = requests_cache.CachedSession('.cache', expire_after = 3600)
    retry_session = retry(cache_session, retries = 5, backoff_factor = 0.2)
    openmeteo = openmeteo_requests.Client(session = retry_session)
    url = "https://api.open-meteo.com/v1/forecast"
    params = {
    "latitude": 57.55668258666992,
    "longitude": 21.849212646484375,
    "current": ["temperature_2m", "relative_humidity_2m", "precipitation", "rain", "showers", "snowfall", "cloud_cover", "surface_pressure", "wind_speed_10m", "wind_direction_10m", "wind_gusts_10m"],
    "wind_speed_unit": "ms"
    }
    responses = openmeteo.weather_api(url, params=params)
    response = responses[0]
    weather_data = response.Current()

    weather_data_json = {
        'temperature_2m': weather_data.Variables(0).Value(),
        'relative_humidity_2m': weather_data.Variables(1).Value(),
        'precipitation': weather_data.Variables(2).Value(),
        'rain': weather_data.Variables(3).Value(),
        'showers': weather_data.Variables(4).Value(),
        'snowfall': weather_data.Variables(5).Value(),
        'cloud_cover': weather_data.Variables(6).Value(),
        'surface_pressure': weather_data.Variables(7).Value(),
        'wind_speed_10m': weather_data.Variables(8).Value(),
        'wind_direction_10m': weather_data.Variables(9).Value(),
        'wind_gusts_10m': weather_data.Variables(10).Value()
    }

    return weather_data_json