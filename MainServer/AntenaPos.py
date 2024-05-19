
import os
import datetime
import re
import requests
import openmeteo_requests
import requests_cache
from retry_requests import retry
from supabase import create_client
from urllib.parse import quote



supabase_url = "https://rcwyiclyrlcrpyvpdrza.supabase.co"
supabase_key = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InJjd3lpY2x5cmxjcnB5dnBkcnphIiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcxMjczNDIwNCwiZXhwIjoyMDI4MzEwMjA0fQ._T5lYpWelk_iPojWZwmc3e7DXE9BxbmKzfB6c8kaKkc"

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



def weatherDataBase():
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
    current = response.Current()
    return current

def get_latest_ptf_id():
    ptf_files = [f for f in os.listdir('.') if f.endswith('.ptf')]
    ptf_files.sort(key=lambda x: os.path.getmtime(os.path.join('.', x)))
    
    if ptf_files:
        latest_ptf = ptf_files[-1]
        
        match = re.search(r"sun_scan_(\d{6})_(\d{4})", latest_ptf)
        if match:
            date_part = match.group(1)
            time_part = match.group(2)
            id = date_part + time_part
            return id
        else:
            return "No match found in filename."
    else:
        return "No PTF files found."
    
def upload_images_to_bucket(directory):
    observation_name = get_latest_ptf_id()
    if observation_name:
        observation_folder = f"{observation_name}"
        image_files = [f for f in os.listdir(directory) if f.endswith('.jpg') or f.endswith('.png')]
        sb = create_client(supabase_url, supabase_key)
        for image_file in image_files:
            file_path = os.path.join(directory, image_file)
            with open(file_path, 'rb') as file:
                response = sb.storage.from_("images").upload(f"{observation_folder}/{image_file}", file)
            os.remove(file_path)
    else:
        print("No observation name found.")
        

def get_latest_ptf_id():
    ptf_files = [f for f in os.listdir('.') if f.endswith('.ptf')]
    ptf_files.sort(key=lambda x: os.path.getmtime(os.path.join('.', x)))
    
    if ptf_files:
        latest_ptf = ptf_files[-1]
        
        match = re.search(r"sun_scan_(\d{6})_(\d{4})", latest_ptf)
        if match:
            date_part = match.group(1)
            time_part = match.group(2)
            id = date_part + time_part
            return id
        else:
            return None
    else:
        return None

import os

def download_images_from_folder():
    sb = create_client(supabase_url, supabase_key)
    folder_name = get_latest_ptf_id()
    bucket_name = "images"
    num_images = 10
    download_dir = "temp_images"
    os.makedirs(download_dir, exist_ok=True)
    for i in range(1, num_images + 1):
        file_name = f"{folder_name}/{folder_name}_{i:02}.jpg"
        destination = os.path.join(download_dir, f"{folder_name}_{i:02}.jpg")
        with open(destination, 'wb+') as f:
            res = sb.storage.from_(bucket_name).download(file_name)
            f.write(res)
    return download_dir

def insert_image_to_table():
    sb = create_client(supabase_url, supabase_key)
    observation_id = get_latest_ptf_id()
    for file_name in os.listdir("temp_images"):
        url = f"https://{supabase_url}/storage/v1/object/public/images/{observation_id}/{quote(file_name)}"
        sb.table("images").insert({"observationid": observation_id, "url": url}).execute()
    filelist = [f for f in os.listdir("temp_images")]
    for f in filelist:
        os.remove(os.path.join("temp_images", f))
    os.rmdir("temp_images")
        
def send_data_to_database():
    weather_data = weatherDataBase()
    
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

    current_date = datetime.datetime.now()
    formatted_date = current_date.strftime('%Y-%m-%d')

    data_to_insert = {
        'id': str(get_latest_ptf_id()),
        'weather': weather_data_json,
        'date': formatted_date
    }

    response = requests.post(
        f"{supabase_url}/rest/v1/testing",
        headers={
            "apikey": supabase_key,
            "Content-Type": "application/json"
        },
        json=data_to_insert
    )

