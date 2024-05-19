import socket
import struct
import threading
import queue
import time
from datetime import datetime
from flask import Flask, render_template, request, jsonify
from flask_socketio import SocketIO, emit
import json
from WeatherUtils import getCurrentWeather
from AntennaUtils import *
import astropy.units as u

app = Flask(__name__)
socketio = SocketIO(app)

# Tamaño del buffer de la cola
BUFFER_SIZE = 2048

# Crear una cola para el buffer
buffer_queue = queue.Queue(BUFFER_SIZE)

# Variables globales para controlar los hilos
receiver_thread = None
processor_thread = None
monitor_thread = None
running = False

output_file = None

# RT32 location (Ventspils, Latvia)
rt32_antenna = RT32()
rt32_antenna.set_location(latitude=57.5535171694, longitude=21.8545525000, elevation=20)


def calculate_antenna_positions(year,month,day,hour,minute):
      
    currentWather = getCurrentWeather()

    # Define constants
    path = ''

    temperature = u.Quantity(currentWather.temperature_2m, unit=u.deg_C)
    pressure = u.Quantity(currentWather.surface_pressure, unit=u.hPa)
    relative_humidity = u.Quantity(currentWather.relative_humidity_2m, unit=u.percent)
    obswl =u.Quantity(50000, unit=u.nm) 

    weather = Weather(temperature, pressure, relative_humidity, obswl)

    observation = SpiralSunObservation(weather,rt32_antenna , year , month , day , hour , minute)

    az_anten, el_anten , az_sun , el_sun , xx1 , yy1, utc = observation.calculatePositions()
    observation.generateFile(path, az_anten , el_anten , utc)  

    return True

def udp_receiver(ip, port, buffer_queue):
    global running
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((ip, port))

    print(f"Esperando mensajes en {ip}:{port}")
    while running:
        data, addr = sock.recvfrom(204)
        if data:
            # Añadir los datos al buffer
            buffer_queue.put((data, addr))
            print(f"Paquete recibido de {addr} y añadido al buffer. Tamaño del buffer: {buffer_queue.qsize()}")
    sock.close()

def data_processor(buffer_queue, output_file):
    global running
    with open(output_file, 'a') as f:
        while running or not buffer_queue.empty():
            if not buffer_queue.empty():
                data, addr = buffer_queue.get()
                print(f"Procesando paquete de {addr}. Tamaño del buffer después de procesar: {buffer_queue.qsize()}")
                # Procesar el paquete recibido
                header = data[:8].decode('utf-8')
                if header.startswith("LNSP"):
                    for i in range(16):
                        timestamp, sensor_data = struct.unpack_from('>L H', data, 8 + i * 12)
                        sensor_data = sensor_data - 0x7FFF  # Ajustar el valor del sensor
                        timestamp_str = datetime.utcfromtimestamp(timestamp / 1000).strftime('%H:%M:%S.%f')[:-3]
                        entry = {"timestamp": timestamp_str, "sensor_data": sensor_data}
                        socketio.emit('new_data', entry)  # Emitir los datos en tiempo real a los clientes
                        f.write(f"{timestamp_str}, {sensor_data}\n")
                    f.flush()  # Asegurar que los datos se escriban en el archivo inmediatamente

def buffer_monitor(buffer_queue):
    global running
    while running:
        print(f"Estado del buffer: {buffer_queue.qsize()} elementos")
        time.sleep(5)  # Imprime el estado del buffer cada 5 segundos

@app.route('/startReceiver', methods=['POST'])
def start_receiver():
    global receiver_thread, processor_thread, monitor_thread, running, output_file

    if not running:
        data = request.get_json()
        date_str = data.get('date')
        time_str = data.get('time')
        if not date_str or not time_str:
            return jsonify({"message": "Fecha u hora no proporcionada"}), 400

        output_file = f"received_data_{date_str}_{time_str.replace(':', '-')}.txt"
        
        running = True
        receiver_thread = threading.Thread(target=udp_receiver, args=("127.0.0.1", 8888, buffer_queue))
        receiver_thread.daemon = True
        receiver_thread.start()

        processor_thread = threading.Thread(target=data_processor, args=(buffer_queue, output_file))
        processor_thread.daemon = True
        processor_thread.start()

        monitor_thread = threading.Thread(target=buffer_monitor, args=(buffer_queue,))
        monitor_thread.daemon = True
        monitor_thread.start()

        return jsonify({"message": f"Receiver started, saving to {output_file}"}), 200
    else:
        return jsonify({"message": "Receiver is already running"}), 200

@app.route('/stopReceiver', methods=['POST'])
def stop_receiver():
    global running
    if running:
        running = False
        # Esperar a que los hilos terminen
        if receiver_thread:
            receiver_thread.join()
        if processor_thread:
            processor_thread.join()
        if monitor_thread:
            monitor_thread.join()

        return jsonify({"message": "Receiver stopped"}), 200
    else:
        return jsonify({"message": "Receiver is not running"}), 200


@app.route('/getCurrentWeather', methods=['POST'])
def getWeatherForecast():
    currentWeather = getCurrentWeather()
  
    return currentWeather,200

@app.route('/calculate_antenna_positions', methods=['POST'])
def calculate_antenna_positions():
    data = request.get_json()
    year = data.get('date')
    month = data.get('time')
    day = data.get('date')
    hour = data.get('time')
    minute = data.get('date')  

    print(year, month, day, hour, minute)
    
    if not all([year, month, day, hour, minute]):
        return jsonify(400, 'Missing date or time information')
    if not all(map(str.isdigit, [year, month, day, hour, minute])):
        return jsonify(400, 'Year, month, day, hour, and minute must be integers')
    if not (len(year) == 4 and len(month) == 2 and len(day) == 2 and len(hour) == 2 and len(minute) == 2):
        return jsonify(400, 'Year must have 4 digits, and month, day, hour, and minute must have 2 digits each')

    antenna_file_thread = threading.Thread(target=calculate_antenna_positions, args=(year,month,day,hour,minute))
    antenna_file_thread.daemon = True
    antenna_file_thread.start()

    return jsonify({'message': 'Antenna positions calculated successfully'})

@app.route('/')
def index():
    return render_template('index.html')

if __name__ == "__main__":
    socketio.run(app, port=5000)
