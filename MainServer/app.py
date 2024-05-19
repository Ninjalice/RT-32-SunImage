from flask import Flask, request, jsonify, render_template, abort
import threading
import queue
import socket
import struct
import time
from datetime import datetime
import AntenaPos

app = Flask(__name__)

# Tamaño del buffer de la cola
BUFFER_SIZE = 2048
OUTPUT_FILE = "received_data.txt"
buffer_queue = queue.Queue(BUFFER_SIZE)
receiver_thread = None
processor_thread = None
running = False

def udp_receiver(ip, port, buffer_queue):
    global running
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((ip, port))

    print(f"Esperando mensajes en {ip}:{port}")
    while running:
        data, addr = sock.recvfrom(204)
        if data:
            buffer_queue.put((data, addr))
            print(f"Paquete recibido de {addr} y añadido al buffer. Tamaño del buffer: {buffer_queue.qsize()}")

def data_processor(buffer_queue, output_file):
    with open(output_file, 'a') as f:
        while running:
            if not buffer_queue.empty():
                data, addr = buffer_queue.get()
                print(f"Procesando paquete de {addr}. Tamaño del buffer después de procesar: {buffer_queue.qsize()}")
                header = data[:8].decode('utf-8')
                if header.startswith("LNSP"):
                    for i in range(16):
                        timestamp, sensor_data = struct.unpack_from('>L H', data, 8 + i * 12)
                        sensor_data = sensor_data - 0x7FFF
                        timestamp_str = datetime.utcfromtimestamp(timestamp / 1000).strftime('%H:%M:%S.%f')[:-3]
                        f.write(f"{timestamp_str}, {sensor_data}\n")
                    f.flush()

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/start_receiver', methods=['POST'])
def start_receiver():
    global running, receiver_thread, processor_thread
    if not running:
        ip = request.json.get('ip')
        port = request.json.get('port')
        if not (ip and port):
            abort(400, 'Missing IP or Port in request')
        try:
            port = int(port)
        except ValueError:
            abort(400, 'Port must be an integer')
        
        running = True
        receiver_thread = threading.Thread(target=udp_receiver, args=(ip, port, buffer_queue))
        receiver_thread.daemon = True
        receiver_thread.start()

        processor_thread = threading.Thread(target=data_processor, args=(buffer_queue, OUTPUT_FILE))
        processor_thread.daemon = True
        processor_thread.start()

        return jsonify({'status': 'started'})
    return jsonify({'status': 'already running'})

@app.route('/stop_receiver', methods=['POST'])
def stop_receiver():
    global running
    running = False
    return jsonify({'status': 'stopped'})

@app.route('/get_buffer_size', methods=['GET'])
def get_buffer_size():
    return jsonify({'buffer_size': buffer_queue.qsize()})


@app.route('/calculate_antenna_positions', methods=['POST'])
def calculate_antenna_positions():
    year = request.form.get('year')
    month = request.form.get('month')
    day = request.form.get('day')
    hour = request.form.get('hour')
    minute = request.form.get('minute')

    if not all([year, month, day, hour, minute]):
        abort(400, 'Missing date or time information')
    if not all(map(str.isdigit, [year, month, day, hour, minute])):
        abort(400, 'Year, month, day, hour, and minute must be integers')
    if not (len(year) == 4 and len(month) == 2 and len(day) == 2 and len(hour) == 2 and len(minute) == 2):
        abort(400, 'Year must have 4 digits, and month, day, hour, and minute must have 2 digits each')

    current = AntenaPos.weatherDataBase()
    AntenaPos.create_antenna_file(year, month, day, hour, minute, current.Variables(7).Value(), 
                                  current.Variables(0).Value(), current.Variables(1).Value(), 50000)

    return jsonify({'message': 'Antenna positions calculated successfully'})

@app.route('/check_weather', methods=['POST'])
def check_weather():
    data = request.json
    option = data.get('option')
    
    # Llamar a la función de tu librería para obtener el pronóstico del tiempo
    if option == "Current":
        weather_output = AntenaPos.weatherForecast(1)
    elif option == "Next Day":
        weather_output = AntenaPos.weatherForecast(2)
    elif option == "Next 3 Days":
        weather_output = AntenaPos.weatherForecast(3)
    else:
        weather_output = "Invalid option selected"

    return jsonify({'weather_output': weather_output})

if __name__ == '__main__':
    app.run(debug=True)
