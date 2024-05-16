import socket
import struct
import threading
import queue
import time
from datetime import datetime
from flask import Flask, render_template, request, jsonify
from flask_socketio import SocketIO, emit
import json

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

@app.route('/')
def index():
    return render_template('index.html')

if __name__ == "__main__":
    socketio.run(app, port=5000)
