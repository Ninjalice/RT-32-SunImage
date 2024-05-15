import socket
import struct
import threading
import queue
import time
from datetime import datetime

# Tamaño del buffer de la cola
BUFFER_SIZE = 2048
OUTPUT_FILE = "received_data.txt"

def udp_receiver(ip, port, buffer_queue):
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((ip, port))

    print(f"Esperando mensajes en {ip}:{port}")
    while True:
        data, addr = sock.recvfrom(204)
        if data:
            # Añadir los datos al buffer
            buffer_queue.put((data, addr))
            print(f"Paquete recibido de {addr} y añadido al buffer. Tamaño del buffer: {buffer_queue.qsize()}")

def data_processor(buffer_queue, output_file):
    with open(output_file, 'a') as f:
        while True:
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
                        f.write(f"{timestamp_str}, {sensor_data}\n")
                    f.flush()  # Asegurar que los datos se escriban en el archivo inmediatamente

def buffer_monitor(buffer_queue):
    while True:
        print(f"Estado del buffer: {buffer_queue.qsize()} elementos")
        time.sleep(5)  # Imprime el estado del buffer cada 5 segundos

if __name__ == "__main__":
    ip = "127.0.0.1"  # Escuchar en la IP local
    port = 8888       # Puerto para recibir mensajes

    # Crear una cola para el buffer
    buffer_queue = queue.Queue(BUFFER_SIZE)

    # Iniciar el hilo del receptor UDP
    receiver_thread = threading.Thread(target=udp_receiver, args=(ip, port, buffer_queue))
    receiver_thread.daemon = True
    receiver_thread.start()

    # Iniciar el procesador de datos
    processor_thread = threading.Thread(target=data_processor, args=(buffer_queue, OUTPUT_FILE))
    processor_thread.daemon = True
    processor_thread.start()

    # Iniciar el monitor del buffer
    monitor_thread = threading.Thread(target=buffer_monitor, args=(buffer_queue,))
    monitor_thread.daemon = True
    monitor_thread.start()

    # Mantener el hilo principal en ejecución
    while True:
        time.sleep(1)
