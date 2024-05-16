import socket
import time
import struct

def sp_channel(channel, pol):
    # Aquí puedes incluir la lógica para configurar el canal y la polarización
    pass

def read_sensor_data():
    # Simula la lectura de un sensor (como ADS1115 en el Arduino)
    # Devuelve un valor ficticio de muestra
    return 12345

def udp_transmitter(ip, port):
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    strh = "LNSP____"
    strh_bytes = strh.encode('utf-8')

    try:
        print("Enviando mensajes...")
        while True:
            packet_buffer = bytearray(204)

            # Leer y enviar datos en modo RCP
            sp_channel(1, 'RCP')
            for i in range(16):
                time_sample_start = int(time.time() * 1000) & 0xFFFFFFFF
                sensor_data = read_sensor_data() & 0xFFFF
                struct.pack_into('>L', packet_buffer, 8 + 12 * i, time_sample_start)
                struct.pack_into('>H', packet_buffer, 8 + 12 * i + 8, sensor_data + 0x7FFF)

            packet_buffer[:8] = strh_bytes
            packet_buffer[4] = ord('R')
            sock.sendto(packet_buffer, (ip, port))

            # Leer y enviar datos en modo LCP
            packet_buffer = bytearray(204)
            sp_channel(1, 'LCP')
            for i in range(16):
                time_sample_start = int(time.time() * 1000) & 0xFFFFFFFF
                sensor_data = read_sensor_data() & 0xFFFF
                struct.pack_into('>L', packet_buffer, 8 + 12 * i, time_sample_start)
                struct.pack_into('>H', packet_buffer, 8 + 12 * i + 8, sensor_data + 0x7FFF)

            packet_buffer[:8] = strh_bytes
            packet_buffer[4] = ord('L')
            sock.sendto(packet_buffer, (ip, port))

            # Esperar 1 milisegundo antes de enviar el siguiente paquete
            time.sleep(1)

    except Exception as e:
        print(f"Error enviando el mensaje: {e}")
    finally:
        sock.close()

if __name__ == "__main__":
    ip = "127.0.0.1"  # Dirección IP del receptor
    port = 8888       # Puerto del receptor
    udp_transmitter(ip, port)