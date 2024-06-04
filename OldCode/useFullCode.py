# Encontrar el índice del elemento en la columna 'columna'
# indice = data_df[data_df['UTC RCP 1'] == 97705372].index

# print(data_df['UTC RCP 1'])
# print(indice)



# # Hora inicial en milisegundos (08:40:00)
# hora_inicial_ms = (8 * 60 * 60 * 1000) + (41 * 60 * 1000) + 0 * 1000

# # Intervalo de tiempo en milisegundos (por ejemplo, 2 segundos)
# intervalo_ms = 2 * 1000

# # Función para calcular los milisegundos desde la hora inicial
# def calcular_milisegundos(fila, hora_inicial_ms, intervalo_ms):
#     return hora_inicial_ms + (fila * intervalo_ms)

# # Crear una nueva columna 'base_miliseconds' basada en el índice de la fila
# data_df['base_miliseconds'] = data_df.index.to_series().apply(calcular_milisegundos, args=(hora_inicial_ms, intervalo_ms))

# # Función para convertir milisegundos a horas
# def milisegundos_a_horas(ms):
#     return ms / (60 * 60 * 1000)

# # Para cada columna que empieza con 'UTC', restar su valor a 'base_miliseconds' y convertir el resultado a horas
# for col in data_df.columns:
#     if col.startswith('UTC'):
#         data_df[col] = ( data_df['base_miliseconds'] - (data_df['UTC LCP 11'] - data_df[col]) ).apply(milisegundos_a_horas)

# # Eliminar la columna 'base_miliseconds'
# data_df = data_df.drop(columns=['base_miliseconds'])

# print(data_df.describe())