#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Ьфн 2 18:36:53 2021

@author: kkn
"""

from os.path import dirname, join as pjoin
from scipy.io import wavfile
import pdb
import scipy.io
import math
import random
import matplotlib.pyplot as plt
import numpy as np

# Частота дискретизации файла
signal_samplerate = 40000

# Частота сигнала
signal_frequency = 4000

# Длительность сигнала секунд
signal_length_second = 0.005

# Требуемое соотношение сигнал/шум
traget_signal_noise_ration_db = 60

# Число отсчётов дискретизации на период сигнала
signal_period_sample_count = signal_samplerate / signal_frequency

# Число отсчётов на период опорного сигнала, больше сигнала
reference_more_samples = 21
reference_signal_period_sample_count = reference_more_samples * signal_period_sample_count

# Число отсчётов дискретизации для длительности сигнала
signal_length_sample_count = signal_samplerate * signal_length_second

# Фаза сигнала и фаза опорного сигнала
signal_phase = 180

# Фаза сигналов в отсчётах дискретизации и разность фаз в отсчётах дискретизации
signal_phase_sample    = signal_phase * signal_period_sample_count / 360

# Вывод параметров
print("Частота дискретизации =", signal_samplerate)
print("Частота сигнала =", signal_frequency)
print("Число отсчётов дискретизации на период сигнала =", signal_period_sample_count)
print("Длительность сигнала секунд = ", signal_length_second)
print("Число отсчётов дискретизации на длину сигнала =", signal_length_sample_count)
print("Фаза сигнала, градусы = ", signal_phase, " отсчёты = ", signal_phase_sample)

# Формируем отсчёты синусов сигнала на всю длину
signal_sin_arange           = np.arange(signal_length_sample_count)
signal_sin_value            = np.sin(2 * np.pi * signal_sin_arange / signal_period_sample_count + signal_phase * 2 * np.pi / 360)

# Мощность сигнала
signal_sin_watts = signal_sin_value ** 2

# Средняя мощность сигнала
signal_sin_avg_watts = np.mean(signal_sin_watts) * 32767
signal_sin_avg_watts_db = 10 * np.log10(signal_sin_avg_watts)

# Вывод параметров сигналов
print("Средняя мощность сигнала, дБ =", signal_sin_avg_watts_db)

# Значение шума
noise_avg_db = signal_sin_avg_watts_db - traget_signal_noise_ration_db
noise_avg_watts = 10 ** (noise_avg_db / 10)

# Вывод параметров шума
print("Средняя мощность шума, дБ =", noise_avg_db)

# Генерация белого шума
mean_noise = 0
noise_volts = np.random.normal(mean_noise, np.sqrt(noise_avg_watts), len(signal_sin_watts))

# Добавим шум к сигналу
signal_noise_data = signal_sin_value + noise_volts

# Масштабируем отсчёты, если превысят 1, что бы не было переполнения в int16
scale_value = 1
for i in range(int(signal_length_sample_count)):
    if abs(signal_noise_data[i]) > scale_value:
       scale_value = signal_noise_data[i]
signal_noise_data = signal_noise_data/scale_value

# Формируем отсчёты опорного сигнала на один период
reference_signal_sin_arange = np.arange(reference_signal_period_sample_count)
reference_signal_sin_value_period  = np.sin(2 * np.pi * reference_signal_sin_arange / reference_signal_period_sample_count)
reference_signal_sin_value = np.linspace(0, 0, int(signal_length_sample_count))

# Фильтр для результата перемножения
filter_result = np.linspace(0, 0, int(signal_length_sample_count))
filter_delay_values = np.linspace(0, 0, int(9))
filter_counter = 0
filter_fir_coefficients = np.linspace(0, 0, int(9))

'''
filter_fir_coefficients[0] = 0.0125
filter_fir_coefficients[1] = 0.0866
filter_fir_coefficients[2] = 0.0989
filter_fir_coefficients[3] = 0.139
filter_fir_coefficients[4] = 0.1638
filter_fir_coefficients[5] = 0.1732
filter_fir_coefficients[6] = 0.1638
filter_fir_coefficients[7] = 0.139
filter_fir_coefficients[8] = 0.0989
filter_fir_coefficients[9] = 0.0866
filter_fir_coefficients[10] = 0.0125
'''

filter_fir_coefficients[0] = 0.039294630167254464
filter_fir_coefficients[1] = 0.08585201167310286
filter_fir_coefficients[2] = 0.14444529924004643
filter_fir_coefficients[3] = 0.19238672735213866
filter_fir_coefficients[4] = 0.21093189447050192
filter_fir_coefficients[5] = 0.19238672735213866
filter_fir_coefficients[6] = 0.14444529924004643
filter_fir_coefficients[7] = 0.08585201167310286
filter_fir_coefficients[8] = 0.039294630167254464

# Перемножение сигнала и опорного сигнала
reference_signal_phase_counter = 0
reference_signal_phase_period_counter = 0
multiplication_result = np.linspace(0, 0, int(signal_length_sample_count))
for i in range(int(signal_length_sample_count)):
    multiplication_result[i] = signal_noise_data[i] * reference_signal_sin_value_period[reference_signal_phase_counter]
    reference_signal_sin_value[i] = reference_signal_sin_value_period[reference_signal_phase_counter]

    # Значения линии задержки фильтра
    for j in reversed(range(1, len(filter_delay_values))):
        filter_delay_values[j] = filter_delay_values[j - 1]
    filter_delay_values[0] = multiplication_result[i]

    # Выход фильтра
    for j in range(len(filter_fir_coefficients)):
        filter_result[i] += filter_delay_values[j] * filter_fir_coefficients[j]

    # Счётчик фазы опорного сигнала
    reference_signal_phase_counter += reference_more_samples
    if reference_signal_phase_counter >= reference_signal_period_sample_count:
        reference_signal_phase_counter = int(reference_signal_phase_counter - reference_signal_period_sample_count)
    
    # Каждый полный период входного сигнала коррекция счётчика фазы
    reference_signal_phase_period_counter += 1
    if reference_signal_phase_period_counter >= signal_period_sample_count:
        reference_signal_phase_period_counter = 0
        if filter_result[i] >= 0:
            reference_signal_phase_counter += int(1 + filter_result[i] * 30)
            if reference_signal_phase_counter >= reference_signal_period_sample_count:
                reference_signal_phase_counter = 0
        if filter_result[i] <= -0:
            reference_signal_phase_counter -= int(1 + filter_result[i] * 30)
            if reference_signal_phase_counter < 0:
                reference_signal_phase_counter = int(reference_signal_period_sample_count - 1)
     
t = np.linspace(0, int(signal_length_sample_count), int(signal_length_sample_count))
plt.subplot(3,1,2)
plt.plot(t, signal_sin_value, 'k', t, signal_noise_data, 'b', t, reference_signal_sin_value, 'g', t, multiplication_result, 'r', t, filter_result, 'm')
plt.title('Signal/Reference signal')
plt.ylabel('Value')
plt.xlabel('Sample')
plt.show()
