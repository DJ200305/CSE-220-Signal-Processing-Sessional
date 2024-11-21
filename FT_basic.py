import numpy as np
import matplotlib.pyplot as plt

# Define the interval and function and generate appropriate x values and y values
def parabola(x):
    return np.where((x >= -2) & (x <= 2), x**2, 0)
def sawtooth(x):
    freq = 1/4
    amp = 1.0
    c = (x+2)*freq
    return np.where((x_values >= -2) & (x_values <= 2),amp*(c-np.floor(c)),0)
def triangle(x):
    freq = 1/4
    amp = 1.0
    c = (x+2)*freq
    return np.where((x_values >= -2) & (x_values <= 2),amp*(2*np.abs(2*(c-np.floor(c)-0.5))),0)
def rectangular(x_values):
    return np.where((x_values >= -2) & (x_values <= 2), 1, 0)
x_values = np.linspace(-10, 10, 1000)
y_values = sawtooth(x_values)

# Plot the original function
plt.figure(figsize=(12, 4))
plt.plot(x_values, y_values, label="Original y = x^2")
plt.title("Original Function (y = x^2)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()


# Define the sampled times and frequencies
sampled_times = x_values
frequencies = np.linspace(-5, 5, 1000)

 # Fourier Transform 
def fourier_transform(signal, frequencies, sampled_times):
    num_freqs = len(frequencies)
    ft_result_real = np.zeros(num_freqs)
    ft_result_imag = np.zeros(num_freqs)
    
    # Store the fourier transform results for each frequency. Handle the real and imaginary parts separately
    # use trapezoidal integration to calculate the real and imaginary parts of the FT
    for i in range(num_freqs):
        ft_result_real[i] = np.trapz(signal * np.cos(2*np.pi*frequencies[i]*sampled_times),sampled_times,sampled_times[1]-sampled_times[0])
        ft_result_imag[i] = -np.trapz(signal * np.sin(2*np.pi*frequencies[i]*sampled_times),sampled_times,sampled_times[1]-sampled_times[0])
    return ft_result_real, ft_result_imag

# Apply FT to the sampled data
ft_data = fourier_transform(y_values, frequencies, sampled_times)
#  plot the FT data
plt.figure(figsize=(12, 6))
plt.plot(frequencies, np.sqrt(ft_data[0]**2 + ft_data[1]**2))
plt.title("Frequency Spectrum of y = x^2")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.show()


# Inverse Fourier Transform 
def inverse_fourier_transform(ft_signal, frequencies, sampled_times):
    n = len(sampled_times)
    reconstructed_signal = np.zeros(n)
    # Reconstruct the signal by summing over all frequencies for each time in sampled_times.
    # use trapezoidal integration to calculate the real part
    # You have to return only the real part of the reconstructed signal
    for i in range(n):
        a = np.trapz(ft_signal[0]*np.cos(2*np.pi*frequencies*sampled_times[i]),frequencies,frequencies[1]-frequencies[0])
        b = -np.trapz(ft_signal[1]*np.sin(2*np.pi*frequencies*sampled_times[i]),frequencies,frequencies[1]-frequencies[0])
        reconstructed_signal[i] = a + b   
    return reconstructed_signal

# Reconstruct the signal from the FT data
reconstructed_y_values = inverse_fourier_transform(ft_data, frequencies, sampled_times)
# Plot the original and reconstructed functions for comparison
plt.figure(figsize=(12, 4))
plt.plot(x_values, y_values, label="Original y = x^2", color="blue")
plt.plot(sampled_times, reconstructed_y_values, label="Reconstructed y = x^2", color="red", linestyle="--")
plt.title("Original vs Reconstructed Function (y = x^2)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()
