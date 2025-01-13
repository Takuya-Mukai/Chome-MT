include("../src/Chemo-MT.jl")

using FFTW
using CairoMakie
sineWave(freq, t) = sin(2π*freq*t)

dt = 0.1
t = 0.0:dt:10.0

freq_of_sine_wave = 2
signal = sineWave.(freq_of_sine_wave, t)

N = length(t)
fft_freq = rfftfreq(N, 1/dt)
F = abs.(rfft(signal))

fig = Figure()
ax = Axis(fig[1, 1])
lines!(fft_freq, F)

fig
