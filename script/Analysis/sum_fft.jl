include("y_oscillation.jl")
using FFTW


function get_sum_for_each_frac(base_dir, file_list::Vector)
  # make sure the vector length is the same for all files
  vector_length = length(load(joinpath(base_dir, file_list[1]))["x_raw"])
  # initialize the vectors and dictionary
  x = zeros(vector_length)
  y = zeros(vector_length)
  dict = Dict()
  frac_list = []

  # loop through all files to get the list of frac and initialize the dictionary with key
  for file in file_list
    frac = split(file, "_")[2]
    if frac ∉ frac_list
      push!(frac_list, frac)
      # initialize the dictionary
      dict[frac] = [x, y]
    end
  end

  for file in file_list
    frac = split(file, "_")[2]
    data = load(joinpath(base_dir, file))
    dict[frac][1] = dict[frac][1] .+ data["x_raw"]/length(file_list)
    dict[frac][2] = dict[frac][2] .+ data["y_raw"]/length(file_list)
  end
  return dict 
end


# write the graph by using data with CairoMakie
function write_sum_graph(base_dir, file_list, smoothing_num)
  data = get_sum_for_each_frac(base_dir, file_list)
  for keys in keys(data)
    x_raw = data[keys][1]
    y_raw = data[keys][2]
    x = smoothing(x_raw, smoothing_num)
    y = smoothing(y_raw, smoothing_num)
    step_list = [i for i in 1:length(x)]

    fig = Figure()
    ax = Axis(fig[1, 1],
              xlabel = "Time step",
              ylabel = "Sum of x",
              title = "Sum of $(keys)")
    lines!(ax, step_list, x)
    # check wheather the directory exists
    if !isdir(joinpath("Data/oscillation/pic"))
      mkdir(joinpath("Data/oscillation/pic"))
    end
    save(joinpath("Data/oscillation/pic", "x_sum_$(keys)_smoothing$(lpad(string(smoothing_num), 3, "0")).png"), fig)
    println("graph of x Sum of $(keys) is saved.")
    fig = Figure()
    ax = Axis(fig[1, 1],
              xlabel = "Time step",
              ylabel = "Sum of y",
              title = "Sum of $(keys)")
    lines!(ax, step_list, y)
    # check wheather the directory exists
    if !isdir(joinpath("Data/oscillation/pic"))
      mkdir(joinpath("Data/oscillation/pic"))
    end
    save(joinpath("Data/oscillation/pic", "y_sum_$(keys)_smoothing$(lpad(string(smoothing_num), 3, "0")).png"), fig)
    println("graph of y Sum of $(keys) is saved.")
  end
end


function get_wave_length_from_fft(str, frequency::Float64, L::Int64)
  # pick up width from load_dir
  if occursin("frac", str)
    idx = findfirst("frac", str)[1]+4
    width = L * parse(Int, str[idx:idx+1])/10
    println("width = ", width)
    T = 1/frequency
    println("T = ", T)
    velocity = width/T
    return velocity
  else
    println("frac is not included")
  end
end




function write_sum_fft_graph(base_dir, file_list, smoothing_num)
  data = get_sum_for_each_frac(base_dir, file_list)

  # prepare the graph
  for keys in keys(data)
    fig = Figure()
    ax = Axis(fig[1, 1],
          xlabel = "Frequency",
          ylabel = "Amplitude",
          title = "Amplitude of $(keys)",
          # limit the range of x-axis 0:0.5
          limits = (0, 10, 0, 35)
         )
    # get the data and smoothing
    y_original = data[keys][2]
    y = smoothing(y_original, smoothing_num)

    # get the frecuency and amplitude
    y_length = length(y)

    fft_freq = rfftfreq(y_length, 1000)
    # get the amplitude
    y_fft = abs.(rfft(y))
    println(keys)
    println("peak frequency: ", fft_freq[argmax(y_fft)])
    velocity = get_wave_length_from_fft(keys, fft_freq[argmax(y_fft)], 30)
    println("velocity = ", velocity)

    # plot the graph with the data
    lines!(ax, fft_freq, y_fft, label = keys*": period = "*string(1/fft_freq[argmax(y_fft)]))


    # check wheather the directory exists
    if !isdir(joinpath("Data/oscillation/pic"))
      mkdir(joinpath("Data/oscillation/pic"))
    end
    # display the label
    axislegend(ax)
    # save the graph
    save_dir = joinpath("Data/oscillation/pic", "sum_fft_smoothing_$(keys)_$(lpad(string(smoothing_num), 3, "0")).png")

    save(save_dir, fig)
    println("graph of Sum of FFT is saved.")
    # initialize the graph
  end
end


base_dir = "./Data/oscillation/smoothed"
file_list = readdir(base_dir)
smoothing_num = 1


write_sum_graph(base_dir, file_list, smoothing_num)
write_sum_fft_graph(base_dir, file_list, smoothing_num)
