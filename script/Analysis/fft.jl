using FFTW
include("../../src/Chemo-MT.jl")

function get_subdirectories(directory::String)
    # Get all entries in the specified directory
    entries = readdir(directory)
    # Filter only directories and return them as a vector
    subdirs = [entry*"/" for entry in entries if occursin("ystep", entry) && isdir(joinpath(directory, entry))]
    if length(subdirs) == 0
      error("No subdirectories found in $(directory). check the function of get_subdirectories")
    end
    return subdirs
end


function get_sorted_file(base_dir::String; load_dir::String="")
    loaddir = ""
    if load_dir ==""
        loaddir = base_dir
    else
        loaddir = base_dir * load_dir
    end
    file_list = get_file_list(loaddir)
    sorted_file = Division_into_smaple(file_list)
    return sorted_file
end


function get_file_list(dir)
    all_file = readdir(dir)
    file_list = []

    for i in eachindex(all_file)
        if split(all_file[i], ".")[end] == "jld2" && split(all_file[i], ".")[1] != ""
            append!(file_list, [all_file[i]])
        else
            nothing
        end
    end
    return file_list
end


function Division_into_smaple(file_list)

  sample_list = []
  for name in file_list
    sample_term = split(name, "_")[1]
    if length(sample_term) == 8
      Nsample = parse(Int, sample_term[end-1:end])
    else
      Nsample = parse(Int, sample_term[end])
    end

    if Nsample ∉ sample_list
      append!(sample_list, Nsample)
    end
  end
  Nsample = length(sample_list)
  sorted_file = [[] for _ in 1:Nsample]
  for name in file_list
    sample_term = split(name, "_")[1]
    if length(sample_term) == 8
      Nsample = parse(Int, sample_term[end-1:end])
    else
      Nsample = parse(Int, sample_term[end])
    end
    append!(sorted_file[Nsample], [name])
  end
  return sorted_file
end


function xy_in_time_all(file::String)
  data = load(file)
  all_ori = data["all_ori"]

  Nstep = length(all_ori)
  all_yop = zeros(Nstep)
  all_xop = zeros(Nstep)
  for i in 1:Nstep
    all_xop[i], all_yop[i] = sum(angle2dir.(all_ori[i]))/length(all_ori[i])
  end
  return all_xop, all_yop
end


function xy_in_time_region(file::String, spatio_func::Function)
  data = load(file)
  all_ori = data["all_ori"]
  all_pos = data["all_pos"]

  Nstep = length(all_ori)
  all_yop = zeros(Nstep)
  all_xop = zeros(Nstep)
  for i in 1:Nstep
    x = getindex.(all_pos[i], 1)
    y = getindex.(all_pos[i], 2)
    num_in_region = sum(spatio_func.(x, y, data["paras"].L))
    all_xop[i], all_yop[i] = sum(angle2dir.(all_ori[i]))/num_in_region
  end
  return all_xop, all_yop
end

function average_loc(file::String)
  data = load(file)
  all_pos = data["all_pos"]
  paras = data["paras"]
  Nstep = length(all_pos)
  average_x_list = zeros(Nstep)
  average_y_list = zeros(Nstep)
  for i in 1:Nstep
    x = getindex.(all_pos[i], 1)
    y = getindex.(all_pos[i], 2)
    average_x = sum(x)/length(x)
    average_y = sum(y)/length(y)
    average_x_list[i] = average_x-paras.L/2
    average_y_list[i] = average_y-paras.L/2
  end
  return average_x_list, average_y_list
end



function flipy(x::Vector, l::Int64)
  if sum(x[1:l]) < 0
    return -x
  else
    return x
  end
end


function flipx(x::Vector, l::Int64)
  if sum(x[end-l:end]) < 0
    return -x
  else
    return x
  end
end


function flipy_loc(x::Vector, l::Int64)
  if sum(x[1:l]) < 0
    return -x
  else
    return x
  end
end


function PickUpNumFromDirName(DirName::String)
  m = match(r"Lfrac(\d+)", DirName)

  if m !== nothing
      value = parse(Float64, m.captures[1])/10
  else
    error("could not pick up number from $(DirName)")
  end
  return value
end


function inTheRegionNum(file::String, spatio_func::Function)
  data = load(file)
  all_pos = data["all_pos"]
  num_in_region_list = zeros(length(all_pos))
  for i in 1:length(all_pos)
    num_in_region = sum(spatio_func.(getindex.(all_pos[i], 1), getindex.(all_pos[i], 2), data["paras"].L))
    num_in_region_list[i] = num_in_region
  end
  return num_in_region_list
end


# define base directory
base_dir = "./Data/oscillation_long/"
# get all subdirectories
load_dir_list = get_subdirectories(base_dir)
flip_num = 200

  f3 = Figure(resolution = (800, 600))

  dir_fft = Dict()
  dir_fft_all = Dict()
  dir_fft_region = Dict()
  dir_fft_averageLoc_x = Dict()
  dir_fft_averageLoc_y = Dict()
  dir_NumInRegion = Dict()
  dir_AlignedVx = Dict()
  dir_yOscillation = Dict()
  dir_averageLoc_x = Dict()
  dir_averageLoc_y = Dict()

  for (color, load_dir) in enumerate(load_dir_list)
    sorted_files = get_sorted_file(base_dir; load_dir = load_dir)

    function spatio_func(x, y, L)
      idx = findfirst("frac", load_dir)[1]+4
      L_frac = parse(Float64, load_dir[idx:idx+1])/10

      up = L / 2 + L * L_frac / 2
      low = L / 2 - L * L_frac / 2

      if y > low && y < up
          return 1.0
      else
          return 0.0
      end
    end

    eta_fft = Dict()
    eta_fft_all = Dict()
    eta_fft_region = Dict()
    eta_fft_averageLoc_x = Dict()
    eta_fft_averageLoc_y = Dict()
    eta_NumInRegion = Dict()
    eta_AlignedVx = Dict()
    eta_yOscillation = Dict()
    eta_averageLoc_x = Dict()
    eta_averageLoc_y = Dict()

    x_length = length(xy_in_time_all(joinpath(base_dir, load_dir, sorted_files[1][2]))[1])
    for eta in 2:10

      sample_num = length(sorted_files)
      time_step = 1:x_length
      num_in_region_list = zeros(x_length)
      x_region, y_region = zeros(x_length), zeros(x_length)
      x_all, y_all = zeros(x_length), zeros(x_length)

      x_average, y_average = zeros(x_length), zeros(x_length)
      for sample_file in sorted_files
        target_file = sample_file[eta]

        all_xop_all, all_yop_all = xy_in_time_all(joinpath(base_dir, load_dir, target_file))
        x_flip_all, y_flip_all = flipx(all_xop_all, flip_num), flipy(all_yop_all, flip_num)
        x_all .+= x_flip_all/sample_num
        y_all .+= y_flip_all/sample_num

        all_xop_region, all_yop_region = xy_in_time_region(joinpath(base_dir, load_dir, target_file), spatio_func)
        x_flip_region, y_flip_region = flipx(all_xop_region, flip_num), flipy(all_yop_region, flip_num)
        x_region .+= x_flip_region/sample_num
        y_region .+= y_flip_region/sample_num

        data = load(joinpath(base_dir, load_dir, target_file))
        all_pos = data["all_pos"]
        num_in_region = inTheRegionNum(joinpath(base_dir, load_dir, target_file), spatio_func)
        num_in_region_list .+= num_in_region/sample_num


        # average location
        x_loc, y_loc =  average_loc(joinpath(base_dir, load_dir, target_file))
        x_average .+= x_loc/sample_num
        y_average .+= flipy_loc(y_loc, flip_num)/sample_num
      end


      # FFT by velocity
      fft_freq = rfftfreq(x_length, 1/(60*0.001))
      y_fft_all = abs.(rfft(y_all))
      y_fft_region = abs.(rfft(y_region))
      eta_fft_all["$eta"] = [fft_freq, y_fft_all]
      eta_fft_region["$eta"] = [fft_freq, y_fft_region]

      # FFT by average location
      fft_freq = rfftfreq(x_length, 1/(60*0.001))
      x_fft_averageLoc = abs.(rfft(x_average))
      y_fft_averageLoc = abs.(rfft(y_average))
      eta_fft_averageLoc_y["$eta"] = [fft_freq, y_fft_averageLoc]
      eta_fft_averageLoc_x["$eta"] = [fft_freq, x_fft_averageLoc]

      eta_NumInRegion["$eta"] = [1:length(num_in_region_list), num_in_region_list]
      # eta_AlignedVx["$eta"] = [1:length(x_region), x_region]
      eta_AlignedVx["$eta"] = [1:length(x_all), x_all]
      eta_yOscillation["$eta"] = [1:length(y_all), y_all]
      eta_averageLoc_x["$eta"] = [1:length(x_average), x_average]
      eta_averageLoc_y["$eta"] = [1:length(y_average), y_average]
      println("band width: $(PickUpNumFromDirName(load_dir)), eta: $eta")
    end
    dir_fft[load_dir] = eta_fft
    dir_fft_all[load_dir] = eta_fft_all
    dir_fft_region[load_dir] = eta_fft_region
    dir_fft_averageLoc_y[load_dir] = eta_fft_averageLoc_y
    dir_fft_averageLoc_x[load_dir] = eta_fft_averageLoc_x
    dir_NumInRegion[load_dir] = eta_NumInRegion
    dir_AlignedVx[load_dir] = eta_AlignedVx
    dir_yOscillation[load_dir] = eta_yOscillation
    dir_averageLoc_x[load_dir] = eta_averageLoc_x
    dir_averageLoc_y[load_dir] = eta_averageLoc_y
  end


  ## pick up the top 3 frequency of each eta and dir_fft
  top_freq_dir = Dict()
  top_freq_eta = Dict()
  for (key_dir, value_dir) in dir_fft_region
    top_freq_eta = Dict()
    for (key_eta, value_eta) in value_dir
      que_for_sort = sortperm(value_eta[2], rev=true)
      sorted_freq = value_eta[1][que_for_sort]
      sorted_freq_amp = value_eta[2][que_for_sort]
      top_freq_eta[key_eta] = [sorted_freq[1:3], sorted_freq_amp[1:3]]
    end
    top_freq_dir[key_dir] = top_freq_eta
  end


  color_list = [:red, :blue, :green, :yellow, :purple, :orange, :black]

  f = Figure(resolution = (1000, 800), fontsize = 13)
  ax = Axis(f[1,1], xlabel = "η", ylabel = "frequency")
  maker_size = 20
  handles = [MarkerElement(marker=:circle, color = color, markersize = maker_size)
             for color in color_list[1:length(top_freq_dir)]]
  labels = string.(PickUpNumFromDirName.(collect(keys(top_freq_dir))))
  legend = Legend(f, handles, labels, "each bandwidth")
  f[1,2] = legend

  Label(f[0, :], "average of top 3 frequency")
  for (index, (key_dict, value_dict)) in enumerate(top_freq_dir)
    for (key_eta, value_eta) in value_dict
      scatter!(ax, parse(Int, key_eta)/10, sum(value_eta[1].*value_eta[2])/sum(value_eta[2]),
               markersize=20, color=color_list[index], label="bandwidth:$( PickUpNumFromDirName(key_dict))")
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/Top_frequency&eta.png", f)

  #plot the graph
  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_fft = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "frequency", ylabel = "amplitude",
                   title = "η = $((i+1)*2/10)", limits = ((0, 0.04), (0, 300))) for i in 1:9]
  Label(f[0, :], "FFT Analysis")
  dir_fft_length = length(dir_fft)
  for (index, (key_dict, value_dict)) in enumerate(dir_fft_all)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_fft)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          # write the legend at the end of the dir_fft
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/fft_all.png",f)

  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_fft = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "frequency", ylabel = "amplitude",
                   title = "η = $((i+1)*2/10)", limits = ((0, 0.07), (0, 160))) for i in 1:9]
  Label(f[0, :], "FFT Analysis")
  dir_fft_length = length(dir_fft_region)
  for (index, (key_dict, value_dict)) in enumerate(dir_fft_region)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_fft)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          # write the legend at the end of the dir_fft
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/fft_region.png",f)

  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_fft = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "frequency", ylabel = "amplitude",
                   title = "η = $((i+1)*2/10)", limits = ((0, 0.04), nothing)) for i in 1:9]
  Label(f[0, :], "FFT Analysis of average location in y direction")
  dir_fft_length = length(dir_fft)
  for (index, (key_dict, value_dict)) in enumerate(dir_fft_averageLoc_y)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_fft)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/fft_averageloc_y.png",f)

  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_fft = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "frequency", ylabel = "amplitude",
                   title = "η = $((i+1)*2/10)", limits = ((0, 0.12), nothing)) for i in 1:9]
  Label(f[0, :], "FFT Analysis of average location in x location")
  dir_fft_length = length(dir_fft)
  for (index, (key_dict, value_dict)) in enumerate(dir_fft_averageLoc_x)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_fft)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/fft_averageloc_x.png",f)


  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_NumInRegion = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Timestep", ylabel = "Number of agent in the band",
                           title = "η = $((i+1)*2/10)", limits=(nothing, (0, 900))) for i in 1:9]
  Label(f[0, :], "Number of agents in the region")
  for (index, (key_dict, value_dict)) in enumerate(dir_NumInRegion)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_NumInRegion)
        if key_eta == string(i+1)
          bandwidth = PickUpNumFromDirName(key_dict)
          lines!(ax, value_eta[1], 900*bandwidth*ones(length(value_eta[1])), color = :gray)

          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/AgentNumber.png",f)

  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_AlignedVx = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "v_x",
                         title = "η = $((i+1)*2/10)", limits = (nothing, (-1.3, 1.3))) for i in 1:9]
  Label(f[0, :], "Aligned v_x")
  for (index,(key_dict, value_dict)) in enumerate(dir_AlignedVx)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_AlignedVx)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/v_x_all.png",f)

  for (index, (key_dict, value_dict)) in enumerate(dir_yOscillation)
    f = Figure(resolution = (1500, 1000), fontsize = 13)
    Axes_yOscillation = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "v_y",
                              title = "η = $((i+1)*2/10)", limits = (nothing, (-0.15, 0.15))) for i in 1:9]
    Label(f[0, :], "Oscillation in y direction")
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_yOscillation)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
    save("/home/muta/Code/slidev/project/b8/last/src/yOscillation_$(PickUpNumFromDirName(key_dict)).png",f)
  end

  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_averageLoc_x = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "x",
                            title = "η = $((i+1)*2/10)", limits = (nothing, nothing)) for i in 1:9]
  Label(f[0, :], "Average location in x direction")
  for (index, (key_dict, value_dict)) in enumerate(dir_averageLoc_x)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_averageLoc_x)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/averageLoc_x.png",f)

  f = Figure(resolution = (1500, 1000), fontsize = 13)
  Axes_averageLoc_y = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "y",
                            title = "η = $((i+1)*2/10)", limits = (nothing, nothing)) for i in 1:9]
  Label(f[0, :], "Average location in y direction")
  for (index, (key_dict, value_dict)) in enumerate(dir_averageLoc_y)
    for (key_eta, value_eta) in value_dict
      for (i, ax) in enumerate(Axes_averageLoc_y)
        if key_eta == string(i+1)
          lines!(ax, value_eta[1], value_eta[2], label = "bandwidth:$( PickUpNumFromDirName(key_dict))")
          if index == dir_fft_length
            Legend(f[((parse(Int,key_eta)-2)÷3+1), (parse(Int,key_eta)-2)%3*2], ax)
          end
        end
      end
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/averageLoc_y.png",f)

function write_graph(file::String)
  data = load(file)
  all_ori = data["all_ori"]

  Nstep = length(all_ori)
  all_yop = zeros(Nstep)
  all_xop = zeros(Nstep)
  for i in 1:Nstep
    all_xop[i], all_yop[i] = sum(angle2dir.(all_ori[i]))/length(all_ori[i])
  end
  f = Figure()
  ax = Axis(f[1,1], xlabel = "x", ylabel = "v")
  lines!(ax, 1:length(all_xop), all_xop, label = "x")
  lines!(ax, 1:length(all_yop), all_yop, label = "y")
  #show legend
  Legend(f[1,2], ax)
  f
end
