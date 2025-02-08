using FFTW
using GLMakie
GLMakie.activate!()
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

  #* find  number of sample
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
  for (i, name) in enumerate(file_list)
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
  # , xmin, xmax, ymin, ymax)
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
    # p_inside = all_ori[i][is_inside_region.(all_pos[i], xmin, xmax, ymin, ymax)]
    # all_xop[i], all_yop[i] = mean(angle2dir.(p_inside))
    all_xop[i], all_yop[i] = sum(angle2dir.(all_ori[i]))/length(all_ori[i])
  end
  return all_xop, all_yop
end


function xy_in_time_region(file::String)
  # , xmin, xmax, ymin, ymax)
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
    # p_inside = all_ori[i][is_inside_region.(all_pos[i], xmin, xmax, ymin, ymax)]
    # all_xop[i], all_yop[i] = mean(angle2dir.(p_inside))
    all_xop[i], all_yop[i] = sum(angle2dir.(all_ori[i]))/num_in_region
  end
  return all_xop, all_yop
end


function aline(x)
  if x < 0
    return -x
  else
    return x
  end
end

function xy_in_time_with_aline_x(file::String, spatio_func::Function)
  # , xmin, xmax, ymin, ymax)
  data = load(file)
  all_ori = data["all_ori"]
  all_pos = data["all_pos"]

  Nstep = length(all_ori)
  all_yop = zeros(Nstep)
  all_xop = zeros(Nstep)
  num_in_region_list = zeros(Nstep)
  for i in 1:Nstep
    x = getindex.(all_pos[i], 1)
    y = getindex.(all_pos[i], 2)
    num_in_region = sum(spatio_func.(x, y, data["paras"].L))
    # all_xop[i] = sum( spatio_func.(x, y, data["paras"].L) .* cos.(all_ori[i]) )/num_in_region
    # all_yop[i] = sum( spatio_func.(x, y, data["paras"].L) .* sin.(all_ori[i]) )/num_in_region
    all_xop[i] = sum( cos.(all_ori[i]) )/num_in_region
    all_yop[i] = sum( sin.(all_ori[i]) )/num_in_region
    num_in_region_list[i] = num_in_region
  end
  return all_xop, all_yop, num_in_region_list
end


function flipy(x::Vector, l::Int64)
  # if sum(x[1:l]) < 0, flip x
  # if sum(x[1:l]) > 0, do nothing
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

function scale_and_count(data, bins)
    hist = zeros(Int, bins)  # initialize
    for x in data
        if x != 0  # 0 を無視
            idx = clamp(floor(Int, 1 + bins * x), 1, bins)  # scaling
            hist[idx] += 1  # 対応するインデックスに+1
        end
    end
    return hist
end
data = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 0.0, 0.0]
histogram = scale_and_count(data, 10)
println(histogram)

function distribution_of_x_in_the_region(file::String, spatio_func::Function, division_num::Int64)
  # , xmin, xmax, ymin, ymax)
  data = load(file)
  all_ori = data["all_ori"]
  all_pos = data["all_pos"]

  Nstep = length(all_ori)
  agent_num = length(all_ori[1])
  num_in_region = length(all_ori[1])
  num_in_region_list = zeros(Nstep)
  vy_list = zeros(Nstep, agent_num)
  vx_list = zeros(Nstep, agent_num)

  # make list of vx and vy in the region
  for i in 1:Nstep
    x = getindex.(all_pos[i], 1)
    y = getindex.(all_pos[i], 2)
    num_in_region = sum(spatio_func.(x, y, data["paras"].L))
    # vx = [if spatio_func(x[j], y[j], data["paras"].L) == 1.0
    #     aline(cos(all_ori[i][j]))
    #   else
    #     0.0
    #   end for j in 1:length(x)]
    # vy = [if spatio_func(x[j], y[j], data["paras"].L) == 1.0
    #     aline(sin(all_ori[i][j]))
    #   else
    #     0.0
    #   end for j in 1:length(x)]
    vx = [cos(all_ori[i][j]) for j in 1:length(x)]
    vy = [sin(all_ori[i][j]) for j in 1:length(x)]
    num_in_region_list[i] = num_in_region
    vx_list[i, :] = vx
    vy_list[i, :] = vy
  end

  # make distribution of vx in the region
  vx_distribution = zeros(Nstep, division_num)
  vy_distribution = zeros(Nstep, division_num)
  for i in 1:Nstep
    vx_distribution[i,:] = scale_and_count(vx_list[i, :], division_num)
    vy_distribution[i,:] = scale_and_count(vy_list[i, :], division_num)
    # for j in 1:agent_num
    #   if vx_list[i, j] != 0.0
    #     idx = Int(floor(vx_list[i, j] * division_num + 1))
    #     vx_distribution[i, idx] += 1
    #   end
    #   if vy_list[i, j] != 0.0
    #     idx = Int(floor(vy_list[i, j] * division_num + 1))
    #     vy_distribution[i, idx] += 1
    #   end
  end
  return vx_distribution, vy_distribution
end

function write_3d_graph(file::String, spatio_func::Function, division_num::Int64)
  vx_distribution, vy_distribution = distribution_of_x_in_the_region(file, spatio_func, division_num)
  Nstep = length(vx_distribution[:, 1])
  division = collect(0:1/(division_num-1):1)
  fig = Figure()
  ax = Axis3(fig[1, 1], xlabel = "time step", ylabel = "alined v_x", zlabel = "number of agent")
  surface!(ax, 1:Nstep, division, vx_distribution, color = :viridis)
  display(fig)
end

function PickUpNumFromDirName(DirName::String)
  m = match(r"Lfrac(\d+)", DirName)

  # 数字を Float64 に変換
  if m !== nothing
      value = parse(Float64, m.captures[1])/10
  else
    println("could not pick up number from $(DirName)")
  end
  return value
end




# define base directory
base_dir = "./Data/oscillation/"
# get all subdirectories
load_dir_list = get_subdirectories(base_dir)
flip_num = 200

# function main()
  f3 = Figure(resolution = (800, 600))

  dir_fft = Dict()
  dir_fft_all = Dict()
  dir_fft_region = Dict()
  dir_NumInRegion = Dict()
  dir_AlignedVx = Dict()
  dir_yOscillation = Dict()

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
    eta_NumInRegion = Dict()
    eta_AlignedVx = Dict()
    eta_yOscillation = Dict()
    x_length = length(xy_in_time_all(joinpath(base_dir, load_dir, sorted_files[1][2]))[1])
    for eta in 2:10

      # ax3 = Axis(f3[(eta-2)÷3+1, (eta-2)%3],title = "eta = $(round((eta-1)*0.2, sigdigits=2)), frac = $(split(load_dir, "_")[2])", xlabel = "Time step", ylabel = "v_x")
      # limit the range of y-axis from 0 to 1
      # ylims!(high = 1.3)
      # ylims!(low = -1.3)

      sample_num = length(sorted_files)
      time_step = 1:x_length
      x, y = zeros(x_length), zeros(x_length)
      x_region, y_region = zeros(x_length), zeros(x_length)
      x_all, y_all = zeros(x_length), zeros(x_length)
      alined_x = zeros(x_length)
      alined_y = zeros(x_length)
      target_file = sorted_files[1][eta]
      alined_xop, alined_yop, num_in_region_list = xy_in_time_with_aline_x(joinpath(base_dir, load_dir, target_file), spatio_func)

      for sample_file in sorted_files
      # sample_file = sorted_files[1]
        target_file = sample_file[eta]

        all_xop_all, all_yop_all = xy_in_time_all(joinpath(base_dir, load_dir, target_file))
        all_xop_region, all_yop_region = xy_in_time_region(joinpath(base_dir, load_dir, target_file))

        alined_xop, alined_yop, num_in_region_list = xy_in_time_with_aline_x(joinpath(base_dir, load_dir, target_file), spatio_func)

        x_flip_all, y_flip_all = flipx(all_xop_all, flip_num), flipy(all_yop_all, flip_num)
        x_flip_region, y_flip_region = flipx(all_xop_region, flip_num), flipy(all_yop_region, flip_num)

        x_all .+= x_flip_all/sample_num
        x_region .+= x_flip_region/sample_num
        y_all .+= y_flip_all/sample_num
        y_region .+= y_flip_region/sample_num

        alined_x .+= alined_xop/sample_num
        alined_y .+= alined_yop/sample_num
      end

      fft_freq = rfftfreq(x_length, 1/(60*0.001))
      y_fft_all = abs.(rfft(y_all))
      y_fft_region = abs.(rfft(y_region))

      frac = split(load_dir, "_")[2]
      # freq = fft_freq[argmax(y_fft)]
      average_alined_x = mean(alined_x[3000:end])
      # T = 1/(2*pi*freq)
      # que_for_sort = sortperm(y_fft, rev=true)
      # sorted_freq = fft_freq[que_for_sort]
      # eta_num = [(eta-1)*0.2 for _ in 1:3]

      # show the graph
      # eta_fft["$eta"] = [fft_freq, y_fft]
  f = Figure(resolution = (1500, 1000), fontsize = 20)
  Axes_yOscillation = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "v_y", title = "η = $((i+1)/10)", limits = (nothing, (-1.3, 1.3))) for i in 1:7]
  Label(f[0, :], "Oscillation in y direction")
  for (index, (key_dict, value_dict)) in enumerate(dir_yOscillation)
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
  end
  save("/home/muta/Code/slidev/project/b8/last/src/yOscillation.png",f)
      eta_fft_all["$eta"] = [fft_freq, y_fft_all]
      eta_fft_region["$eta"] = [fft_freq, y_fft_region]

      eta_NumInRegion["$eta"] = [1:length(num_in_region_list), num_in_region_list]
      eta_AlignedVx["$eta"] = [1:length(alined_x), alined_x]
      eta_yOscillation["$eta"] = [1:length(y_region), y_region]
    end
    dir_fft[load_dir] = eta_fft
    dir_fft_all[load_dir] = eta_fft_all
    dir_fft_region[load_dir] = eta_fft_region
    dir_NumInRegion[load_dir] = eta_NumInRegion
    dir_AlignedVx[load_dir] = eta_AlignedVx
    dir_yOscillation[load_dir] = eta_yOscillation
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

  f = Figure(resolution = (1000, 800), fontsize = 20)
  ax = Axis(f[1,1], xlabel = "η", ylabel = "frequency")
  maker_size = 20
  handles = [MarkerElement(marker=:circle, color = color, markersize = maker_size) for color in color_list[1:length(top_freq_dir)]]
  labels = string.(PickUpNumFromDirName.(collect(keys(top_freq_dir))))
  legend = Legend(f, handles, labels, "each bandwidth")
  f[1,2] = legend

  Label(f[0, :], "average of top 3 frequency")
  for (index, (key_dict, value_dict)) in enumerate(top_freq_dir)
    for (key_eta, value_eta) in value_dict
      scatter!(ax, parse(Int, key_eta)/10, sum(value_eta[1].*value_eta[2])/sum(value_eta[2]), markersize=20, color=color_list[index], label="bandwidth:$( PickUpNumFromDirName(key_dict))")
    end
  end
  save("/home/muta/Code/slidev/project/b8/last/src/Top_frequency&eta.png", f)

  #plot the graph
  f = Figure(resolution = (1500, 1000), fontsize = 20)
  Axes_fft = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "frequency", ylabel = "amplitude", title = "η = $((i+1)/10)", limits = ((0, 0.12), (0, 160))) for i in 1:9]
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

  f = Figure(resolution = (1500, 1000), fontsize = 20)
  Axes_fft = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "frequency", ylabel = "amplitude", title = "η = $((i+1)/10)", limits = ((0, 0.12), (0, 160))) for i in 1:9]
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

  f = Figure(resolution = (1500, 1000), fontsize = 20)
  Axes_NumInRegion = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Timestep", ylabel = "Number of agent in the band", title = "η = $((i+1)/10)", limits=(nothing, (0, 900))) for i in 1:9]
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

  f = Figure(resolution = (1500, 1000), fontsize = 20)
  Axes_AlignedVx = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "v_x", title = "η = $((i+1)/10)", limits = (nothing, (-1, 1.0))) for i in 1:9]
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
  save("/home/muta/Code/slidev/project/b8/last/src/v_x.png",f)

  f = Figure(resolution = (1500, 1000), fontsize = 20)
  Axes_yOscillation = [Axis(f[((i-1)÷3+1), (i-1)%3*2-1], xlabel = "Time step", ylabel = "v_y", title = "η = $((i+1)/10)", limits = (nothing, (-1.3, 1.3))) for i in 1:7]
  Label(f[0, :], "Oscillation in y direction")
  for (index, (key_dict, value_dict)) in enumerate(dir_yOscillation)
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
  end
  save("/home/muta/Code/slidev/project/b8/last/src/yOscillation.png",f)




# save("/home/muta/slidev/project/b8/discussion/src/frequency_frac04.png", f)
  # save("/home/muta/slidev/project/b8/discussion/src/Top_frequency&eta_frac04.png", f1)
  # save("/home/muta/slidev/project/b8/discussion/src/number_of_agents_frac04.png", f2)
  # save("/home/muta/slidev/project/b8/discussion/src/alined_vx_$(split(load_dir, "_")[2]).png", f3)
  # save("/home/muta/slidev/project/b8/discussion/src/vx_$(split(load_dir, "_")[2]).png", f3)
# end

load_dir = load_dir_list[2]
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
file_list = get_file_list(joinpath(base_dir,load_dir))
file = joinpath(base_dir, load_dir, file_list[10])
idx = findfirst("frac", load_dir)[1]+4
# print(parse(Float64, load_dir[idx:idx+1])/10)
# write_3d_graph(file, spatio_func, 100)
# main()
