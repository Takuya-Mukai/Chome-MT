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


function xy_in_time(file::String)
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
    all_xop[i] = sum( spatio_func.(x, y, data["paras"].L) .* aline.(cos.(all_ori[i])) )/num_in_region
    all_yop[i] = sum( spatio_func.(x, y, data["paras"].L) .* aline.(sin.(all_ori[i])) )/num_in_region
    # all_xop[i] = sum( cos.(all_ori[i]) )/num_in_region
    # all_yop[i] = sum( sin.(all_ori[i]) )/num_in_region
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



# define base directory
base_dir = "./Data/oscillation/"
# get all subdirectories
load_dir_list = get_subdirectories(base_dir)
flip_num = 200

function main()
  # f = Figure(resolution = (800, 600))
  # f1 = Figure(resolution = (800, 600))
  # f2 = Figure(resolution = (800, 600))
  f3 = Figure(resolution = (800, 600))
  # define color list for graph
  color_list = ["red", "blue", "green", "orange", "purple", "black", "brown", "pink", "gray", "cyan", "magenta", "yellow"]

  for (color, load_dir) in enumerate(load_dir_list)
  # color = 1
  # i = 1
  # load_dir = load_dir_list[i]
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

    # ax1 = Axis(f1[1, 1],title = "frac = $(split(load_dir, "_")[2]), frac = $(split(load_dir, "_")[2])", xlabel = "frequency", ylabel = "amplitude")
    for eta in 2:8
      # eta = 3
      # ax = Axis(f[(eta-2)÷3+1, (eta-2)%3],title = "eta = $(round((eta-1)*0.2, sigdigits=2)), frac = $(split(load_dir, "_")[2])", xlabel = "frequency", ylabel = "amplitude")
      # ylims!(high = 55)

      # ax2 = Axis(f2[(eta-2)÷3+1, (eta-2)%3],title = "eta = $(round((eta-1)*0.2, sigdigits=2)), frac = $(split(load_dir, "_")[2])", xlabel = "Timestep", ylabel = "Number of agent in the band")

      ax3 = Axis(f3[(eta-2)÷3+1, (eta-2)%3],title = "eta = $(round((eta-1)*0.2, sigdigits=2)), frac = $(split(load_dir, "_")[2])", xlabel = "Time step", ylabel = "v_x")
      # limit the range of y-axis from 0 to 1
      ylims!(high = 1.3)
      ylims!(low = -1.3)

      x_length = length(xy_in_time(joinpath(base_dir, load_dir, sorted_files[1][eta]))[1])
      sample_num = length(sorted_files)
      time_step = 1:x_length
      x, y = zeros(x_length), zeros(x_length)
      alined_x = zeros(x_length)
      alined_y = zeros(x_length)

      # for sample_file in sorted_files
      sample_file = sorted_files[1]
        target_file = sample_file[eta]
        all_xop, all_yop = xy_in_time(joinpath(base_dir, load_dir, target_file))
        alined_xop, alined_yop, num_in_region_list = xy_in_time_with_aline_x(joinpath(base_dir, load_dir, target_file), spatio_func)
        x_flip, y_flip = flipx(all_xop, flip_num), flipy(all_yop, flip_num)
        x .+= x_flip/sample_num
        y .+= y_flip/sample_num
        alined_x .+= alined_xop
        alined_y .+= alined_yop
        # lines!(time_step, y_flip)
      # end
      fft_freq = rfftfreq(x_length, 1/(60*0.001))
      y_fft = abs.(rfft(y))
      frac = split(load_dir, "_")[2]
      freq = fft_freq[argmax(y_fft)]
      average_alined_x = mean(alined_x[3000:end])
      T = 1/(2*pi*freq)
      que_for_sort = sortperm(y_fft, rev=true)
      sorted_freq = fft_freq[que_for_sort]
      eta_num = [(eta-1)*0.2 for _ in 1:3]
      # println("| $(round((eta-1)*0.2, sigdigits=2)) | $(round(sorted_freq[1], sigdigits=3)) | $(round(sorted_freq[2], sigdigits=3)) | $(round(sorted_freq[3], sigdigits=3)) |")
      # println("frac: $(frac), η: $(round(eta, sigdigits=3)), freq: $(round(freq, sigdigits=3)), T= $(round(T, sigdigits=3)), average_alined_x: $(round(average_alined_x, sigdigits=3)), wave_length: $(round(T*average_alined_x, sigdigits=3))")
      # println("| $(frac) | $(round((eta-1)*0.2, sigdigits=3)) | $(round(freq, sigdigits=3)) | $(round(T, sigdigits=3)) | $(round(average_alined_x, sigdigits=3)) | $(round(T*average_alined_x, sigdigits=3)) |")

      # lines!(ax, fft_freq[1:25], y_fft[1:25], color = color_list[color], label = frac)
      # scatter!(ax1, eta_num, sorted_freq[1:3], color = color_list[eta-1])
      # lines!(ax2, 1:length(num_in_region_list), num_in_region_list, label = "eta = $(eta_num[1])")
      # lines!(ax3, time_step, alined_x, color = color_list[color], label = frac*" alined v_x")
      lines!(ax3, time_step, alined_x, color = color_list[color], label = frac*" v_x")
      # show the graph
    end
    display(f3)
  end
  # save("/home/muta/slidev/project/b8/discussion/src/frequency_frac04.png", f)
  # save("/home/muta/slidev/project/b8/discussion/src/Top_frequency&eta_frac04.png", f1)
  # save("/home/muta/slidev/project/b8/discussion/src/number_of_agents_frac04.png", f2)
  # save("/home/muta/slidev/project/b8/discussion/src/alined_vx_$(split(load_dir, "_")[2]).png", f3)
  # save("/home/muta/slidev/project/b8/discussion/src/vx_$(split(load_dir, "_")[2]).png", f3)
end

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
main()
