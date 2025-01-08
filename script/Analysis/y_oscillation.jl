include("../../src/Chemo-MT.jl")
using CairoMakie

"""
get all file with .jld2
"""
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


"""
Division into sample X
"""
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


function get_sorted_file(base_dir::String, load_dir::String)
    loaddir = base_dir * load_dir
    file_list = get_file_list(loaddir)
    sorted_file = Division_into_smaple(file_list)
    return sorted_file
end




function xy_in_time(file::String)
    # , xmin, xmax, ymin, ymax)
    data = load(file)
    all_ori = data["all_ori"]
    # all_pos = data["all_pos"]

    Nstep = length(all_ori)
    all_yop = zeros(Nstep)
    all_xop = zeros(Nstep)
    for i in 1:Nstep
        # p_inside = all_ori[i][is_inside_region.(all_pos[i], xmin, xmax, ymin, ymax)]
        # all_xop[i], all_yop[i] = mean(angle2dir.(p_inside))
        all_xop[i], all_yop[i] = mean(angle2dir.(all_ori[i]))
    end

    return all_xop, all_yop
end


function spatio_func(x, y, L)
    L_frac = 0.2
    up = L / 2 + L * L_frac / 2
    low = L / 2 - L * L_frac / 2

    if y > low && y < up
        return 1.0
    else
        return 0.0
    end
end


function smoothing(x::Vector, n::Int64)
  # make smoothing data
  N = length(x)
  x_smooth = zeros(N-n+1)
  # change method if n is even or not
  if n % 2 == 1
    for i in 1:N-n+1
      if (i == 1)
        x_smooth[i] = sum(x[1:n]) / n
      else
        x_smooth[i] = x_smooth[i-1] + (x[i+n-1] - x[i-1]) / n
      end
    end
  else
    error("make smoothing parameter to odd number")
  end
  return x_smooth
end


function xy_in_time_with_smoothing(file_dir::String; smoothing_parameter=11)
  x, y = xy_in_time(file_dir)
  x_smooth = smoothing(x, smoothing_parameter)
  y_smooth = smoothing(y, smoothing_parameter)
  return x_smooth, y_smooth
end


# main
L = 30.0
L_frac = 0.2
up = L / 2 + L * L_frac / 2
low = L / 2 - L * L_frac / 2

# define file to pick up data
base_dir = "./Data/oscillation/"
load_dir = "ystep_Lfrac02_L30_nstep_2M_v1/"
sorted_file = get_sorted_file(base_dir, load_dir)
target_file = sorted_file[1][4]
file_path = base_dir * load_dir * target_file

# make data to write graph from file_path
# s_p means smoothing parameter
s_p = 11
x, y = xy_in_time_with_smoothing(file_path, smoothing_parameter = s_p)
step_list = [i for i in 1:length(x)]

# make graph from data
f = Figure()
ax = Axis(f[1, 1],
          xlabel = "Time step",
          ylabel = "Mean of all agents",
          title = "Oscillation in y and x direction",)
lines!(ax, step_list, x, label = "x direction")
lines!(ax, step_list, y, label = "y direction")
axislegend(ax)
pic_name = "Data/oscillation/pic/oscillation_"* reduce(*, split(target_file, ".jld2")) *"_$(smoothing_parameter).png"
save(pic_name, f)
println("saved graph to $pic_name")


# make animation
# data = load(file_path)
# all_ori = data["all_ori"]
# all_pos = data["all_pos"]
# @time Animattion(all_pos, all_ori, spatio_func,join(split(loaddir * target_file, ".")[1:end-1], ".") * ".mp4"; framerate=60, L=30)
