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


function main()
  # make graph for each directory
  for (index, load_dir) in enumerate(load_dir_list)
    # define file to pick up data
    sorted_file = get_sorted_file(base_dir; load_dir)
    sample_num = length(sorted_file)

      for i in 1:sample_num
        target_file = sorted_file[i][4]
        file_path = base_dir * load_dir * target_file

        # pick up data for animation
        data = load(file_path)
        all_ori = data["all_ori"]
        all_pos = data["all_pos"]

        # make data to write graph from file_path
        x_raw, y_raw = xy_in_time(file_path)
        x_raw, y_raw = flipx(x_raw, flip_num), flipy(y_raw, flip_num)
        x = smoothing(x_raw, s_p)
        y = smoothing(y_raw, s_p)
        x, y = flipx(x, flip_num), flipy(y, flip_num)
        step_list = [i for i in 1:length(x)]

        # write x, y to file with jld2
        save_file_path = "Data/oscillation/smoothed/oscillation_$(split(load_dir, "_")[2])_$(target_file)"
        @save save_file_path x_raw y_raw
        println("saved data to $save_file_path")

        # make animation from data
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


        # save_animation_dir = joinpath("Data/oscillation/pic/animation/", load_dir)
        # animation_file = joinpath( target_file*".mp4")
        # # check if the directory is already exist
        # if isdir(save_animation_dir) == false
        #   mkdir(save_animation_dir)
        # end
        # @time Animattion(all_pos, all_ori, spatio_func, joinpath(save_animation_dir, animation_file); framerate=60, L=30)
        # println("saved animation to $(joinpath(save_animation_dir, animation_file))")



        # make graph from data
        f = Figure()
        ax = Axis(f[1, 1],
                xlabel = "Time step",
                ylabel = "Mean of all agents",
                title = "Oscillation in y direction",)
        # lines!(ax, step_list, x, label = "x direction")
        lines!(ax, step_list, y, label="y direction $(s_p)")
        axislegend(ax)
        pic_name = "Data/oscillation/pic/y_oscillation_$(split(load_dir, "_")[2])_"*
                          reduce(*, split(target_file, ".jld2")) *"_smoothing$(s_p).png"
        save(pic_name, f)
        println("saved graph to $pic_name")

      end
  end
end


# main
L = 30.0
L_frac = 0.2
up = L / 2 + L * L_frac / 2
low = L / 2 - L * L_frac / 2

# define base directory
base_dir = "./Data/oscillation/"
# get all subdirectories
load_dir_list = get_subdirectories(base_dir)

# s_p means smoothing parameter
s_p = 41
# flip_num means the length of the sum of x and y
flip_num = 200


# main()
