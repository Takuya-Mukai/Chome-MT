include("y_oscillation.jl")



function get_sum_for_each_frac(base_dir, file_list::Vector)
  # make sure the vector length is the same for all files
  vector_length = length(load(joinpath(base_dir, file_list[1]))["data"]["x"])
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
    data = load(joinpath(base_dir, file))["data"]
    dict[frac][1] = dict[frac][1] .+ data["x"]
    dict[frac][2] = dict[frac][2] .+ data["y"]
  end
  return dict 
end


base_dir = "./Data/oscillation/smoothed"
file_list = readdir(base_dir)
data = get_sum_for_each_frac(base_dir, file_list)


smoothing_num = 1
# write the graph by using data with CairoMakie
for keys in keys(data)
  x_raw = data[keys][1]
  y_raw = data[keys][2]
  x = smoothing(x_raw, smoothing_num)
  y = smoothing(y_raw, smoothing_num)
  step_list = [i for i in 1:length(x)]
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
  save(joinpath("Data/oscillation/pic", "sum_$(keys)_smoothing$(smoothing_num).png"), fig)
  println("Sum of $(keys) is saved.")
end
