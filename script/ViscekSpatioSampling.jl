include("../src/Chemo-MT.jl")


# Extract the code of a function from a file
function extract_function_code(function_name, file_path)
  # check if the specified file exists
  if !isfile(file_path)
    println("file $file_path does not exist.")
    # return nothing
  end

  # open the file and read its contents
  content = readlines(file_path)

  # a flag to track if we are inside the function
  inside_function = false
  function_code = String[]  # to store the function code

  # iterate through each line in the file
  for line in content
    # check for function definition
    if occursin("#", line)
      continue
    elseif occursin("function $function_name", line)
      inside_function = true
      push!(function_code, line)  # add the function definition line
    elseif inside_function
      # add the function code until we encounter the end of the function
      push!(function_code, line)
      # check for the end of the function (i.e., 'end')
      if occursin("end", line)
        break
      end
    end
  end
  return join(function_code, "\n")
end


function phase_transition_spatio(eta_list, paras::Paras, spatio_func::Function, spatio_func_string::String, p; save_dir="")
    #* check if folder exist
    # op_list = similar(eta_list)
    LEN = length(eta_list)
    @show LEN
    for i in 1:LEN
        eta = round(eta_list[i], digits=3)
        println("")
        @show eta
        paras.Dr = eta
        system = initSystem(paras, N=N)
        

        @time all_pos, all_ori = simulation_visceks_spatio(system, paras, spatio_func;
            nsteps=50_000,
            dt=0.001,
            isave=60)
        
        if length(save_dir) != 0
            data = Dict("all_pos"=>all_pos, "all_ori"=>all_ori, "paras"=>paras)
            @time save(save_dir*"eta_$(eta).jld2", data)

            obj2txt(paras, save_dir*"eta_$(eta)_input.txt")
            fio = open(save_dir*"eta_$(eta)_input.txt", "a")
            write(fio, spatio_func_string)

            @time save_bg_field(paras, spatio_func, save_dir)
        end

        next!(p)
        # op = get_order_parameter(all_ori)
        # println("ordre parameter = $op")
    end
end

function pt_sampling(N_sample::Int32, eta_list::Vector, spatio_func::Function, spatio_func_string::String; save_dir::String="")
    println("---- Simulation Started ----")

    p = Progress(N_sample*length(eta_list); dt=1.0)
    Threads.@threads for Nsamp in 1:N_sample
        local save_dir_ = save_dir*"sample$(Nsamp)_"
        @show save_dir_
        @time phase_transition_spatio(eta_list, paras, spatio_func, spatio_func_string, p; save_dir=save_dir_)
        # data = Dict("op" => op_list, "eta" => eta_list)
        # save("Data/viscek_spatio/phase_transistion_spatio_05_sample_$i.jld2", data)
        # next!(p)
    end
end

function save_bg_field(paras, spatio_func, save_dir)
    xs = 0:0.2:paras.L
    ys = 0:0.2:paras.L
    bg_field = [spatio_func(x, y, paras.L) for x in xs, y in ys]

    fig = Figure()

    ax, hm = heatmap(fig[1, 1][1, 1], xs, ys, bg_field)
    Colorbar(fig[1, 1][1, 2], hm)

    ax.aspect = DataAspect()
    save(save_dir*"backend_field.png",fig)
end

# Main

#=
rho: density of the particles
L: size of the box
j: interaction strength
eta: noise strength
N: number of particles
=#

rho = 1.0
L = 100
j = 1.0
eta = 1.0
N = round(Int, L^2*rho)

paras = initParas(
    Dr = eta, 
    vel=1.0, 
    J=j, 
    epsilon=0.0, 
    cutoff=2, 
    L=L,
    dx=0.1, 
    sigma=1
)


function f10(x,y,L)
  return 1.0
end

function f8(x,y,L)
  if y < L*9/10 && y > L*1/10
    return 1.0
  else
    return 0.0
  end
end

function f6(x,y,L)
  if y < L*8/10 && y > L*2/10
    return 1.0
  else
    return 0.0
  end
end

function f4(x,y,L)
  if y < L*7/10 && y > L*3/10
    return 1.0
  else
    return 0.0
  end
end

function f2(x,y,L)
  if y < L*6/10 && y > L*4/10
    return 1.0
  else
    return 0.0
  end
end

function f0(x,y,L)
  return 0.0
end

spatio_func_list = [f0, f2, f4, f6, f8, f10]
spatio_func_list_str = ["f0", "f2", "f4", "f6", "f8", "f10"]

const spatio_func_code = [extract_function_code(spatio_func_list_str[i], "script/ViscekSpatioSampling.jl")
                          for i in 1:6]

save_root = "Data/"
fname_list = ["f0", "f2", "f4", "f6", "f8", "f10"]
eta_list = collect(0.0:0.2: (Jmax + 2.0))
N_sample::Int32 = 24
for i in 1:6
  
  local save_dir = save_root*fname_list[i]*"/"
  Jmax = maximum(spatio_func_list[i].(0:L, 0:L, paras.L))* π
  Jmin = minimum(spatio_func_list[i].(0:L, 0:L, paras.L)) * π

  pt_sampling(N_sample, eta_list, spatio_func_list[i], spatio_func_list_str[i]; save_dir=save_dir)

  println("finished all data saved to $(save_dir)")
end
