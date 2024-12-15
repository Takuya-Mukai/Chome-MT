include("../src/Chemo-MT.jl")
using Base.Threads

"""
rho: density of particle
L: length
j: strength of interaction
N: number of particle
"""
rho = 1.0
L = 10
j = 2.0
eta = 1.0
N = round(Int, L^2*rho)
# alpha = eta - J*rho
# alpha = -1
# j = 
# j = (eta - alpha) / rho
# J = j / π
# # @show alpha

paras = initParas(
                Dr = eta,
                vel=1.0,
                J=j,
                epsilon=0.0,
                cutoff=2,
                L=L,
                dx=0.1,
                sigma=1)



# function spatio_func(x,y, L)
#     c0 = SA[L/2, L/2]
#     r = SA[x,y]
#     dr = norm(c0 - r)
    
#     if dr <= L/4
#         return 1.0
#     else
#         return 0.0
#     end
# end

# function spatio_func(x,y,L)
#     L_frac = 0.6
#     up = L/2 + L*L_frac/2
#     low = L/2 - L*L_frac/2
#
#     if y > low && y < up
#         return 1.0
#     else
#         return 0.0
#     end
# end

# distribution of interaction
spatio_func(x,y,L) = 1.0
# const spatio_func_string = " if r <= L/4,  j =1.0, else j = 0.0"
const spatio_func_string = "spatio_func(x,y,L) = 1.0"

Jmax = maximum(spatio_func.(0:L, 0:L, paras.L))* π
Jmin = minimum(spatio_func.(0:L, 0:L, paras.L)) * π

# list of time to make data
eta_list = 0.0:0.2: (Jmax + 2.0)
# eta_list = [0.0]
N_sample = 5

# save_root = "/media/kylechan/Extreme SSD/Kyle"
# save_root = "/mnt/SSD_RAID0/kyle/viscek_spatio/" 
# save_root = "D:/kylechan/tmp_data/"
# save_root = "E:/Kyle/Simulation_Data/"
# fname = "Const_j1_L30_v3/"


function phase_transition_spatio(eta_list, paras, spatio_func, p; save_dir="")
    #* check if folder exist
    # op_list = similar(eta_list)
    LEN = length(eta_list)
    @show LEN
    # simulate and save data
    for i in 1:LEN
        eta = round(eta_list[i], digits=3)
        println("")
        @show eta
        paras.Dr = eta
        system = initSystem(paras, N=N)
        

        # make data
        @time all_pos, all_ori = simulation_visceks_spatio(system, paras, spatio_func;
            nsteps=50_000,
            dt=0.001,
            isave=60)
        
        # save data
        if length(save_dir) != 0
            # data to save
            data = Dict("all_pos"=>all_pos, "all_ori"=>all_ori, "paras"=>paras)
            @time save(save_dir*"eta_$(eta).jld2", data)

            obj2txt(paras, save_dir*"eta_$(eta)_input.txt")
            fio = open(save_dir*"eta_$(eta)_input.txt", "a")
            write(fio, spatio_func_string)

            # save pitcture of each background field
            @time save_bg_field(paras, spatio_func, save_dir)
        end

        next!(p)
        # op = get_order_parameter(all_ori)
        # println("ordre parameter = $op")
    end
end

function pt_sampling(N_sample, paras, eta_list; save_dir="")
    println("---- Simulation Started ----")
    # L = 20
    # Jmax = maximum(spatio_func.(0:L, 0:L)) * π
    # Jmin = minimum(spatio_func.(0:L, 0:L)) * π

    # eta_list = 0.0:0.2:(Jmax+2.0)
    p = Progress(N_sample*length(eta_list); dt=1.0)
    for Nsamp in 1:N_sample
        save_dir_ = save_dir*"sample$(Nsamp)_"
        # @show save_dir_regularity
        @time phase_transition_spatio(eta_list, paras, spatio_func, p; save_dir=save_dir_)
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

# Data directory
save_root = "Data/basic_Vicsek/"
savedir_list = String[]  # Initializing an empty array of Strings
j_list = []

# make list for rho
for i in 1:10
    local fname = "constant_L10_ρ$(i)/"
    local j = i*1.0
    push!(savedir_list, save_root * fname)  # Adding the constructed filename to the list
    push!(j_list, j)  # Adding the rho value to the list
end

# iterate for each ρ
println("Number of threads: $(Threads.nthreads())")
@threads for i in 1:10
    local save_dir = savedir_list[i]
    local copy_paras = deepcopy(paras)
    copy_paras.J = j_list[i]
    local field_strength = round(Int, L^2*rho)
    pt_sampling(N_sample, copy_paras, eta_list; save_dir=save_dir)
    println("finished all data saved to $(save_dir)")
end
