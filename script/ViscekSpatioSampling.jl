include("../src/Chemo-MT.jl")


rho = 1.0
L = 30
j = 1.0
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

function spatio_func(x,y,L)
    L_frac = 0.6
    up = L/2 + L*L_frac/2
    low = L/2 - L*L_frac/2
    
    if y > low && y < up
        return 1.0
    else
        return 0.0
    end
end

# spatio_func(x,y,L) = 1.0
# const spatio_func_string = " if r <= L/4,  j =1.0, else j = 0.0"
const spatio_func_string = "function spatio_func(x,y,L)
    L_frac = 0.6
    up = L/2 + L*L_frac/2
    low = L/2 - L*L_frac/2
    
    if y > low && y < up
        return 1.0
    else
        return 0.0
    end 
end"

Jmax = maximum(spatio_func.(0:L, 0:L, paras.L))* π
Jmin = minimum(spatio_func.(0:L, 0:L, paras.L)) * π

eta_list = 0.0:0.2: (Jmax + 2.0)
# eta_list = [0.0]
N_sample = 20

# save_root = "/media/kylechan/Extreme SSD/Kyle"
# save_root = "/mnt/SSD_RAID0/kyle/viscek_spatio/" 
save_root = "Data/"
# save_root = "D:/kylechan/tmp_data/"
# save_root = "E:/Kyle/Simulation_Data/"

fname = "ystep_Lfrac06_L30_v3/"
# fname = "Const_j1_L30_v3/"

save_dir = save_root*fname
function phase_transition_spatio(eta_list, paras, spatio_func, p; save_dir="")
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

function pt_sampling(N_sample, eta_list; save_dir="")
    println("---- Simulation Started ----")
    # L = 20
    # Jmax = maximum(spatio_func.(0:L, 0:L)) * π
    # Jmin = minimum(spatio_func.(0:L, 0:L)) * π

    # eta_list = 0.0:0.2:(Jmax+2.0)
    p = Progress(N_sample*length(eta_list); dt=1.0)
    for Nsamp in 1:N_sample
        save_dir_ = save_dir*"sample$(Nsamp)_"
        @show save_dir_
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


pt_sampling(N_sample, eta_list; save_dir=save_dir)

println("finished all data saved to $(save_dir)")