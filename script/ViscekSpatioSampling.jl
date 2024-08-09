include("../src/Chemo-MT.jl")


rho = 1
L = 5
j = 1.0
eta = 1.0
N = round(Int, L^2*rho)
# alpha = eta - J*rho
# alpha = -1
# j = 
# j = (eta - alpha) / rho
J = j / π
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

            
spatio_func(x, y) = sin(y * π / L) + 1.0

Jmax = maximum(spatio_func.(0:L, 0:L))*π
Jmin = minimum(spatio_func.(0:L, 0:L))*π

eta_list = 0.0:0.2: (Jmax + 2.0)
N_sample = 2

# save_root = "/mnt/SSD_RAID0/kyle/viscek_spatio/" 
save_root = "Data/test/"
fname = "VS_sin1_v1_"

save_dir = save_root*fname
function phase_transition_spatio(eta_list, paras, spatio_func; save_dir="")
    #* check if folder exist
    # op_list = similar(eta_list)
    LEN = length(eta_list)
    @show LEN
    for i in 1:LEN
        eta = eta_list[i]
        @show eta
        paras.Dr = eta
        system = initSystem(paras, N=N)

        all_pos, all_ori = simulation_visceks_spatio(system, paras, spatio_func;
            nsteps=50_000,
            dt=0.001,
            isave=60)
        
        if length(save_dir) != 0
            data = Dict("all_pos"=>all_pos, "all_ori"=>all_ori)
            save(save_dir*"sample$i.jld2", data)
        end
        # op = get_order_parameter(all_ori)
        # println("ordre parameter = $op")
    end
end

function pt_sampling(N_sample, eta_list; save_dir="")
    # L = 20
    # Jmax = maximum(spatio_func.(0:L, 0:L)) * π
    # Jmin = minimum(spatio_func.(0:L, 0:L)) * π

    # eta_list = 0.0:0.2:(Jmax+2.0)
    p = Progress(N_sample; dt=1.0)
    for i in 1:N_sample
        @time phase_transition_spatio(eta_list, paras, spatio_func; save_dir=save_dir)
        # data = Dict("op" => op_list, "eta" => eta_list)
        # save("Data/viscek_spatio/phase_transistion_spatio_05_sample_$i.jld2", data)
        next!(p)
    end
end

pt_sampling(N_sample, eta_list; save_dir=save_dir)

