include("../src/Chemo-MT.jl")

rho = 1
L = 20
j = 1.0
eta = 0.0
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
                J=J, 
                epsilon=0.0, 
                cutoff=2, 
                L=L,
                dx=0.1, 
                sigma=1)

system = initSystem(paras, N=N)

# spatio_func(x, y) = sin(y * π / L) + 0.5
spatio_func(x, y) = + 1

# alpha_y(eta, rho, x,y) = eta - rho*j*spatio_func(x,y)/π

@time all_pos, all_ori = simulation_visceks_spatio(system, paras, spatio_func;
                                                    nsteps=50_000,
                                                    dt=0.001,
                                                    isave=60)



@time Animattion(all_pos, all_ori, "Videos/Viscek/positiveJ" * ".mp4"; framerate=60, L=paras.L)




function phase_transition_spatio(eta_list, spatio_func)

    op_list = similar(eta_list)
    LEN = length(eta_list)
    @show LEN
    for i in 1: LEN

        eta = eta_list[i]
        @show eta
        # @show α
        rho = 1
        L = 20
        J = 1
        # eta = 1
        N = round(Int, L^2 * rho)
        # j = (eta - α) / rho
        # J = j / π
        # alpha = eta - J * rho

        paras = initParas(
            Dr=eta,
            vel=1.0,
            J=J,
            epsilon=0.0,
            cutoff=2,
            L=L,
            dx=0.1,
            sigma=1)

        system = initSystem(paras, N=N)

        @time all_pos, all_ori = simulation_visceks_spatio(system, paras, spatio_func;
            nsteps=50_000,
            dt=0.001,
            isave=60)
        op = get_order_parameter(all_ori)
        println("ordre parameter = $op")


        op_list[i] = op
    end

    return op_list
end


# spatio_func(x, y) = -1

# Jmax = maximum(spatio_func.(0:L, 0:L))*π
# Jmin = minimum(spatio_func.(0:L, 0:L))*π

# eta_list = 0.0:0.2: (Jmax + 2.0)


# @time op_list = phase_transition_spatio(eta_list, spatio_func)

# data = Dict("eta"=>eta_list, "op"=> op_list)
# save("Data/viscek/phase_transition_spatio_v3.jld2", data)

"""
sampling
"""
# N_sample = 20

# function pt_sampling(N_sample)
#     L = 20
#     Jmax = maximum(spatio_func.(0:L, 0:L))*π
#     Jmin = minimum(spatio_func.(0:L, 0:L))*π

#     eta_list = 0.0:0.2: (Jmax + 2.0)
#     p = Progress(N_sample; dt=5.0)
#     for i in 1:N_sample
#         op_list = phase_transition_spatio(eta_list, spatio_func)
#         data = Dict("op" => op_list, "eta" => eta_list)
#         save("Data/viscek_spatio/phase_transistion_spatio_05_sample_$i.jld2", data)

#         next!(p)
#     end
# end

# @time @time pt_sampling(N_sample)


"""
plotting
"""
# XS = 0:0.2:L
# YS = 0:0.2:L
# bg = [spatio_func(x,y) for x in XS, y in YS]
# # heatmap(bg)