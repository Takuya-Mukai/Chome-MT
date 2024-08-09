include("../src/Chemo-MT.jl")

# rho = 1
# L = 10
# J = 0.9
# eta = 1
# N = round(Int, L^2*rho)
# alpha = eta - J*rho

# @show alpha

# paras = initParas(
#                 Dr = eta, 
#                 vel=1.0, 
#                 J=J, 
#                 epsilon=0.0, 
#                 cutoff=2, 
#                 L=L,
#                 dx=0.1, 
#                 sigma=1)

# system = initSystem(paras, N=N)

# @time all_pos, all_ori = simulation_visceks(system, paras; 
#                                                     nsteps=50_000,
#                                                     dt=0.001,
#                                                     isave=60)



# @time Animattion(all_pos, all_ori, "Videos/Viscek/test_v1" * ".mp4"; framerate=60, L=paras.L)

function get_order_parameter(all_ori)
    # sample = 1:100:10_000
    # n = length(sample)
    op_list = zeros(100)
    N = length(all_ori)
    for i in 1:100
        op_list[i] = norm(mean(angle2dir.(all_ori[N-i])))
    end
    return mean(op_list)
end

function phase_transition(alpha_list)

    op_list = similar(alpha_list)
    for i in 1:length(alpha_list)

        α = alpha_list[i]
        @show α
        rho = 1
        L = 20
        j = 1
        # eta = 
        N = round(Int, L^2 * rho)
        @show N
        # j = (eta - α) / rho
        eta = α + j*rho
        J = j/π
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

        @time all_pos, all_ori = simulation_visceks(system, paras;
            nsteps=50_000,
            dt=0.001,
            isave=60)
        op = get_order_parameter(all_ori)
        println("ordre parameter = $op")


        op_list[i] = op
    end

    return op_list
end


# alpha_list = range(-2,1,40)
# op_list = phase_transition(alpha_list)
# data = Dict("op"=>op_list, "alpha" => alpha_list)
# save("Data/viscek/phase_transistion_v2.jld2", data)



"""
sampling
"""

N_sample = 20

function pt_sampling(N_sample)
    for i in 1:N_sample
        alpha_list = range(-1, 1, 40)
        op_list = phase_transition(alpha_list)
        data = Dict("op" => op_list, "alpha" => alpha_list)
        save("Data/viscek/phase_transistion_v5_sample_$i.jld2", data)
    end
end

@time pt_sampling(N_sample)

data_dir = "Data/viscek/phase_transistion_v5_sample"
function out_put_order(data_dir)
    data = load(data_dir * "_1.jld2")
    alpha = data["alpha"]
    ops = [data["op"]]
    for i in 2:10
        dir = data_dir*"_$i.jld2"
        data = load(data_dir * "_$i.jld2")
        op = data["op"]
        push!(ops, op)
    end
    return mean(ops)
end

op = out_put_order(data_dir)
fig = Figure()
ax = Axis(fig[1,1], xlabel="α", ylabel="|<v>|")
scatter!(alpha,op)
save("phase_transotion_v5.png", fig)

# scatter(alpha_list, op_list)

# println("ordre parameter = $(get_order_parameter(all_ori))")
"""
plot phase phase_transition
"""

# load("Data/viscek/phase_transistion_v2.jld2")
# op = data["op"]
# alpha = data["alpha"]

# fig = Figure()
# ax = Axis(fig[1,1], xlabel="α", ylabel="|<v>|")
# scatter!(alpha,op)
# save("phase_transotion_v2.png", fig)

