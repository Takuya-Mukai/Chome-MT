include("../src/Chemo-MT.jl")

rho = 1
L = 50
J = 0.9
eta = 1
N = round(Int, L^2*rho)
alpha = eta - J*rho

@show alpha

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

@time all_pos, all_ori = simulation_visceks(system, paras; 
                                                    nsteps=50_000,
                                                    dt=0.001,
                                                    isave=60)



@time Animattion(all_pos, all_ori, "Videos/Viscek/test_v1" * ".mp4"; framerate=60, L=paras.L)

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


alpha_list = range(-2,1,40)
op_list = phase_transition(alpha_list)
data = Dict("op"=>op_list, "alpha" => alpha_list)
save("Data/viscek/phase_transistion_v2.jld2", data)


