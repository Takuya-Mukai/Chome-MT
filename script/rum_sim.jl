include("../src/Chemo-MT.jl")


angle2dir(ϕ) = SA[cos(ϕ), sin(ϕ)]
function vec2comp(vec)
    v1 = [v[1] for v in vec]
    v2 = [v[2] for v in vec]

    return v1, v2
end



function run_Nsample(nsample = 5, ϵ=0.1)
    r = 1.0
    density = 0.1
    A = π*r^2
    # N = round(Int, density*paras.L / A)
    N = 1000
    L = sqrt(N*A / density)

    paras = initParas(Dr=0.01, vel=1.0, J=1.0, epsilon=0.05, cutoff=10, L = L, sigma=r)
    dt = 5e-4
    isave = 400
    n_step = 5_000_000
    for i in 1:nsample 
        system = initSystem(paras, N=N)
        
        println("Relaxation start")
        Relaxation!(system, paras)

        println("Simulation Start")
        @time all_pos, all_ϕ = simulate_midStep(system, paras, nsteps=n_step, dt=dt, isave=isave)
        println("---- Finished ----")

        #* save data
        save_dir = "Data/horizontal_N$(N)_eps$(ϵ)_nstep$(n_step)_dt$(dt)__isave$(100)_sample$(i+2)"
        data = Dict(
            "pos" => all_pos,
            "angle" => all_ϕ
        )
        save(save_dir * ".jld2", data)
        println("Data saved to $save_dir")
    end
end


run_Nsample(1)



function run_eps_list()
    eps_list = (0.1)
    for ϵ in eps_list
        run_Nsample(1, ϵ)
    end
end
# run_eps_list()


# dir = "Data/radial_500_eps0.1_nstep60000_dt0.002_sample19.jld2"
# dir = "Data/horizontal__N500_eps0.01_nstep60000_dt0.002_sample1.jld2"
# dir = "Data/radial_500_eps0.01_nstep100000_dt0.0001__isave100_sample1.jld2"
# data = load(dir)
# all_pos = data["pos"]
# all_ϕ = data["angle"]
# nt = [angle2dir.(ϕs) for ϕs in all_ϕ]

# xs, ys = vec2comp.(all_pos)
# us, vs = vec2comp.(nt)
# savew_dir = join(split(split(dir, "/")[2], ".")[1:end-1], ".")
# XS = 0:L
# YS = 0:L
# # field = [x/L for x in XS, y in YS]
# field = [Gaussian2D(x, y, L / 2, L / 2, L / 4) for x in XS, y in YS]


# @time Animattion(all_pos, nt, field, "Videos/" * savew_dir * ".mp4"; framerate=60, isave=1, L=L)
# @time Animattion(all_pos, nt, "test_radius.mp4"; framerate=600, isave=10, L = L)