include("../src/Chemo-MT.jl")


angle2dir(ϕ) = SA[cos(ϕ), sin(ϕ)]
function vec2comp(vec)
    v1 = [v[1] for v in vec]
    v2 = [v[2] for v in vec]

    return v1, v2
end




r = 1.0
density = 0.05
A = π*r^2
# N = round(Int, density*paras.L / A)
N = 500

L = sqrt(N / density)

# paras = initParas(Dr=0.1, vel=1.0, J=1.0, epsilon=0.01, cutoff=10, L = L, sigma=r)

# system = initSystem(paras,N = N)
# Relaxation!(system, paras)
# @time all_pos, all_ϕ = simulate(system, paras, nsteps=20_000, dt=2e-3)
# println("---- Finished ----")


# dir = "Data/radial_500_eps0.1_nstep60000_dt0.002_sample19.jld2"
dir = "Data/horizontal__N500_eps0.01_nstep60000_dt0.002_sample1.jld2"

data = load(dir)
all_pos = data["pos"]
all_ϕ = data["angle"]
nt = [angle2dir.(ϕs) for ϕs in all_ϕ]

xs, ys = vec2comp.(all_pos)
us, vs = vec2comp.(nt)
savew_dir = join(split(split(dir, "/")[2], ".")[1:end-1], ".")
XS = 0:L
YS =0:L
field = [x/L for x in XS, y in YS]
# field = [Gaussian2D(x, y, L/2, L/2, L/4) for x in XS, y in YS]


@time Animattion(all_pos, nt, field, "Videos/"*savew_dir*".mp4"; framerate=600, isave=10, L=L)
# @time Animattion(all_pos, nt, "test_radius.mp4"; framerate=600, isave=10, L = L)

function run_Nsample(nsample = 5, ϵ=0.1)
    r = 1.0
    density = 0.05
    A = π*r^2
    # N = round(Int, density*paras.L / A)
    N = 500

    L = sqrt(N / density)

    paras = initParas(Dr=0.1, vel=1.0, J=1.0, epsilon=ϵ, cutoff=10, L=L, sigma=r)
    dt = 2e-3
    for i in 1:nsample 
        system = initSystem(paras, N=N)
        # Relaxation!(system, paras)
        @time all_pos, all_ϕ = simulate(system, paras, nsteps=60_000, dt=dt)
        

        #* save data
        save_dir = "Data/radial_$(N)_eps$(ϵ)_nstep$(60_000)_dt$(dt)_sample$(i)"
        data = Dict(
            "pos" => all_pos,
            "angle" => all_ϕ
        )
        save(save_dir * ".jld2", data)
    end
end


# run_Nsample(20)



function run_eps_list()
    eps_list = (1.0, 0.1, 0.01)
    for ϵ in eps_list
        run_Nsample(20, ϵ)
    end
end
# run_eps_list()