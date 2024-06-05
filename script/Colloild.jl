include("../src/Chemo-MT.jl")

paras = initParas(
                Dr=1/500, 
                vel=1.0, 
                J=1.0, 
                epsilon=1., 
                cutoff=2, 
                L=50,
                dx=0.1, 
                sigma=1)

system = initSystem(paras, N=50)

@time all_pos, all_ori, all_field = simulation_chemotaxis(system, paras; 
                                                    nsteps=100_000,
                                                    dt=0.001,
                                                    isave=60)



@time Animattion(all_pos, all_ori, all_field, "Videos/test/colloid" * ".mp4"; framerate=30, L=paras.L)
