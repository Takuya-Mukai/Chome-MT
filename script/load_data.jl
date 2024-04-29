include("../src/Chemo-MT.jl")

angle2dir(ϕ) = SA[cos(ϕ), sin(ϕ)]
function vec2comp(vec)
    v1 = [v[1] for v in vec]
    v2 = [v[2] for v in vec]

    return v1, v2
end

dir = "Data/horizontal_N1000_eps0.1_nstep5000000_dt0.0005__isave100_sample3.jld2"

data = load(dir)
all_pos = data["pos"]
all_ϕ = data["angle"]
nt = [angle2dir.(ϕs) for ϕs in all_ϕ]

xs, ys = vec2comp.(all_pos)
us, vs = vec2comp.(nt)
savew_dir = join(split(split(dir, "/")[2], ".")[1:end-1], ".")
XS = 0:L
YS = 0:L

# field = [Gaussian2D(x, y, L / 2, L / 2, L / 4) for x in XS, y in YS]
# field = [(.517 - Gaussian1D(y, sigma=L / 4, mu=L / 2) * 100)  for x in XS, y in YS]
field = [( .05 + Gaussian1D(y, sigma=L / 4, mu=L / 2) * 100) for x in XS, y in YS]
@time Animattion(all_pos, nt, field, "Videos/" * savew_dir * ".mp4"; framerate=60, isave=1, L=L)


"""
angle analyise
"""

function MeanOrientation(particle_orientations, particles_pos, L, ;num_layers)
    #* Define the size of your simulation box
    L

    #* Define the number of layers
    num_layers = num_layers

    #* Calculate the height of each layer
    layer_height = L / num_layers

    #* Initialize an array to store mean orientations for each layer
    # mean_orientations = zeros(num_layers)
    num_paricles_in_layer = zeros(num_layers)
    order_parameter = zeros(num_layers)
    # mean_vec = [SA[0.0 .0] for v in 1:length(particles_pos)]

    #* Iterate through each layer
    for i in 1:num_layers
        # @show i
        #* Calculate the boundaries of the current layer
        layer_bottom = (i - 1) * layer_height
        layer_top = i * layer_height

        #* Find particles within the current layer
        particles_in_layer = [(layer_bottom .<= v[2]) .& (v[2] .<= layer_top) for v in particles_pos]
        # particles_in_layer = (layer_bottom .<= particle_orientations) .& (particle_orientations .<= layer_top)

        # Calculate the mean orientation for particles in the current layer
        # mean_vec = mean(particle_orientations[particles_in_layer])
        order_parameter[i] = mean(cos.(particle_orientations[particles_in_layer]))
        num_paricles_in_layer[i] = length(particle_orientations[particles_in_layer])

        # mean_orientations[i] = atand(mean_vec[2], mean_vec[1])
    end

    return order_parameter, num_paricles_in_layer
end


# mean_angle, num_in_layer = MeanOrientation(all_ϕ[5000], all_pos[5000],L, ;num_layers=30)

# fig = Figure()
# ax1 = Axis(fig[1, 1], yticklabelcolor=:blue, ylabel="polar parameter", xlabel="y-axis / L", ylabelcolor="blue", 
#     xgridvisible=false, ygridvisible=false,
#     xlabelsize=28, ylabelsize=28, yticklabelsize=22, xticklabelsize=22)
# ax2 = Axis(fig[1, 1], yticklabelcolor=:red, yaxisposition=:right, ylabel="portation of MTs (%)", 
#     ylabelcolor="red", xgridvisible=false, ygridvisible=false,
#     xlabelsize=28, ylabelsize=28, yticklabelsize=22, xticklabelsize=22)
# hidespines!(ax2)
# hidexdecorations!(ax2, grid=false)
# # hidexdecorations!(ax1, grid=false)

# lines!(ax1, 0..1, mean_angle, color=:blue)
# lines!(ax2,  100*num_in_layer/1000, color=:red)
# fig
# lines(mean_angle)
