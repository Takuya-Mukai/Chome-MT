circle = Makie.Polygon([Point2f(cos(a), sin(a)) for a in range(0, 2π, length=30)])

# function Animattion(pos, dir, out_filepath::AbstractString = "test.mp4";
#     framerate::Integer=30, isave=1, L = 1
# )
 
#     # fig, ax = wireframe(Sphere(Point3f(0), scale), linewidth=0.5, color=(:grey, 0.2))
#     fig = Figure()
#     ax = Axis(fig[1, 1])
#     ax.aspect = DataAspect()
#     limit = (-L, L, -L, L, -L, L)
#     #ax = Axis3(fig[1, 1], limits=limit, aspect=:data)
#     # pos = Observable(Point2f.(first(rt)))
#     xs, ys = vec2comp.(pos)
#     x = Observable(vec2comp(pos[1])[1])
#     y = Observable(vec2comp(pos[1])[2])
#     u = Observable(vec2comp(dir[1])[1])
#     v = Observable(vec2comp(dir[1])[2])
#     # pos = Observable(first(rt))
#     # dir = Observable(first(nt))

#     # strength = (sqrt.(u .^ 2 .+ v .^ 2))

#     arrows!(ax, x, y, u, v, 
#         linewidth=1, lengthscale= L / 40,
#         align=:center)

#     # scatter!(ax, pos)

#     CairoMakie.record(fig, out_filepath, eachindex(dir); framerate=framerate) do frame_i
#         # pos[] = rt[frame_i]
#         # dir[] = nt[frame_i]
#         if frame_i % isave == 0
#             x[], y[] = vec2comp(pos[frame_i])
#             #  = ys[frame_i]
#             # u_tmp, v_tmp = vec2comp(dir[frame_i])
#             u[], v[] = vec2comp(dir[frame_i])
#             # strength[] = Observable(sqrt.(u_tmp .^ 2 .+ v_tmp .^ 2))
#         end
#         # v[] = vs[frame_i]
#     end
# end

function Animattion(pos, dir, out_filepath::AbstractString="test.mp4";
    framerate::Integer=30, isave=1, L=1)


    # fig, ax = wireframe(Sphere(Point3f(0), scale), linewidth=0.5, color=(:grey, 0.2))
    fig = Figure()
    ax = Axis(fig[1, 1])
    ax.aspect = DataAspect()
    limit = (-L, L, -L, L, -L, L)
    XS = 0:L
    YS = 0:L 
    #ax = Axis3(fig[1, 1], limits=limit, aspect=:data)
    # pos = Observable(Point2f.(first(rt)))
    # xs, ys = vec2comp.(pos)
    x = Observable(vec2comp(pos[1])[1])
    y = Observable(vec2comp(pos[1])[2])
    if dir[1] isa Vector{Float64}
        dir = [angle2dir.(ϕs) for ϕs in dir]
    end
    u = Observable(vec2comp(dir[1])[1])
    v = Observable(vec2comp(dir[1])[2])

    # pos = Observable(first(rt))
    # dir = Observable(first(nt))

    # strength = (sqrt.(u .^ 2 .+ v .^ 2))
    # C = Observable(field[1])
    # heatmap!(ax, XS, YS, C)
    arrows!(ax, x, y, u, v,
        linewidth=1, lengthscale=L / 40,
        align=:center)

    # scatter!(ax, pos)

    CairoMakie.record(fig, out_filepath, eachindex(dir); framerate=framerate) do frame_i
        # pos[] = rt[frame_i]
        # dir[] = nt[frame_i]
        if frame_i % isave == 0
            x[], y[] = vec2comp(pos[frame_i])
            #  = ys[frame_i]
            # u_tmp, v_tmp = vec2comp(dir[frame_i])
            u[], v[] = vec2comp(dir[frame_i])
            # strength[] = Observable(sqrt.(u_tmp .^ 2 .+ v_tmp .^ 2))
            # C[] = field[frame_i]
        end
        # v[] = vs[frame_i]
    end
end


function Animattion(pos, dir, spatio_func, out_filepath::AbstractString="test.mp4";
    framerate::Integer=30, isave=1, L=1, backend=GLMakie)

    # fig, ax = wireframe(Sphere(Point3f(0), scale), linewidth=0.5, color=(:grey, 0.2))
    fig = backend.Figure()
    limit = (-L, L, -L, L)
    ax = backend.Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)
    ax.aspect = backend.DataAspect()
    # limits!(ax, limit)
    ylims!(ax,0,L)
    xlims!(ax, 0,L)

    XS = 0:L
    YS = 0:L 

    field = [spatio_func.(x,y, L) for x in 0:0.2:L, y in 0:0.2:L]
    #ax = Axis3(fig[1, 1], limits=limit, aspect=:data)
    # pos = Observable(Point2f.(first(rt)))
    # xs, ys = vec2comp.(pos)
    x = backend.Observable(vec2comp(pos[1])[1])
    y = backend.Observable(vec2comp(pos[1])[2])
    if dir[1] isa Vector{Float64}
        dir = [angle2dir.(ϕs) for ϕs in dir]
    end
    u = backend.Observable(vec2comp(dir[1])[1])
    v = backend.Observable(vec2comp(dir[1])[2])

    # pos = Observable(first(rt))
    # dir = Observable(first(nt))

    # strength = (sqrt.(u .^ 2 .+ v .^ 2))
    C = backend.Observable(field)
    backend.heatmap!(ax, XS, YS, C)
    backend.arrows!(ax, x, y, u, v,
        linewidth=1, lengthscale=L / 40,
        align=:center)

    # scatter!(ax, pos)

    backend.record(fig, out_filepath, eachindex(dir); framerate=framerate) do frame_i
        # pos[] = rt[frame_i]
        # dir[] = nt[frame_i]
        if frame_i % isave == 0
            x[], y[] = vec2comp(pos[frame_i])
            #  = ys[frame_i]
            # u_tmp, v_tmp = vec2comp(dir[frame_i])
            u[], v[] = vec2comp(dir[frame_i])
            # strength[] = Observable(sqrt.(u_tmp .^ 2 .+ v_tmp .^ 2))
            # C[] = field[frame_i]
        end
        # v[] = vs[frame_i]
    end
end


function snapshort(pos, dir, spatio_func; backend=CairoMakie, colormap=:viridis,L=30)
    fig = backend.Figure()
    ax = backend.Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)
    ax.aspect = backend.DataAspect()
    ylims!(ax, 0, L)
    xlims!(ax, 0, L)
    XS = 0:0.2:L
    YS = 0:0.2:L
    field = [spatio_func.(x, y, L) for x in XS, y in YS]
    x, y = vec2comp(pos)
    if dir isa Vector{Float64}
        dir = [angle2dir.(ϕs) for ϕs in dir]
    end
    u,v= vec2comp(dir)
    # v = vec2comp(dir)

    # backend.heatmap!(ax, XS, YS, field, colormap=colormap, alpha=.1)
    upb = 2L/3
    lowb = L/2
    backend.poly!(ax, Point2f[(0, lowb), (L, lowb), (L, upb), (0, upb)], color=(colormap, 0.5))
    backend.arrows!(ax, x, y, u, v,
        linewidth=1, lengthscale=L / 40,
        align=:center)

    # hlines!([(2L / 3 - L / 2) / 2 + L / 2], linewidth=30, color=:blue)
    return fig
end
