using Pkg ; Pkg.activate(".")
using JLD2
using StaticArrays
using Statistics
using LinearAlgebra
using GLMakie

file_list = [
    "Data/horizontal_N3000_nstep10000_sample1.jld2",
    "Data/horizontal_N3000_nstep10000_sample2.jld2",
    "Data/horizontal_N3000_nstep10000_sample3.jld2",
    "Data/horizontal_N3000_nstep10000_sample4.jld2",
    "Data/horizontal_N3000_nstep10000_sample5.jld2",
    "Data/horizontal_N3000_nstep10000_sample6.jld2",
    "Data/horizontal_N3000_nstep10000_sample7.jld2",
    "Data/horizontal_N3000_nstep10000_sample8.jld2",
    "Data/horizontal_N3000_nstep10000_sample9.jld2",
    "Data/horizontal_N3000_nstep10000_sample10.jld2",
    "Data/horizontal_N3000_nstep10000_sample11.jld2",
    "Data/horizontal_N3000_nstep10000_sample12.jld2",
    "Data/horizontal_N3000_nstep10000_sample13.jld2",
    "Data/horizontal_N3000_nstep10000_sample14.jld2",
    "Data/horizontal_N3000_nstep10000_sample15.jld2",
    "Data/horizontal_N3000_nstep10000_sample16.jld2",
    "Data/horizontal_N3000_nstep10000_sample17.jld2",
    "Data/horizontal_N3000_nstep10000_sample18.jld2",
    "Data/horizontal_N3000_nstep10000_sample19.jld2",
    "Data/horizontal_N3000_nstep10000_sample20.jld2",
    "Data/horizontal_N3000_nstep10000_sample21.jld2",
    "Data/horizontal_N3000_nstep10000_sample22.jld2",
    "Data/horizontal_N3000_nstep10000_sample23.jld2",
    "Data/horizontal_N3000_nstep10000_sample24.jld2",
    "Data/horizontal_N3000_nstep10000_sample25.jld2",
    "Data/horizontal_N3000_nstep10000_sample26.jld2",
    "Data/horizontal_N3000_nstep10000_sample27.jld2",
    "Data/horizontal_N3000_nstep10000_sample28.jld2",
    "Data/horizontal_N3000_nstep10000_sample29.jld2",
    "Data/horizontal_N3000_nstep10000_sample30.jld2",
    "Data/horizontal_N3000_nstep10000_sample31.jld2",
    "Data/horizontal_N3000_nstep10000_sample32.jld2",
    "Data/horizontal_N3000_nstep10000_sample33.jld2",
    "Data/horizontal_N3000_nstep10000_sample34.jld2",
    "Data/horizontal_N3000_nstep10000_sample35.jld2",
    "Data/horizontal_N3000_nstep10000_sample36.jld2",
    "Data/horizontal_N3000_nstep10000_sample37.jld2",
    "Data/horizontal_N3000_nstep10000_sample38.jld2",
    "Data/horizontal_N3000_nstep10000_sample39.jld2",
    "Data/horizontal_N3000_nstep10000_sample40.jld2",
]

angle2dir(ϕ) = SA[cos(ϕ), sin(ϕ)]

function vec2comp(vec)
    v1 = [v[1] for v in vec]
    v2 = [v[2] for v in vec]

    return v1, v2
end


function plot_mean_dir(file_list)
    fig = Figure()
    ax = Axis(fig[1,1], yticks=[-180,-90,0,90,180])
    for dir in file_list
        data = load(dir)
        # rs = data["pos"]
        ϕs = data["angle"][1:1:end]
        # N = length(ϕs)
        nt = [angle2dir.(ϕ) for ϕ in ϕs]
        nt_bar = mean.(nt)
        ϕ_bar = [atand(v[2], v[1]) for v in nt_bar]
        # ϕ_deg = [rad2deg.(phi) for phi in ϕs]

        # ts = 1:2e-4:(N-1)*2e-4
        scatter!(ax, ϕ_bar[1:1:end])
    end
    # (ax, [-180, 180])
    hlines!(ax, [90, -90])
    # yticks!(ax, -180:90:180)
    display(fig)
end


function plot_mean_dir_auto_read_dir(Nsample)
    fig = Figure()
    ax = Axis(fig[1, 1], yticks=[0, 90, 180], xlabel="step", ylabel="<θ>")
    # dir_base = "Data/horizontal_N300_nstep10000_sample"
    dir_base = "Data/horizontal_vel0210_N500_eps0.01_nstep60000_dt0.002_sample"
    # dir_base = "Data/horizontal__N500_eps1.0_nstep60000_dt0.002_sample"
    for i in 1:Nsample
        dir = dir_base * "$i.jld2"
        data = load(dir)
        # rs = data["pos"]
        ϕs = data["angle"][1:100:end]
        N = length(ϕs)
        ts = collect(1:N)
        nt = [angle2dir.(ϕ) for ϕ in ϕs]
        nt_bar = mean.(nt)
        ϕ_bar = [atand(v[2], v[1]) for v in nt_bar]
        # ϕ_deg = [rad2deg.(phi) for phi in ϕs]

        # ts = 1:2e-4:(N-1)*2e-4
        lines!(ax, ts, abs.(ϕ_bar[1:1:end]))

    end
    # (ax, [-180, 180])
    hlines!(ax, [90])
    # yticks!(ax, -180:90:180)
    display(fig)
    return fig
end

# plot_mean_dir(file_list)
fig = plot_mean_dir_auto_read_dir(20)
# fig
# dir = file_list[1]
# data = load(dir)
# pos = data["pos"]
# all_ϕ = data["angle"]
# nt = [angle2dir.(ϕs) for ϕs in all_ϕ]

# xs, ys = vec2comp.(pos)
# us, vs = vec2comp.(nt)
# savew_dir = join(split(split(dir, "/")[2], ".")[1:end-1], ".")
# @time Animattion(pos, nt, savew_dir*".mp4"; framerate=300, isave=60)