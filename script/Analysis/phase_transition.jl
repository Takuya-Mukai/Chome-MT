include("../../src/Chemo-MT.jl")
using DelimitedFiles

"""
get all file with .jld2
"""
function get_file_list(dir)
    all_file = readdir(dir)
    file_list = []
    
    for i in eachindex(all_file)
        if split(all_file[i], ".")[end] == "jld2" && split(all_file[i], ".")[1] != ""
            append!(file_list, [all_file[i]])
        else
            nothing
        end
    end
    return file_list 
end


"""
Division into sample X
"""
function Division_into_smaple(file_list)

    #* find  number of sample
    sample_list = []
    for name in file_list
        sample_term = split(name, "_")[1]
        if length(sample_term) == 8
            Nsample = parse(Int, sample_term[end-1:end])
        else
            Nsample = parse(Int, sample_term[end])
        end
        
        if Nsample ∉ sample_list
            append!(sample_list, Nsample)
        end
    end
    Nsample = length(sample_list)
    sorted_file = [[] for _ in 1:Nsample]
    for (i, name) in  enumerate(file_list)
        sample_term = split(name, "_")[1]
        if length(sample_term) == 8
            Nsample = parse(Int, sample_term[end-1:end])
        else
            Nsample = parse(Int, sample_term[end])
        end
        append!(sorted_file[Nsample], [name])

    end
    return sorted_file
end



function phase_transistion(sorted_file, dir)
    Nsample = length(sorted_file)
    Neta = length(sorted_file[1])
    eta_list = zeros(Neta)
    op_list = []
    Threads.@threads for sample in sorted_file
        op = zeros(Neta)
        for i_eta in 1:Neta
            data = load(dir*sample[i_eta])
            all_ori = data["all_ori"]

            op[i_eta] = get_order_parameter(all_ori)
            eta_list[i_eta] = parse(Float64, split(sample[i_eta], "_")[end][1:3])
        end
        
        append!(op_list, [op])
    end

    return eta_list, mean(op_list)
end

function XYoder(sorted_file, dir)
    Neta = length(sorted_file[1])
    eta_list = zeros(Neta)
    xorder_list = []
    yorder_list = []
        op_list = []

    Threads.@threads for sample in sorted_file
        xop = zeros(Neta)
        yop = zeros(Neta)
        op = zeros(Neta)
        for i_eta in 1:Neta
            data = load(dir * sample[i_eta])
            all_ori = data["all_ori"]
            op[i_eta] = get_order_parameter(all_ori)
            xop[i_eta], yop[i_eta] = get_XYorder_parameter(all_ori)
            eta_list[i_eta] = parse(Float64, split(sample[i_eta], "_")[end][1:3])
        end
        append!(op_list, [op])
        append!(xorder_list, [xop])
        append!(yorder_list, [yop])
    end

    return mean(op_list), mean(xorder_list), mean(yorder_list), eta_list
end



function XYorder_inside(sorted_file, dir; xmin=0, xmax=30, ymin=0, ymax=30)
    Neta = length(sorted_file[1])
    eta_list = zeros(Neta)
    xorder_list = []
    yorder_list = []
    op_list = []
    Threads.@threads for sample in sorted_file
        xop = zeros(Neta)
        yop = zeros(Neta)
        op = zeros(Neta)
        for i_eta in 1:Neta
            data = load(dir * sample[i_eta])
            all_ori = data["all_ori"]
            op[i_eta] = get_order_parameter_inside(all_pos, all_ori, xmin, xmax, ymin, ymax)
            xop[i_eta], yop[i_eta] = get_XYorder_parameter_inside(all_pos, all_ori, xmin, xmax, ymin, ymax)
            eta_list[i_eta] = parse(Float64, split(sample[i_eta], "_")[end][1:3])
        end
        append!(op_list, [op])
        append!(xorder_list, [xop])
        append!(yorder_list, [yop])
    end

    return mean(op_list), mean(xorder_list), mean(yorder_list), eta_list
end

function output_order_parameters_to_txt(eta, op, xo, yo, js, dir = loaddir)
    open(dir*"order_parameters.txt", "w") do file
        writedlm(file, [eta, op, xo, yo, js])
    end
end
# loaddir = "/Volumes/Extreme SSD/Kyle/Simulation_Data/VS_sin_v1/"
# loaddir = "F:/Kyle/Simulation_Data/VS_step_funcJ1_v1/"
# loaddir = "//NAS-Kawamata/kyoto/Kyle Chan/Simulation_Data/Anisotropic_MTs_Swarm/Data/VS_sin1_v2/"
# loaddir = "NAS-Kawamata._smb._tcp.local/kyoto/Kyle Chan/Simulation_Data/Anisotropic_MTs_Swarm/Data/VS_sin05_v1/"
# loaddir = "data/VS_sin2pi_v1/"

# loaddir = "/home/kylechan/Dropbox/Research/Project/Chome-MT/Data/VS_RL_4_v1/"
# loaddir = "/media/kylechan/Extreme SSD/Kyle/Simulation_Data/VS_RL_4_L50_v2/"
# loaddir = "Data/post_data/yStep_Lfrac02_L30_v1/"
loaddir = "Data/Const_J1_L30_v1/"
# loaddir = "/Volumes/Extreme SSD/Kyle/Simulation_Data/VS_step_funcJ1_v1/"

# loaddir = "/mnt/SSD_RAID0/kyle/Viscek_anisotropic/yStep_v1/"

file_list = get_file_list(loaddir)
sorted_file = Division_into_smaple(file_list)

@time eta_list, op_list = phase_transistion(sorted_file, loaddir)
println("analysing...")
# @time op_list, xorder_list, yorder_list, eta_list = XYoder(sorted_file, loaddir)
# @time op_list, xorder_list, yorder_list, eta_list = XYorder_inside(sorted_file, loaddir; xmin=0, xmax=L, ymin=L/2, ymax=2L/3)



# S1 = sorted_file[5]
# data = load(loaddir*S1[2])
# all_ori = data["all_ori"]
# all_pos = data["all_pos"]


# L = 30
# function spatio_func(x,y, L)
#    (L/2 < y < 2L/3) ? 1.0 : 0.0
# end
# spatio_func(x,y,L) = 1.0

function spatio_func(x,y,L)
    L_frac = .2
    up = L/2 + L*L_frac/2
    low = L/2 - L*L_frac/2
    if y > low && y < up
        return 1.0
    else
        return 0.0
    end 
end

# Jmax = maximum(spatio_func.(0:L, 0:L, L)) 
# Jmin = minimum(spatio_func.(0:L, 0:L, L))
# J_field = [spatio_func(x,y,L) for x in 0:0.2:L, y in 0:0.2:L]
# Jmean = mean(J_field)

# output_order_parameters_to_txt(eta_list, op_list, xorder_list, yorder_list, [Jmean for i in 1:length(op_list)])

# @time Animattion(all_pos, all_ori, spatio_func,join(split(loaddir * S1[3], ".")[1:end-1], ".") * ".mp4"; framerate=60, L=30)

# frame = 500
# fig = snapshort(all_pos[frame], all_ori[frame], spatio_func; L=30, colormap=:blue)
# display(fig)
# save(loaddir * "snapshort.pdf", fig)

# fig = Figure()
# ax = Axis(fig[1, 1], xlabel=L"\eta/\pi", ylabel="order parameter",
#     xgridvisible=false, ygridvisible=false)
# scatterlines!(eta_list ./ π, op_list, color=Makie.wong_colors()[3], label=L"\left< v \right>")
# scatter!(eta_list ./ π, xorder_list, label=L"\left< v_x \right>")
# scatter!(eta_list ./ π, yorder_list, label=L"\left< v_y \right>")

# vlines!([Jmax], linestyle=:dash, color="black", label=L"\text{Max}(J)")
# vlines!([Jmean], linestyle=:dot, color="black", label=L"\left< J \right>")

# axislegend(framevisible=false)
# display(fig)
# save(loaddir*"phase_transition_in_XY.png",fig)
