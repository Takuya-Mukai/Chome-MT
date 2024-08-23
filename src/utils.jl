"""
    obj2txt(object, file_name)
Write a julia object to txt file with field name and value
"""
function obj2txt(object, file_name)
    open(file_name, "w") do io
        T = typeof(object)
        for name in fieldnames(T)
            write(io, "$name = $(getfield(object, name))\n")
        end
    end
end

println("loaded")