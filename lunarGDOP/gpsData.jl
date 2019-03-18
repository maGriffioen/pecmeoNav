function openOrbitData(filename::String)
    x = Array{Float64, 2}(undef, 0, 32)
    y = Array{Float64, 2}(undef, 0, 32)
    z = Array{Float64, 2}(undef, 0, 32)

    display("Warning: Loading ECEF data without conversion to non-rotating reference")
    open(filename) do file
        i_time = 0
        for ln in eachline(file)

            if ln[1:2] == "* "
                i_time += 1
                x = vcat(x, zeros(1, 32))
                y = vcat(y, zeros(1, 32))
                z = vcat(z, zeros(1, 32))

            elseif ln[1:2] == "PG"
                prn = parse(Int, ln[3:4])
                x[i_time , prn] = parse(Float64, ln[6:18])*1e3
                y[i_time , prn]= parse(Float64, ln[20:32])*1e3
                z[i_time , prn] = parse(Float64, ln[34:46])*1e3
            end
        end
    end

    return (x=x, y=y, z=z)
end
