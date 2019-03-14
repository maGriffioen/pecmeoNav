function findNavCovmat(
    refPoint::Tuple{Number, Number, Number},
    constel::KeplerConstellation )

    n_navsats = size(constel)       #Number of navigation satellites
    mat_A = zeros(n_navsats, 4)     #Create empty A matrix
    mat_A[:, 4] .= -1.0;            #Set -1 for time positions

    #Loop over the navigation satellites in the constellation
    for isat in 1:n_navsats
        posnavsat = position(constel[isat])         #Nav satellite position
        posdiff = collect(posnavsat .- refPoint)    #Position difference
        r  = norm(posdiff)             #Total distance
        mat_A[isat, 1:3] = posdiff ./ r             #Write unit vector distances
    end

    return inv(mat_A' * mat_A)
end

function findNavGDOP(
    refPoint::Tuple{Number, Number, Number},
    constel::KeplerConstellation )

    mat_Q = findNavCovmat(refPoint, constel)
    return sum(diag(mat_Q))     #GDOP
end
