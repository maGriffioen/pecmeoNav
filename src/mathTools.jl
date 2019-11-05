# Calculate Rodrigues rotation of 3d vectors
function vector_rotation(orgvec, refvec, angle)
    if length(orgvec) > 3 || length(refvec) > 3
        error("Rodrigues rotation needs 3d vectors")
    end
    k = normalize(refvec)
    v = orgvec
    v_new = v*cos(angle) + cross(k, v) * sin(angle) + k*dot(k, v) * (1-cos(angle))

    return v_new
end

# Calculate the moving average over a data series.
# Input:
# Data  Input data vector
# n::int
function movingAverage(data; n=3)
    n_max = length(data)
    global n_diff = ((n-1)÷2)
    if (n%2 ==1 && n_max >= n)
        new_data = zeros(length(data))

        for i in 1:n_max
            if ( n/2 > i)
                new_data[i] = mean([data[q] for q in 1:i+n_diff])
            elseif (i+n_diff <= n_max)
                new_data[i] = mean([data[q] for q in i-n_diff:i+n_diff])
            else
                new_data[i] = mean([data[q] for q in i-n_diff:n_max])
            end
        end


    elseif (n_max>=n)
        error("No even moving average implemented, or ")
    else
        error("Too little data")
    end

    return new_data
end

function sampleStandardDeviation(data)
    mn = mean(data)
    samples = length(data)
    ssd = sqrt((1/(samples-1)) * sum((data.-mn).^2))
    return ssd
end

# Calculate the mean of a data series
mean(data) = sum(data) / length(data)
rms(data) = sqrt(mean(data.^2))
rss(data) = sqrt(sum(data.^2))
ssd(data) = sampleStandardDeviation(data)

lin2log(linValue::Number) = 10 * log10(linValue)    #Linear to (10)logarithmic scale
log2lin(logValue::Number) = 10^(logValue / 10)      #(10)logarithmic to linear scale
