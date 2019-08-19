# Calculate the moving average over a data series.
# Input:
# Data  Input data vector
# n::int
function movingAverage(data; n=3<:Int)
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

# Calculate the mean of a data series
mean(data) = sum(data) / length(data)

lin2log(linValue::Number) = 10 * log10(linValue)    #Linear to (10)logarithmic scale
log2lin(logValue::Number) = 10^(logValue / 10)      #(10)logarithmic to linear scale
