n = 10
a = zeros((n, n))

Threads.@threads for i in 1:n
    for j in 1:n
        c = i*j
        a[i, j] = c
    end
end

for i in 1:n
    println(a[i, :])
end