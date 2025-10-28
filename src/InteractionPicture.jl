module InteractionPicture

struct IPfy{lT, tT, N}
    linop::lT
    transform::tT
    L::Array{ComplexF64, N}
    y::Array{ComplexF64, N}
    f::Array{ComplexF64, N}
end

function IPfy(linop::AbstractArray, transform)
    IPfy(linop, transform, similar(linop), similar(linop), similar(linop))
end

function update_L!(L, linop::AbstractArray, t)
    for ii in eachindex(L)
        L[ii] = linop[ii]*t
    end
end

function (ip::IPfy)(fbar, ybar, _, t)
    update_L!(ip.L, ip.linop, t)
    @. ip.y = ybar * exp(ip.L)
    ip.transform(ip.f, ip.y, t)
    @. fbar = exp(-ip.L) * ip.f 
    nothing
end

end