using LinearAlgebra
using Richardson: extrapolate
using DelimitedFiles
#
module constant
    Ω = 2^2 * 5.368153453*1e-5 # eV*Å^3
#    !!!!!!!!!!!!!!!!!!
    kB = 8.617333262e-5 # eV/K
    μB = 5.7883818060e-5 # eV/T
end
#
module latt  #CrSBr
    a = 3.54780687
    b = 4.740272313
    𝛿x = a/2
    𝛿y = b/2 
    𝛿z = 2.065678385
end
#
# DipoleFT is a function that calculates FT of the dipole energy per atom for a cluster (n * n)
# k is the 2-dim k-vector
# M is the vector pointing in the magnetization direction
function DipoleFT(n, M, k)
    energy1 = 0+0im
    energy2 = 0+0im
    k_scale = [ k[1]/latt.a, k[2]/latt.b ] # input k is in units of 1/a (1/b)
    K = [k_scale; 0.0]
    for ix = -n:n, iy = -n:n
        if (ix==0)&&(iy==0)
            R2 = [ latt.𝛿x, latt.𝛿y, latt.𝛿z]
            KdotR2 = K⋅R2
            fac2 = 2π*KdotR2*1im
            cosij2 = (M ⋅ R2) / (norm(R2)*norm(M))
            energy2 += (1 - 3*cosij2^2) / norm(R2)^3 *exp(fac2)
            continue
        end
        R1 = [ latt.a*ix, latt.b*iy, 0.0]
        KdotR1 = K⋅R1
        fac1 = 2π*KdotR1*1im
        R2 = [ latt.a*ix + latt.𝛿x, latt.b*iy + latt.𝛿y, latt.𝛿z]
        KdotR2 = K⋅R2
        fac2 = 2π*KdotR2*1im
        if (ix==0)&&(iy==0)
            continue
        end
        cosij1 = (M ⋅ R1) / (norm(R1)*norm(M))
        cosij2 = (M ⋅ R2) / (norm(R2)*norm(M))
        energy1 += (1 - 3*cosij1^2) / norm(R1)^3 *exp(fac1)
        energy2 += (1 - 3*cosij2^2) / norm(R2)^3 *exp(fac2)
    end
    return real(energy1) + real(energy2)
    end
#        
function pxx(k::Vector{<:Real})
        if size(k)[1] != 2
            error("Dimension of the pxx argument k is different from 2")
        end
        return extrapolate(N -> DipoleFT(Int(N),[1 0 0],k), 8, contract=0.5, atol=1e-8, x0=Inf, maxeval=4)[1]  
end
# 
function pyy(k::Vector{<:Real})
        if size(k)[1] != 2
            error("Dimension of the pyy argument k is different from 2")
        end
        return extrapolate(N -> DipoleFT(Int(N),[0 1 0],k), 8, contract=0.5, atol=1e-8, x0=Inf, maxeval=4)[1]  
end
#   
function pzz(k::Vector{<:Real})
        if size(k)[1] != 2
            error("Dimension of the pzz argument k is different from 2")
        end
        return extrapolate(N -> DipoleFT(Int(N),[0 0 1],k), 8, contract=0.5, atol=1e-8, x0=Inf, maxeval=4)[1]  
end
#
#
    function write_pxxyyzz(filename,dk)
    file = open(filename,"w")
    for kx = -0.5:dk:0.5, ky = -0.5:dk:0.5
    writedlm(file, [ kx ky pxx([kx,ky]) pyy([kx,ky]) pzz([kx,ky]) ])
#    @printf("%14.10f %12.10f %16.12f %16.12f %16.12f \n", kx, ky, pxx([kx,ky]), pzz([kx,ky]), pzz([kx,ky]))
#    println(i," ",j," ",pxx([i/100, j/100])*S^2*constant.Ω/2, " ",pyy([i/100, j/100])*S^2*constant.Ω/2, " ", pzz([i/100, j/100])*S^2*constant.Ω/2)
    end
    close(file)
    end 
#
#
function read_pxxyyzz(filename,dk)
        file = open(filename,"r")
        dlm = readdlm(file,Float64)
        close(file)
        nk = floor(Int,1/dk) + 1
        pxx = transpose(reshape(dlm[:,3],nk,nk))
        pyy = transpose(reshape(dlm[:,4],nk,nk))
        pzz = transpose(reshape(dlm[:,5],nk,nk))
        return [pxx,pyy,pzz]
end   
#
#
