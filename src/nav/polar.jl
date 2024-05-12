module polar

    using StaticArrays
    using ..geom: radian, degree, scalar 
    using ..geom.r2 #: vector, dcm, rot_mtx 

    export pos, vel
    
    struct pos
        θ::radian
        r::scalar 
    end
    struct vel
        γ::radian # elevación
        v::scalar # velocidad  
    end

        pos(p::vector)::pos  = pos(radian(atan(p.y,p.x)), norm(p))
     vector(p::pos )::vector = vector(p.r*cos(p.θ), p.r*sin(p.θ))
    horizon(p::pos )::radian = p.θ + radian(π/2)

    # ECEF a NED / NED a ECEF    
    Base.:<<(p::pos, q::radian)::radian = horizon(p) - q
    # NED a ECEF / ECEF a NED
    Base.:>>(p::pos, q::radian)::radian = horizon(p) - q     
    # Mapeo al horizonte local (módulo y ángulo respecto del horizonte)
    # (v, γ) = p << v
    function Base.:<<(p::pos, v::vector)::vel
        vp = pos(v)
        β  = horizon(p)
        γ  = β - vp.θ 
        if  scalar(γ) > 2*π 
            γ = γ - radian(2*π) 
        end
        vel(γ, vp.r)
    end  

end

