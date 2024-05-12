#include("geom.jl")

module nav

    using ..geom
    using StaticArrays
    using Printf
    using LinearAlgebra

    export geo_pos , nav_pos  , 
           ecef_pos, ecef_vel , 
           eci_pos , eci_vel  ,  
           ned_pos , ned_vel  ,  
           enu_pos , enu_vel  ,
           WGS84 , ECI  

    const scalar = geom.scalar

    # posición sobre la superficie en coordenadas de navegación
    struct geo_pos
        λ::radian
        φ::radian
    end
    # posición en coordenadas de navegación (LLA)
    struct nav_pos
        λ::radian
        φ::radian
        h::scalar 
    end
    # velocidad en terna LGV
    struct lgv_vel
        γ::radian  # elevación
        Ψ::radian  # rumbo
        v::scalar  # velocidad
    end

    # posición y velocidad terna Earth Centered Earth Fixed
    struct ecef_pos <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end
    struct ecef_vel <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end    

    # posición y velocidad terna Earth Centered Inertial
    struct eci_pos <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end
    struct eci_vel <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end  

    # posición y velocidad terna local East North Up
    struct enu_pos <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end
    struct enu_vel <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end  

    # posición y velocidad terna local North East Down
    struct ned_pos <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end
    struct ned_vel <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end  

    geo_pos(p::nav_pos)::geo_pos = geo_pos(p.λ, p.φ)
    nav_pos(p::geo_pos)::nav_pos = nav_pos(p.λ, p.φ, 0)
    ned_pos(p::enu_pos)::ned_pos = ned_pos(p.y, p.x, -p.z)
    ned_vel(v::enu_vel)::ned_vel = ned_vel(v.y, v.x, -v.z)
    enu_pos(p::enu_pos)::enu_pos = enu_pos(p.y, p.x, -p.z)
    enu_vel(v::enu_vel)::enu_vel = enu_vel(v.y, v.x, -v.z)

    include("nav/WGS84.jl")
    include("nav/sphere.jl")
    include("nav/polar.jl")
    include("nav/ECI.jl")

     nav_pos(p::ecef_pos)::nav_pos = WGS84.to_lgv(p)
    ecef_pos(p::nav_pos)::ecef_pos = WGS84.to_ecef(p)

end

