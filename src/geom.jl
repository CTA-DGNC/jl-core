module geom

    using StaticArrays
    using Printf

    export radian, degree, scalar,
           vector, 
           matrix_3x3, dcm, rot_mtx,
           quaternion, quaternion_raw, roll, pitch, yaw, norm, raw, C_ω,
           euler_angles   

    const scalar = Float64

    # https://en.wikibooks.org/wiki/Introducing_Julia/print
    # https://stackoverflow.com/questions/47452353/implementing-custom-primitive-types-in-julia
    # https://blog.glcs.io/julia-type-system
    # https://blog.glcs.io/julia-types-advanced

    include("geom/angles.jl") 

    # Vector en R³ 
    # https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/
    struct vector <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end

    # Vector en R² 
    module r2
        export vector, norm, matrix_2x2, dcm, rot_mtx 
        using StaticArrays
        using ..geom: radian, scalar

        struct vector <: FieldVector{2, scalar}
            x::scalar
            y::scalar
        end
        norm(v::vector)::scalar = sqrt(v.x^2+v.y^2) 
        const matrix_2x2 = SMatrix{2,2,scalar}        
        struct dcm <: FieldMatrix{2, 2, scalar}
            xx::scalar
            yx::scalar
            xy::scalar
            yy::scalar
        end  
        function rot_mtx(θ::radian)::dcm 
            s = sin(θ)
            c = cos(θ)
            dcm([c  s 
                -s  c])
        end              
    end

    # Matriz de rotación
    const matrix_3x3 = SMatrix{3,3,scalar}
    struct dcm <: FieldMatrix{3, 3, scalar}
        xx::scalar
        yx::scalar
        zx::scalar
        xy::scalar
        yy::scalar
        zy::scalar
        xz::scalar
        yz::scalar
        zz::scalar    
    end

    # Cuaternión
    const quaternion_raw = SVector{4,scalar}

    struct quaternion 
        re::scalar
        im::vector
    end

    quaternion(q::quaternion_raw)::quaternion = quaternion(q[1], q[2:4])

    #=
        Actitud en términos de ángulos de Euler con métodos de conversión
        Se adopta la secuencia intrínseca de Tait-Bryan: ψ:z, θ:y', ϕ:x"
    =#
    struct euler_angles
        ϕ::radian  # roll 
        θ::radian  # pitch
        ψ::radian  # yaw
    end

    # This allows the REPL to show values of type degree.
    Base.show(io::IO, e::euler_angles) = @printf("ϕ=%0.3f⁰, θ=%0.3f⁰, ψ=%0.3f⁰", float(e.ϕ)*r2d, float(e.θ)*r2d, float(e.ψ)*r2d)

    struct euler_vector <: FieldVector{3, scalar}
        x::scalar
        y::scalar
        z::scalar
    end

    include("geom/rot_mtx.jl")    
    include("geom/quaternion.jl")
    include("geom/euler_angles.jl")

end # module geom





