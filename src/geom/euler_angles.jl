function Base.:+(a::euler_angles, b::euler_angles)::euler_angles
    euler_angles(a.ϕ + b.ϕ, a.θ + b.θ, a.ψ + b.ψ)
end 
function Base.:-(a::euler_angles, b::euler_angles)::euler_angles
    euler_angles(a.ϕ - b.ϕ, a.θ - b.θ, a.ψ - b.ψ)
end  
function Base.:≈(a::euler_angles, b::euler_angles)::Bool
   # (abs(a.ϕ - b.ϕ) < 1e-15 && abs(a.θ - b.θ) < 1e-15 && abs(a.ψ - b.ψ) < 1e-15)
    a.ϕ ≈ b.ϕ && a.θ ≈ b.θ && a.ψ ≈ b.ψ 
end  

function quaternion(e::euler_angles)::quaternion
    sϕ = sin(scalar(e.ϕ)/2)
    sθ = sin(scalar(e.θ)/2)
    sψ = sin(scalar(e.ψ)/2)

    cϕ = cos(scalar(e.ϕ)/2)
    cθ = cos(scalar(e.θ)/2)
    cψ = cos(scalar(e.ψ)/2)

    r = cϕ*cθ*cψ + sϕ*sθ*sψ
    x = sϕ*cθ*cψ - cϕ*sθ*sψ 
    y = cϕ*sθ*cψ + sϕ*cθ*sψ
    z = cϕ*cθ*sψ - sϕ*sθ*cψ

    quaternion(r, [x, y, z]) 
end
