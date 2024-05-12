#= 
    axis: 1/2/3 o 'x'/'y'/'z'
    angle: [rad]
    S proyecta vectores de la terna inicial en la final 
=#
function rot_mtx(axis, angle::radian)::dcm
    c = cos(scalar(angle));
    s = sin(scalar(angle));
    if     axis == 1 || axis =='x' 
        return dcm([1 0 0 ; 0 c s ; 0 -s c])
    elseif axis == 2 || axis =='y' 
        return dcm([c 0 -s ; 0 1 0 ; s 0 c])
    elseif axis == 3 || axis =='z' 
        return dcm([c s 0 ; -s c 0 ; 0 0 1])
    end
end  

rot_mtx(axis, angle::degree) = rot_mtx(axis, radian(angle))

function rot_mtx(q::quaternion)::dcm
    r = q.re
    x = q.im[1]
    y = q.im[2]
    z = q.im[3]

    dcm([1-2*(y^2 + z^2)   2*(x*y + r*z)   2*(x*z - r*y) 
           2*(x*y - r*z) 1-2*(x^2 + z^2)   2*(y*z + r*x) 
           2*(x*z + r*y)   2*(y*z - r*x) 1-2*(x^2 + y^2)])
end 

function rot_mtx(e::euler_angles)::dcm
    sϕ = sin(scalar(e.ϕ))
    sθ = sin(scalar(e.θ))
    sψ = sin(scalar(e.ψ))

    cϕ = cos(scalar(e.ϕ))
    cθ = cos(scalar(e.θ))
    cψ = cos(scalar(e.ψ))
    
    dcm([cθ*cψ           cθ*sψ          -sθ 
         sϕ*sθ*cψ-cϕ*sψ  sϕ*sθ*sψ+cϕ*cψ  sϕ*cθ 
         cϕ*sθ*cψ+sϕ*sψ  cϕ*sθ*sψ-sϕ*cψ  cϕ*cθ])
end 


