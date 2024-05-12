using LinearAlgebra

# Suma y resta de cuaterniones
function Base.:+(a::quaternion, b::quaternion)::quaternion
    quaternion( a.re + b.re, a.im + b.im)
end 
function Base.:-(a::quaternion, b::quaternion)::quaternion
    quaternion( a.re - b.re, a.im - b.im)
end 
function Base.:≈(a::quaternion, b::quaternion)::Bool
    if sign(a.re) == sign(b.re)
        return a.re ≈ b.re && a.im ≈ b.im
    end
    return a.re ≈ -(b.re) && a.im ≈ -(b.im)
end  
function Base.:conj(q::quaternion)::quaternion
    quaternion(q.re, -q.im)
end  
Base.:adjoint(q::quaternion)::quaternion = conj(q)

# Propiedades del cuaternión
Base.:real(q::quaternion)::scalar = q.re
Base.:imag(q::quaternion)::vector  = q.im
      nrm2(q::quaternion)::scalar = q.re*q.re + q.im' * q.im
      norm(q::quaternion)::scalar = sqrt(nrm2(q))

raw(q::quaternion)::quaternion_raw = quaternion_raw([q.re, q.im[1], q.im[2], q.im[3]])

# Producto de Hamilton
function Base.:*(a::quaternion, b::quaternion)::quaternion
    quaternion(  a.re   *b.re - a.im[1]*b.im[1] - a.im[2]*b.im[2] - a.im[3]*b.im[3], 
                    [ a.im[1]*b.re + a.re   *b.im[1] - a.im[3]*b.im[2] + a.im[2]*b.im[3]
                    , a.im[2]*b.re + a.im[3]*b.im[1] + a.re   *b.im[2] - a.im[1]*b.im[3]
                    , a.im[3]*b.re - a.im[2]*b.im[1] + a.im[1]*b.im[2] + a.re   *b.im[3]])
end 
×(a::quaternion, b::quaternion)::quaternion = a * b
⋅(a::quaternion, b::quaternion)::quaternion = a * b

# Conversiones a ángulos de Euler 
function euler_angles(q::quaternion)::euler_angles
    r = q.re;
    x = q.im[1];
    y = q.im[2];
    z = q.im[3];

    c31 =     2*(x.*z - r.*y);
    c11 = 1 - 2*(y.^2 + z.^2);
    c21 =     2*(x.*y + r.*z);
    c32 =     2*(y.*z + r.*x);
    c33 = 1 - 2*(x.^2 + y.^2);           

    euler_angles(atan(c32, c33), asin(-c31), atan(c21, c11));
end
function roll(q::quaternion)::radian
    r = q.re
    x = q.im[1]
    y = q.im[2]
    z = q.im[3]
    c32 =     2*(y.*z + r.*x)
    c33 = 1 - 2*(x.^2 + y.^2)          
    radian(atan(c32, c33))
end  
function pitch(q::quaternion)::radian
    r = q.re
    x = q.im[1]
    y = q.im[2]
    z = q.im[3]
    c31 = 2*(x.*z - r.*y)  
    radian(asin(-c31))
end 
function yaw(q::quaternion)::radian
    r = q.re
    x = q.im[1]
    y = q.im[2]
    z = q.im[3]
    c11 = 1 - 2*(y.^2 + z.^2)
    c21 =     2*(x.*y + r.*z)
    radian(atan(c21, c11))
end  

# Conversiones cuaternión / DCM
function quaternion(c::dcm)::quaternion
    t = c[1,1] + c[2,2] + c[3,3] 
    if t > -0.99999999
        n = 0.50 * sqrt(1 + t)
        f = 0.25 / n
        x = f*(c[3,2]-c[2,3])
        y = f*(c[1,3]-c[3,1]) 
        z = f*(c[2,1]-c[1,2])
    elseif (c[1,1] > c[2,2] && c[1,1] > c[3,3]) 
        p = sqrt(1.0 + c[1,1] - c[2,2] - c[3,3])
        s = 0.5 / p
        x = 0.5 * p
        y = (c[1,2] + c[2,1]) * s
        z = (c[3,1] + c[1,3]) * s
        n = (c[2,3] - c[3,2]) * s
    elseif (c[2,2] > c[3,3])
        p = sqrt(1.0 + c[2,2] - c[1,1] - c[3,3])
        s = 0.5 / p
        x = (c[1,2] + c[2,1]) * s
        y = 0.5 * p
        z = (c[2,3] + c[3,2]) * s
        n = (c[3,1] - c[1,3]) * s
    else 
        p = sqrt(1.0 + c[3,3] - c[1,1] - c[2,2])
        s = 0.5 / p
        x = (c[3,1] + c[1,3]) * s
        y = (c[2,3] + c[3,2]) * s
        z = 0.5 * p
        n = (c[1,2] - c[2,1]) * s
    end
    # lo conjugamos para que coincida con as_dcm = Cba
    quaternion(n, [-x, -y, -z])
end

# TRANSFORMACIONES

# dq/dt = C_ω(q)⋅ω
function C_ω(q::quaternion)::SMatrix{4,3,scalar}
    SMatrix{4,3,scalar}(
        [-q.im[1] -q.im[2] -q.im[3] 
          q.re    -q.im[3]  q.im[2] 
          q.im[3]  q.re    -q.im[1] 
         -q.im[2]  q.im[1]  q.re   ] * 0.5)
end

# proyección del vector 
function Base.:(<<)(q::quaternion, v::vector)::vector
    r = q.re
    x = q.im[1]; x2 = x*x;
    y = q.im[2]; y2 = y*y;
    z = q.im[3]; z2 = z*z;
    a = 2*((0.5-y2-z2) * v[1] + (x*y - r*z) * v[2] + (x*z + r*y) * v[3])
    b = 2*((x*y + r*z) * v[1] + (0.5-x2-z2) * v[2] + (y*z - r*x) * v[3])
    c = 2*((x*z - r*y) * v[1] + (y*z + r*x) * v[2] + (0.5-x2-y2) * v[3])

    vector(a, b, c)
end 
