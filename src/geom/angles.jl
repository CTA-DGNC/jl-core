primitive type radian <: Real 64 end
primitive type degree <: Real 64 end

# A constructor to create values of the type radian
radian(x::scalar)    = reinterpret(radian, x)
radian(x::Int64)      = reinterpret(radian, convert(scalar,x))
radian(x::Irrational) = reinterpret(radian, convert(scalar,x))

# A constructor to convert back.
scalar(x::radian) = reinterpret(scalar, x)
    Int64(x::radian) = reinterpret(Int64, x)
    Real(x::radian) = reinterpret(Real, x)
Base.:AbstractFloat(x::radian) = reinterpret(scalar, x)

Base.:≈(a::radian, b::radian)::Bool   = scalar(a) ≈ scalar(b) 
Base.:+(a::radian, b::radian)::radian = radian(scalar(a)+scalar(b))
Base.:-(a::radian, b::radian)::radian = radian(scalar(a)-scalar(b))
Base.:*(a::radian, b::Real  )::radian = (scalar(a)*scalar(b))
Base.:/(a::radian, b::Real  )::radian = (scalar(a)/scalar(b))
Base.:≈(a::radian, b::Real  )::Bool   = scalar(a) == b
Base.:-(a::radian)::radian = radian(-scalar(a))
Base.:abs(a::radian)::scalar = abs(scalar(a)) 

 degree(x::scalar) = reinterpret(degree, x)
 degree(x::Int64  ) = reinterpret(degree, convert(scalar,x))
scalar(x::degree ) = reinterpret(scalar, x)
  Int64(x::degree ) = reinterpret(Int64, x)

Base.:≈(a::degree, b::degree)::Bool = scalar(a) ≈ scalar(b) 
Base.:+(a::degree, b::degree)::degree = degree(scalar(a)+scalar(b))
Base.:-(a::degree, b::degree)::degree = degree(scalar(a)-scalar(b))
Base.:-(a::degree)::degree = degree(-scalar(a))
Base.:abs(a::degree)::scalar = abs(scalar(a)) 

const global d2r = π/180
const global r2d = 180/π

radian(x::degree) = reinterpret(radian, reinterpret(scalar, x)*d2r)
degree(x::radian) = reinterpret(degree, reinterpret(scalar, x)*r2d)

# This allows the REPL to show values of type radian and degree.
Base.show(io::IO, x::radian) = print(io, scalar(x))
Base.show(io::IO, x::degree) = print(io, scalar(x))
