
module WGS84

    using ..geom
    using ..nav

    #=
    Info:
        https://icgem.gfz-potsdam.de/calcpoints
        https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84
    
    Elementos geodráficos para la navegación, usando el modelo WGS84 
    Nomenclatura:
    C_ba: (a -> b) DCM para proyección desde terna {a} en la terna {b}.
        Vb = Cba·Va 
    ECEF: (e) earth centered earth fixed
    LGV : local geodesic vertical 
    LGCV: local geocentric vertical 
    =#  

    const global e   = 0.08181919084
    const global e²  = 0.006694379990140            # f·(2 - f)        Square of (first) eccentricity
    const global ep² = e²/(1-e²)                    # Square of second eccentricity 
    const global f   = 1/298.257223563
    const global a   = 6378137                      # Radio ecuatorial
    const global a²  = a^2
    const global b   = a*√(1-e²)                    # Radio polar: 6356752.314245 
    const global M   = 5.9722e24                    # ±0.0006 
    const global G   = 6.674e-11
    const global GM  = 3.986004418000000e+014       # Constante gravitacional
    const global go  = 9.80665
    const global gp  = 9.832                        # gravitación en los polos
#   const global ge  = 9.7982876                    # gravitación en el ecuador - Atronautics, p.173
    const global ge  = 9.8142                       # gravitación en el ecuador - WGS84
    const global J2  = 1.082629821313305e-3
    const global Ωe  = 7.2921151467e-5              # rotación terrestre en rad/s
    const global Ωe² = Ωe^2

    const SWe = matrix_3x3([0  -Ωe  0
                            Ωe  0   0 
                            0   0   0])
    const we  = vector([0 , 0 , Ωe])
    const c   = 299792458 # velocidad de la luz en el vacío  
    # obliquity: inclinación del eje de rotación terrestre recpecto de su plano orbital  
    # T: Julian centuries from J2000.0
    ϵ(T::Int16)::degree = degree(23+26/60+(21.26 - 46.845*T- 0.0059*T²+ 0.00181*T³)/3600)
    const tlt = degree(23.44)      
    # DCM del eje terrestre
    const Ctl = geom.rot_mtx('z', tlt)
    # Rn: radio normal
    # λ: latitud geodésica
    Rn(λ::radian)::scalar = scalar(a/√(1-e²*sin(λ)^2))
    # ps: punto normal sobre el elipsoide en el plano meridiano (2D)
    # λ : latitud geodésica
    function Ps(λ::radian)::vector_r2 
        sr = √(1-e²*sin(λ)^2)
        vector_r2(a*cos(λ)/sr, a*(1-e²)*sin(λ)/sr)
    end
    # λ: latitud geocéntrica
    # r  : radio geocéntrico
    Rc(λ::radian)::scalar = scalar(a/√(1-sin(λ)^2*e²))

    function to_lgv(pe::ecef_pos)::nav_pos # [verificado]
        R2 = pe' * pe 
        if R2 < 1
            return nav_pos(0, 0, 0)
        end
        φ  = atan(pe[2], pe[1])
        Rp = √(pe[1]^2 + pe[2]^2)
        if Rp > 0.01
            # 	            lgc = asin(pe[3]/R)        # geocéntrica
            #                 h   = 0
            #                 R   = a
            #                 for n = 1:5
            #                     λ = atan(tan(lgc)*(R+h)/(R*(1-e²)+h))
            #                     h1  = h
            #                     h   = Rp/cos(λ) - a/√(1-e²*sin(λ)^2)
            #                     if abs(h - h1) < 0.001
            #                         break
            #                     end
            #                 end  
            # Spheroid properties
            z = pe[3]
            # Bowring's formula for initial parametric (β) and geodetic
            # (phi) latitudes
            β = atan(z, (1-f) * Rp)
            λ = atan(z + b*ep²*sin(β)^3, Rp - a*e²*cos(β)^3)
            # Fixed-point iteration with Bowring's formula
            # (typically converges within two or three iterations)
            βn = atan((1-f)*sin(λ), cos(λ))
            count = 0
            while β ≠ βn && count < 5
                β  = βn
                λ  = atan(z + b*ep²*sin(β)^3, Rp - a*e²*cos(β)^3)
                βn = atan((1-f)*sin(λ), cos(λ))
                count = count + 1
            end
            # Ellipsoidal height from final value for latitude
            slat = sin(λ)
            N = a/√(1 - e²*slat^2)
            h = Rp*cos(λ) + (z + e²*N*slat)*slat - N
        else
            λ = π/2*sign(pe[3])
            h = abs(pe[3]) - b 
        end
        nav_pos(λ,φ,h)
    end 

    function to_ecef(p::nav_pos)::ecef_pos # [verificado]
        rn  = Rn(p.λ); r = (rn+p.h)*cos(p.λ)
        ecef_pos(r*cos(p.φ), r*sin(p.φ), ((1-e²)*rn+p.h)*sin(p.λ))
    end 

    # q^g_e: ecef_to_lgv_enu_quat: rotación desde terna ECEF a ENU 
    function q_ge(p::geo_pos)::quaternion
        #     qa = quaternion(cos(λ/2), 0, -sin(λ/2), 0)
        #     ql = quaternion(cos(φ/2), [0 0 sin(φ/2)]) #sin(φ/2) * [sin(λ) 0 cos(λ)])
        #     q = qa * ql

        sl = sin(p.φ)
        cl = cos(p.φ)
        sq = sin(p.λ)
        cq = cos(p.λ)

        c11 =    -sl; c12 =     cl; c13 =  0 
        c21 = -sq*cl; c22 = -sq*sl; c23 = cq
        c31 =  cq*cl; c32 =  cq*sl; c33 = sq    
        r = 0.5*√(1+c11+c22+c33)
        quaternion(r, 0.25/r*[c32-c23 c13-c31 c21-c12])
    end
    # C^g_e ecef_to_lgv_enu_dcm
    function C_ge(p::geo_pos)::dcm
        sl = sin(p.φ)
        cl = cos(p.φ)
        sq = sin(p.λ)
        cq = cos(p.λ)

        dcm([-sl     cl    0 
             -sq*cl -sq*s cq
              cq*cl  cq*sl sq ])
    end   
    function gravitation(p::ecef_pos)::vector 
        r1 = norm(p)
        r2 = r1*r1
        r3 = r2*r1
        r4 = r3*r1
        m1 = GM/r3
        m2 = 3*J2*a²/(2*r2)
        m3 = 15*J2*a²*(p.z^2)/(2*r4)
        vector(-(p.x*m1)*(1+  m2-m3), 
               -(p.y*m1)*(1+  m2-m3),
               -(p.z*m1)*(1+3*m2-m3))
    end
    # gravitación + centrífuga
    function gravity(p::ecef_pos)::vector 
        g_e = gravitation(p)
        xy = vector(p.x, p.y, 0)
        ac  = xy * We²
        g_e + ac
    end
    function gravity_jacobian(p::ecef_pos)::matrix_3x3
        r1 = norm(p)
        r2 = r1*r1
        r3 = r2*r1
        r4 = r3*r1
        m2 = 3*J2*a2/(2*r2)
        m3 = 15*J2*a2*(p.y^2)/(2*r4)
        p2 = p.*p / r2
        
        J  = zeros(3,3)
        matrix_3x3([
        1 - 3*p2[1] + m2*(1-5*p2[1]) - m3*(1-7*p2[1])  p.x*Pe[2]*(-3/r2 - m3*(1-7*p2[3]))             p.x*Pe[3]*(-3/r2 - m3*(1+(2-7*p2[3])))        
        p.x*Pe[2]*(-3/r2 - m3*(1-7*p2[3]))             1 - 3*p2[2] + m2*(1-5*p2[2]) - m3*(1-7*p2[2])  p.y*Pe[3]*(-3/r2 - m3*(1+(2-7*p2[3])))
        p.x*Pe[3]*(-3/r2 - m3*(3-7*p2[3]))             p.y*Pe[3]*(-3/r2 - m3*(3-7*p2[3]))             1 - 3*p2[3] + m2*(1-5*p2[3]) - m3*(3-7*p2[3])
        ])/r3
    end
    function coriolis(v::ecef_vel)::vector # -2(Ωe x Ve)
        vector(2*Ωe*v.y, -2*Ωe*v.x, 0)
    end
    # v = ω × r = [-ry + qz
    #               rx - pz
    #              -qx + py]
    function earth_speed(pe)::ecef_vel  
        ecef_vel(-pe.y*Ωe, pe.x*Ωe, 0) 
    end
    # [elv, yaw]  
    function gravity_anomaly(p::nav_pos)
        p   =  to_ecef(p) 
        g   = -gravity(p)
        m   =  norm(g)
        elv =  asin(g.z/m)
        yaw =  atan(g.y, g.x)
        (elv - p.λ, yaw - p.φ)
    end  
    
    # Navegación inercial en coordenadas ECEF
    function ∫(pe::vector, ve::vector, qeb::quaternion, fb::vector, ωb::vector)
        fe = qeb << fb
        ge = gravity(pe) 
        ae = fe + ge
        cr = coriolis(ve) 

        dp = ve
        dv = ae + cr

        qbe = qeb'
        ωeb = ωb - (qbe << vector(0, 0, Ωe)) 
        dq = Cω(qeb) * ωeb

        (dp, dv, dq)
    end   

end
