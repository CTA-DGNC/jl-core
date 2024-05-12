module sphere

    using ..geom
    using ..nav

    export pos , geo_ref

    struct pos
        λ::radian
        φ::radian
        r::scalar 
    end

    function pos(p::vector)::pos
        re2 = p.x^2 + p.y^2      # radio ecuatorial            
        pos(atan(pc.z/ √re2), atan(p.y, p.x), sqrt(re2+p.z^2))
    end        

    #=
    Conversiones entre actitud en terna ECEF y terna geográfica local
    ENU/NED
    
    ECEF: Earth Centered Earth Fixed
    ENU : (g/u) terna geográfica { east north up }
    NED : (l/n) terna local { north east down } 
    =#
    
    # q^e_g: rotación desde terna ENU a ECEF  
    # Cuartenion para proyectar vectores ECEF a terna ENU
    function enu_quat(p::geo_pos)::quaternion
        sl = sin(p.φ); cl = cos(p.φ)
        sq = sin(p.λ); cq = cos(p.λ)
        r  = 0.5*sqrt(1-sl-sq.*sl+sq)
        ni = 0.25./r
        quaternion(r, [-ni*cq*(sl-1), ni*cq*cl, ni*cl*(sq+1)])
    end        
    # C^u_e 
    function enu_dcm(p::geo_pos)::dcm
        sz = sin(p.φ); cz = cos(p.φ)
        sy = sin(p.λ); cy = cos(p.λ)
        dcm([-sz     cz     0 
             -sy*cz -sy*sz  cy
              cy*cz  cy*sz  sy])
    end         
    #= -----------------------------------------------------------------
        q^e_n: rotación desde terna NED a ECEF  
        Cuartenion para proyectar vectores ECEF a terna NED
    =#
    function ned_quat(p::geo_pos)::quaternion
        sl = sin(p.φ/2); cl = cos(p.φ/2)
        sq = sin(p.λ/2); cq = cos(p.λ/2)
        r2 = 0.5*sqrt(2);
        quaternion(r2*cl*(cq-sq), [r2*sl*(cq+sq), -r2*cl*(cq+sq), r2*sl*(cq-sq)])         
    end        
    # C^n_e 
    function ned_dcm(p::geo_pos)::dcm
        sz = sin(p.φ); cz = cos(p.φ);
        sy = sin(p.λ); cy = cos(p.λ);
        dcm([-sy*cz -sy*sz  cy
             -sz     cz     0 
             -cy*cz -cy*sz -sy])
    end
    #= -----------------------------------------------------------------
        C = rot_mtx('y', ll1(1) - ll2(1)) * rot_mtx('z', ll2(2) - ll1(2));
        tr = sum(diag(C))
        cos(q/2) = 1/2 * sqrt(1+tr)
    
        multiplicar q por el radio de la esfera
    =#
    function distance(p1::geo_pos, p2::geo_pos)
        lat = p2.λ - p1.λ
        lng = p2.φ - p1.φ
        cy  = cos(lat)
        cz  = cos(lng)
        tr  = cy.*cz + cz + cy
        acos(0.5*sqrt(1+tr))*2
    end 

    # ------------------------------------------------------------------------- 
    # rotación entre terna ENU (g) y terna NED (l)
    const C_lg = dcm([0  1  0 ; 
                      1  0  0 ; 
                      0  0 -1])
    const q_nu = quaternion(0, 0.5*sqrt(2)*[-1, -1, 0])      
    # Referencia geográfica para mapear coordenadas cartesianas
    struct geo_ref
        p::geo_pos
        q::quaternion # rotación NED → ECEF
    end

    geo_ref(p::geo_pos)::geo_ref = geo_ref(p, ned_quat(p))    
    # ECEF a NED    
    Base.:<<(r::geo_ref, q::quaternion)::quaternion = r.q' * q 
    # NED a ECEF 
    Base.:>>(r::geo_ref, q::quaternion)::quaternion = r.q * q 
    #
    to_enu(q_ned::quaternion)::quaternion = q_ned * q_nu
    to_ned(q_enu::quaternion)::quaternion = q_enu * q_nu
    enu_ned(v::vector)::vector = vector(v.y, v.x, -v.z)

end

