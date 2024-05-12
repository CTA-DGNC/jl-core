module ECI

    module R2
        using  ...geom: radian, degree 
        using  ...geom.r2 #: vector, dcm, rot_mtx 
        import ...polar    
        using  ...nav: WGS84, scalar
        export nav_pos, lgv_vel, gcn_pos, eci_vel, gravitation

        # posición geocéntrica
        const gcn_pos = polar.pos
        # velocidad en terna local
        const lgv_vel = polar.vel

        struct nav_pos
            θ::radian
            h::scalar # altura
        end        

        gcn_pos(p::nav_pos)::gcn_pos = gcn_pos(p.θ, p.h + WGS84.a) 
        nav_pos(p::gcn_pos)::nav_pos = nav_pos(p.θ, p.r - WGS84.a) 
        nav_pos(p::vector )::nav_pos = nav_pos(gcn_pos(p)) 
        eci_vel(θ::radian, ve::scalar)::vector = vector(-ve*sin(θ) , ve*cos(θ)) 
        
        gravitation(p::gcn_pos)::gcn_pos = pos(-p.θ, WGS84.ge * WGS84.a² / p.r^2) 
        gravitation(p::vector )::vector  = -p * (WGS84.ge * WGS84.a² / (p.x^2+p.y^2)^(1.5))    

    end

    module R3
        import ...sphere
        using  ...nav , ...geom
        export gcn_pos, gravitation     

        # posición geocéntrica
        const gcn_pos = sphere.pos        

        gravitation(p::vector)::vector = -p * (WGS84.ge * WGS84.a² / (p.x^2+p.y^2)^(1.5)) 

        # Navegación inercial en coordenadas ECI
        function ∫(pc::vector, vi::vector, qib::quaternion, fb::vector, ωb::vector)
            fi = qib << fb
            gi = gravitation(pc) 
            dp = ve 
            dv = fi + gi
            dq = Cω(qib) * ωb
    
            (dp, dv, dq)
        end 
    end

end
