% coordenadas inerciales esféricas 
classdef ICRF
    properties (Constant)
        pe  = (23+27/60)*pi/180;         % inclinación del eje polar
        Ce  = rot_dcm('y', -ICRF.pe);
        Ds  = 2*pi/WGS84.We;             % día sidéreo 23:56:04.09; 
        Ys  = 365*86400+5*3600+48*60+45; % año
    end
    properties (GetAccess = public, SetAccess = protected)
        ut1;
        Cei;
    end
    methods
        %% ICRF / ECEF map
        function [this] = ICRF(tm)
            this.ut1 = tm;
            rz  = tm * WGS84.We; % rotación terrestre  
            Cz  = rot_dcm('z', rz);
            this.Cei = Cz*ICRF.Ce;
        end       
        % ve : vector ECEF 
        % vc : vector ICRF 
        function vc = to_eci(this, ve)            
            vc  = this.Cei'*ve;
        end
        % pc : posición ICRF 
        % vc : velocidad ICRF 
        % pe : posición ECEF
        % ve : velocidad ECEF 
        function [pe, ve] = to_ecef(this, pc, vc)            
            pe  = this.Cei*pc;
            ve  = WGS84.earth_speed(pe); 
            ve  = this.Cei*vc - ve;
        end
        % pc : posición ICRF 
        % pe : posición ECEF
        function pe = ecef_pos(this, pc)            
            pe  = this.Cei*pc;
        end
        % vc : velocidad ICRF 
        % ve : velocidad ECEF 
        function ve = ecef_vel(this, pe, vc)            
            ve  = WGS84.earth_speed(pe); 
            ve  = this.Cei*vc - ve;
        end
        % qib : quaternion ICRF 
        % qeb : quaternion ECEF 
        function qeb = ecef_att(this, qib)            
            qei = quaternion(this.Cei);
            qeb = qei * qib;
        end
    end
    methods(Static)
        %% Coordenadas esféricas {ρ, θ, φ} ISO 80000-2:2019
        %
        % ρ: radial distance (slant distance to origin)
        % θ: polar angle (angle with respect to positive polar axis)  
        % φ: azimuthal angle (angle of rotation from the initial meridian plane). 
        
        % Radio geocéntrico
        function [R, R2] = Rc(pc)
            R2  = pc'*pc;
            R   = sqrt(R2);
        end
        % pc: posición en coordenadas cartesianas
        % sp: posición en coordenadas esféricas
        %   .rad: radio geocéntrico (ρ)
        %   .lng: (θ)
        %   .lat: (pi/2 - φ)
        %   .Cli: matriz de rotación a coordenadas inerciales locales
        function sp = to_nav(pc)
            sp.rad = ICRF.Rc(pc);
            re     = sqrt(pc(1:2)'*pc(1:2)); % radio ecuatorial            
            sp.lat = atan(pc(3)/re);         % latitud
            sp.lng = atan2(pc(2), pc(1));    % longitud 
            
            Cz  = rot_dcm('z', sp.lng);
            Cy  = rot_dcm('y',-sp.lat);
            sp.Cli = Cy*Cz;
        end
        function [vl, qlb, sp] = to_loc(pc, vc, qib)
            sp.rad = ICRF.Rc(pc);
            re     = sqrt(pc(1:2)'*pc(1:2)); % radio ecuatorial            
            sp.lat = atan(pc(3)/re);         % latitud
            sp.lng = atan2(pc(2), pc(1));    % longitud 
            
            Cz  = rot_dcm('z', sp.lng);
            Cy  = rot_dcm('y',-sp.lat);
            sp.Cli = Cy*Cz;
            
            vl = sp.Cli * vc;
            qli = quaternion(sp.Cli);
            qlb = qli * qib;
        end
        %% Utilidades
        % yyyy, mm, dd
        function JD = julian_date(yyyy, mm, dd)
            JD = 367 * yyyy - floor(1.75 * (yyyy + floor((mm + 9)/12))) ...
               + floor(275/9 * mm) + dd + 1721013.5;            
        end       
        function [vf, g, R] = v_orbit(h)
            R  = WGS84.a + h;
            g  = WGS84.ge * (WGS84.a./R).^2;
            vf = sqrt(g.*R);
        end       
        function h = height(pc)
            h = ICRF.Rc(pc) - WGS84.a;
        end        
        function y = in_orbit(pc, vc)
            sp = ICRF.to_nav(pc);
            vl = sp.Cli * vc;
            vt = sqrt(vl(1:2)'*vl(1:2));
            vo = sqrt(WGS84.ge/sp.rad)*WGS84.a;
            y  = vt >= vo;
        end
        %%
        % h  : alturas
        % r  : rango
        % p  : ángulo de cabeceo en terna geo loc
        % y  : desvío la trayectoria respecto de la tangente a una
        %      orbita circular 
        % a  : ángulo entre empuje y velocidad
        % g  : aceleración gravitacional        
        % R  : radio geocéntrico
        % vb : ground speed in body coordinates
        function [h, r, p, y, a, g, R, vb] = geo_map(t, pc, vc, qc)
            r2 = pc(:,1).^2 + pc(:,2).^2 ;
            R2 = r2 + pc(:,3).^2;
            R  = sqrt(R2);
            h  = R - WGS84.a;
            
            re  = sqrt(r2); % radio ecuatorial            
            lat = atan(pc(:,3)./re);          % latitud
            lng = atan2(pc(:,2), pc(:,1));    % longitud 
            
            
            % Cgi = rot_dcm('y', pi/2-lat)*rot_dcm('z', lng);
            % sy = sym('sy'); cy = sym('cy'); sz = sym('sz'); cz = sym('cz');
            % [cy 0 -sy ; 0 1 0 ; sy 0 cy] * [cz sz 0 ; -sz cz 0 ; 0 0 1] 
            % SEU -> ENU
            % [0 -1 0 ; 1 0 0 ; 0 0 1]
            % ENU
            % Cgi = [   sz   -cz   0
            %        cy*cz cy*sz -sy
            %        cz*sy sy*sz  cy];
            sy = sin(-lat); cy = cos(lat);
            sz = sin( lng); cz = cos(lng);
            u  = vc(:,1); v = vc(:,2); w = vc(:,3);
            vl = [(    sz.*u -     cz.*v) ... 
                  (cy.*cz.*u + cy.*sz.*v - sy.*w) ...
                  (cz.*sy.*u + sy.*sz.*v + cy.*w)]; 
            vt = sqrt(vl(:,1).^2 + vl(:,2).^2);
            y  = atan(vc(:,3)./vt);
            
            
            %%
            vb = zeros(length(t),3);
            r  = qc(:,1);
            x  = qc(:,2); x2 = x.*x;
            y  = qc(:,3); y2 = y.*y;
            z  = qc(:,4); z2 = z.*z;
            vb(:,1) = 2*((0.5-y2-z2).* u + (x.*y-r.*z).* v + (x.*z+r.*y).* w);
            vb(:,2) = 2*((x.*y+r.*z).* u + (0.5-x2-z2).* v + (y.*z-r.*x).* w);
            vb(:,3) = 2*((x.*z-r.*y).* u + (y.*z+r.*x).* v + (0.5-x2-y2).* w);

            
            
            

            p = zeros(size(y));
            a = zeros(size(y));
            g = zeros(size(y));
            r = zeros(size(y));
        end
        
        % ECI - ECEF
        % utc: hora UTC en segundos
%         function ve = ecef_to_eci(vi, utc)
%             ut = ephemeris.UT(utc);
%             r = -WGS84.We * ut;
%             s = sin(r);
%             c = cos(r);
%             if isrow(vi)
%                 vi = vi';
%             end
%             vi = WGS84.Ctl * vi;
%             ve = [vi(1)*c - vi(2)*s  
%                   vi(2)*c + vi(1)*s  
%                   vi(3)];
%         end
%         function vi = eci_to_ecef(ve, utc)
%             ut = ephemeris.UT(utc);
%             r = WGS84.We * ut;
%             s = sin(r);
%             c = cos(r);
%             vi = [ve(1)*c - ve(2)*s  
%                   ve(2)*c + ve(1)*s  
%                   ve(3)];
%         end        
            
%             ry = pi/2 - nav.lat;
%             rz = nav.lng;
%             sy  = sin(ry); cy = cos(ry);
%             sz  = sin(rz); cz = cos(rz);
%             
%             % ENU
%             Cgi = [   sz   -cz   0
%                    cy*cz cy*sz -sy
%                    cz*sy sy*sz  cy];
%             vg  = Cgi * vc;
%             vt  = sqrt(vg(1:2)'*vg(1:2));
%             y   = atan(vg(3)/vt); % inclinación de la velocidad
%             vg  = Cgi * vc;
%             qgi = quaternion(Cgi);
%             qg  = qgi * qib;
            
  
%         function y = in_orbit(pc, vc)
%             
%             [R, R2] = ICRF.Rc(pc);
%             g = ICRF.go * ICRF.a2/R2;
%             
%             re  = sqrt(pc(1:2)'*pc(1:2)); % radio ecuatorial            
%             lat = atan(pc(3)/re);        % latitud
%             lng = atan2(pc(2), pc(1));   % longitud 
%             
%             ry = pi/2 - lat;
%             rz = lng;
%             sy = sin(ry); cy = cos(ry);
%             sz = sin(rz); cz = cos(rz);
%             
%             % ENU
%             Cgi = [   sz   -cz   0
%                    cy*cz cy*sz -sy
%                    cz*sy sy*sz  cy];
%             vg  = Cgi * vc;
%             vt2 = vg(1:2)'*vg(1:2);
%             
%             y = vt2/R >= g;
%             
%         end 
     end
end


