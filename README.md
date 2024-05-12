# DGNC
## Dinámica, Guiado, Navegación y Control

- mtlb: código Matlab reutilizable
    - geom [tipos, estructuras y métodos para cómputo geométrico y manejo de rotaciones]
        - degree : tipo primitivo ángulo en grados sexagecimales
        - radian : tipo primitivo ángulo en radianes
        - euler_angles : ángulos de euler  
        - quaternion   : cuaternión 
        - dcm y rot_mtx : matrices de rotación 
        + R2 [especializaciones para geometría en el plano] 

    - nav [tipos, estructuras y métodos para cómputos de navegación]
        - xxx_pos : estructura para datos de posición
        - xxx_vel : estructura para datos de velocidad
        + polar : manejo de coordenadas polares (en R²)
        + sphere: manejo de coordenadas esféricas (en R³)
        + WGS84 : conversiones y datos con el modelo WGS84 
        + ECI : coordenadas inerciales con origen arbitrario
            - R2 : coordenadas inerciales en un plano xy R²
            - R3 : coordenadas inerciales en R³
        


___

## CASO DE ESTUDIO:  

***




