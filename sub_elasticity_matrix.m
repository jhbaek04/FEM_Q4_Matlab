function [ D ] = sub_elasticity_matrix ( mat )

lambda  =   mat.Lames_first_parameter;
mu      =   mat.shear_modulus;
M       =   2*mu+lambda;

D       =   [      M , lambda ,      0  ;
              lambda ,      M ,      0  ;
                   0 ,      0 ,     mu  ] ;

end