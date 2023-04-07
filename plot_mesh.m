function plot_mesh ( x_node , connectivity )

x_min  =  min ( x_node(:,1) );
x_max  =  max ( x_node(:,1) );
y_min  =  min ( x_node(:,2) );
y_max  =  max ( x_node(:,2) );

L1 = x_max - x_min;
L2 = y_max - y_min;

if L1 >= L2
    f = figure ( 'Position' , [0 0 500 L2/L1*500]);
else
    f = figure ( 'Position' , [0 0 L1/L2*500 500]);
end
movegui ( f , 'center' );


x = x_node(:,1);
y = x_node(:,2);

patch ( x(connectivity) , y(connectivity) , [.9 .9 .9] ), hold on
plot  ( x , y , '.', 'MarkerSize',15 , 'Color',[.1 .4 .7] )

axis equal
drawnow


end