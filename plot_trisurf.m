function plot_trisurf ( x , val )

% x: N-by-2 position array with the number of points N. Each row is (x1,x2)
% of a evaluation point. N is the size of data.
% val: field to plot. a column vector of size N

x_min  =  min ( x(:,1) );
x_max  =  max ( x(:,1) );
y_min  =  min ( x(:,2) );
y_max  =  max ( x(:,2) );

L1 = x_max - x_min;
L2 = y_max - y_min;

if L1 >= L2
    f = figure ( 'Position' , [0 0 1.3*500 L2/L1*500]);
else
    f = figure ( 'Position' , [0 0 1.3*L1/L2*500 500]);
end
movegui ( f , 'center' );


t = delaunay(x(:,1), x(:,2));
trisurf(t, x(:,1), x(:,2), val);
colormap(gca,'jet'); shading interp;
set(gca,'TickLabelInterpreter','latex','fontsize',15)
cb = colorbar; cb.TickLabelInterpreter = 'latex';
set(gca,'YDir','normal');
axis equal
grid off
view(0,90);
set(gca,'XColor','none','YColor','none');


% color bar limit
dval    =   max ( val )  -  min ( val ) ;
dval    =   0.1 * dval ;
if dval < 1e-6
    dval = 1e-6;
end
caxis (  [ min(val) , max(val) ]  +  dval * [ -1 , 1 ]  )

drawnow



end