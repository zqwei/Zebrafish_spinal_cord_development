function drawArrow(x,y, z,props) 
    quiver3( x(1),y(1), z(1),x(2)-x(1),y(2)-y(1), z(2)-z(1),0, props{:} );
end