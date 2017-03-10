function [result, mean] = integral2D(planar, nx, x, ny,y, phi)

result = 0;
for i=2:1:nx-1
    
    ri = x(i) - 0.5*(x(i)-x(i-1));
    re = x(i) + 0.5*(x(i+1)-x(i));
    
    if (planar == false)
        A = pi*(re^2-ri^2);
    else
        A = re-ri;
    end
    
    for j=2:1:ny-1
        
        hm = y(j) - 0.5*(y(j)-y(j-1));
        hp = y(j) + 0.5*(y(j+1)-y(j));
        
        V = (hp-hm)*A;
        
        result = result + V*phi(i,j);
        
    end
    
    %North
    deltay = 0.50*(y(ny)-y(ny-1));
    V = A*deltay;
    result = result + V*phi(i,ny);
    
    %South
    deltay = 0.50*(y(2)-y(1));
    V = A*deltay;
    result = result + V*phi(i,1);
    
end

for j=2:1:ny-1
    
    hm = y(j) - 0.5*(y(j)-y(j-1));
    hp = y(j) + 0.5*(y(j+1)-y(j));
    deltay = hp-hm;
    
    %West
    ri = x(1);
    re = x(1)+0.50*(x(2)-x(1));
    A = pi*(re^2-ri^2);
    V = A*deltay;
    result = result + V*phi(1,j);;
    
    %East
    ri = x(nx-1)+0.50*(x(nx)-x(nx-1));
    re = x(nx);
    A = pi*(re^2-ri^2);
    V = A*deltay;
    result = result + V*phi(nx,j);;
end

%Corners S
deltay = 0.50*(y(2)-y(1));

%W
ri = x(1);
re = x(1)+0.50*(x(2)-x(1));
A = pi*(re^2-ri^2);
V = A*deltay;
result = result + V*phi(1,1);
%E
ri = x(nx-1)+0.50*(x(nx)-x(nx-1));
re = x(nx);
A = pi*(re^2-ri^2);
V = A*deltay;
result = result + V*phi(nx,1);

%Corners N
deltay = 0.50*(y(nx)-y(nx-1));

%W
ri = x(1);
re = x(1)+0.50*(x(2)-x(1));
A = pi*(re^2-ri^2);
V = A*deltay;
result = result + V*phi(1,ny);
%E
ri = x(nx-1)+0.50*(x(nx)-x(nx-1));
re = x(nx);
A = pi*(re^2-ri^2);
V = A*deltay;
result = result + V*phi(nx,ny);

if (planar == true)
    Vtot = (x(nx)-x(1))*(y(ny)-y(1));
    mean = result/Vtot;
else
    Vtot = pi*(x(nx)^2-x(1)^2)*(y(ny)-y(1));
    mean = result/Vtot;
end
