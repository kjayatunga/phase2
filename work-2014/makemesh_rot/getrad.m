function [r,Qout,temp,dev] = getrad(x1,x2)

  ## Function to get radius of curvature for an element edge, given the
  ## extent of the element, and the function (surffunc) that defines
  ## where the surface should be

  ## OUTPUTS: 

  ## r = the radius of curvature 

  ## Qout = integrated difference between surface and constant radius curve

  ## temp = whatever

  ## dev = ratio of max deviation of the constant radius curve from a
  ## straight line between the two end points to the distance between
  ## the points (this is what goes in the mesh file

  if(x1>x2)
    x = x1;
    x1 = x2;
    x2 = x;
  endif
  
  y1 = surffunc(x1);
  y2 = surffunc(x2);
  
  func = inline('(circarc(x,rc,x1,x2,y1,y2) - surffunc(x)).^2', 'x','rc','x1','x2','y1','y2');
  
  funQ = inline('quadv(func,x1,x2,tol,trace,rc,x1,x2,y1,y2)', 'func', 'x1','x2', 'tol', 'trace', 'rc', 'y1', 'y2');
  
  ## Need to minimise funQ wrt rc
  rc = 0.5*sqrt((x2-x1)^2 + (y2-y1)^2);
  rc_inc = rc/10;
  rc_base = rc;
  tol = 10^-10;
  
  Q = zeros(3,1);
  Q(1) = funQ(func, x1, x2, tol, 0, rc, y1, y2);
  Q(2) = funQ(func, x1, x2, tol, 0, rc+rc_inc, y1, y2);
  Q(3) = funQ(func, x1, x2, tol, 0, rc+2*rc_inc, y1, y2);
  
  rc_base = rc_base+rc_inc;
  m = 1;
  while(((!((Q(1)>=Q(2))&(Q(3)>Q(2))))&(m<2000))||(m<2))
    Q(1) = Q(2);
    Q(2) = Q(3);
    Q(3) = funQ(func, x1, x2, tol, 0, rc_base+2*rc_inc, y1, y2);
    rc_base = rc_base+rc_inc;
    temp(m,1) = rc_base;
    temp(m,2) = Q(1);
    m++;
  endwhile
  temp(m,1) = rc_base+rc_inc;
  temp(m,2) = Q(2);
  temp(m+1,1) = rc_base+2*rc_inc;
  temp(m+1,2) = Q(3);
  
  ## Now have minimum bracketed. Refine solution
  Qchange = 1;
  
  while(Qchange>0.000001)
    rc_inc = rc_inc/2;
    Qhalf = funQ(func, x1, x2, tol, 0, rc_base+rc_inc, y1, y2);
    Qminold = Q(2);
    
    if((Q(1)>Qhalf) & (Q(2)>Qhalf))
      ## Minimum in the first half of the interval
      Q(3) = Q(2);
      Q(2) = Qhalf;
    else
      ## Minimum in the second half of the interval
      Q(1) = Q(2);
      Q(2) = funQ(func, x1, x2, tol, 0, rc_base+3*rc_inc, y1, y2);
      rc_base = rc_base+2*rc_inc;
    endif
    Qchange = abs(Qminold - Q(2))/Qminold;
  endwhile
  
  r = rc_base+rc_inc;
  Qout = funQ(func, x1, x2, tol, 0, r, y1, y2);

  ## Calculate dev
  [y,xc,yc,xdev] = circarc(x1,r,x1,x2,y1,y2);
  dev = xdev/sqrt((x2-x1)^2+(y2-y1)^2);
  
end

