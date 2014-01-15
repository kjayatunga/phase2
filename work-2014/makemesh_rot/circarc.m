## Function describing the circular arc of one edge

function [y,xc,yc,xdev] = circarc(x,rc,x1,x2,y1,y2)

## OUTPUTS

##  y = position on arc at x

##  xc,yc = coordinates of centre of curvature

##  xdev = max deviation of arc from straight line joining points 1 and 2

  d = sqrt((x2-x1)^2+(y2-y1)^2);
  l = sqrt(rc^2-d^2/4);
  xdev = abs(rc-l);

  if(x1>x2)
    temp = x1;
    x1 = x2;
    x2 = temp;
    temp = y1;
    y1 = y2;
    y2 = temp;
  endif

  xc = x1+(x2-x1)/2 + l*cos(atan2((x2-x1),(y2-y1)));
  yc = y1+(y2-y1)/2 - l*sin(atan2((x2-x1),(y2-y1)));
  y = sqrt(rc^2 - (x-xc).^2) + yc;
  
# function [y,xc,yc] = circarc(x1,x2,y1,y2,rc,x)

# np = length(x);
# for n=(1:np)
#   if((x(n)<x1&x(n)<x2)||(x(n)>x1&x(n)>x2))
#     printf("x(%d) outside element edges: stopping\n", n);
#     break;
#   else
#     xc = x1+(x2-x1)/2 + rc*cos(atan2((x2-x1),(y2-y1)));
#     yc = y1+(y2-y1)/2 - rc*sin(atan2((x2-x1),(y2-y1)));
#     l = sqrt(rc^2 + (0.5*(x2-x1))^2 + (0.5*(y2-y1))^2);
#     y(n) = sqrt(l^2 - (x(n)-xc)^2) + yc;
#   endif
# endfor

