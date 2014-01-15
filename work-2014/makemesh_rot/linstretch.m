function xpos = linstretch(dx1,n,L)

  ##This function distributes points over an interval, with a linear
  ##increase in the distance between the points

  ##INPUTS:

  ## dx1, the first increment
  ## n, the number of points
  ## L, the distance over which to distribute the points

  ##OUTPUTS:

  ## x, the location of the points

  sumn = 0;
  for k=(1:(n-1))
    sumn = sumn+(k-1);
  endfor

  sumnn = 0;
  for k=(1:(n-1))
    sumnn = sumnn+k;
  endfor
  ddx = dx1/sumnn;

  m = (L-n*dx1)/sumn;
  
  xpos = zeros((n),1);
  dxall = zeros((n-1),1);
  dxall(1) = dx1;

  for k=(2:(n-1))
    dxall(k) = m*(k-1)+dx1 + k*ddx;
  endfor
  
  for k=(2:(n))
    xpos(k)=sum(dxall(1:(k-1)));
  endfor
