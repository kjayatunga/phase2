## File to make bluff-body meshes, where bodies are cylindrical, or have
## a requirement to be rotated, or something like that

## General mesh structure will be body in the centre of a circle, circle
## in a square, that matches with upstream and downstream rectangular
## mesh. This is all imbedded in a lo-res outer mesh

## Generate surface points of body. Order points anti-clockwise, first
## point closest to centre axis downstream. Number of surface points
## must be divisible by 8

## REQUIRED INPUTS CONTROLLING DOMAIN SIZE AND MESH DENSITY AND WHETHER
## TO ADD CURVATURE TO THE SURFACE: EDIT AS NEEDED

## Radius of circle surrounding the body
rcircle = 3.0;

## Half the side length of the square surrounding the circle
lsquare = 6.0;

## Control mesh between the body and the circle
bl_inc = 0.25;   # thickness of first boundary layer of elements
nbclayers = 6;  # number of layers of elements from body to circle

## Control mesh between the circle and the square
ncslayers = 6;  # number of layers of elements from circle to square

## Mesh downstream of the square
loutlet = 50;   # length of the domain downstream of square rear edge
noutlet = 30;   # number of elements from square rear edge to outlet

## Mesh upstream of the square
linlet = 14;    # length of the domain upstream of the square front edge
ninlet = 3;     # number of elements from square front edge to inlet

## Coarse mesh to either side
ltrans = 24;    # length of the domain transverse of the square
				# transverse edge
ntrans = 2;     # number of elements from square transverse edge to
				# transverse boundary

## Dictate whether curvature is required
curves = 0      # 0 = no curves, 1 = curves

## END OF REQUIRED INPUTS

## PUT CODE TO GENERATE SURFACE POINTS HERE

## Square
# nsurf = 24;
# xsurf = zeros(nsurf,1);
# ysurf = zeros(nsurf,1);

# for n=(1:nsurf/8)
#   xsurf(n) = 0.5;
#   ysurf(n) = (n-1)*4/nsurf;
# endfor

# m=1;
# for n=(nsurf/8+1:3*nsurf/8)
#   xsurf(n) = 0.5-(m-1)*4/nsurf;
#   ysurf(n) = 0.5;
#   m++;
# endfor

# m=1;
# for n=(3*nsurf/8+1:5*nsurf/8)
#   xsurf(n) = -0.5;
#   ysurf(n) = 0.5-(m-1)*4/nsurf;
#   m++;
# endfor

# m=1;
# for n=(5*nsurf/8+1:7*nsurf/8)
#   xsurf(n) = -0.5+(m-1)*4/nsurf;
#   ysurf(n) = -0.5;
#   m++;
# endfor

# m=1;
# for n=(7*nsurf/8+1:nsurf)
#   xsurf(n) = 0.5;
#   ysurf(n) = -0.5+(m-1)*4/nsurf;
#   m++;
# endfor

## Square with an angle of attack
#aoa_degrees=0
#alpha = -((aoa_degrees)/(180))*pi
#nsurf = 16;
#xsurf = zeros(nsurf,1);
#ysurf = zeros(nsurf,1);

#for n=(1:nsurf/8)
# xsurf(n) = 0.5;
#  ysurf(n) = (n-1)*4/nsurf;
#endfor

#m=1;
#for n=(nsurf/8+1:3*nsurf/8)
#  xsurf(n) = 0.5-(m-1)*4/nsurf;
#  ysurf(n) = 0.5;
#  m++;
#endfor

#m=1;
#for n=(3*nsurf/8+1:5*nsurf/8)
#  xsurf(n) = -0.5;
#  ysurf(n) = 0.5-(m-1)*4/nsurf;
#  m++;
#endfor

#m=1;
#for n=(5*nsurf/8+1:7*nsurf/8)
#  xsurf(n) = -0.5+(m-1)*4/nsurf;
#  ysurf(n) = -0.5;
#  m++;
#endfor

#m=1;
#for n=(7*nsurf/8+1:nsurf)
#  xsurf(n) = 0.5;
#  ysurf(n) = -0.5+(m-1)*4/nsurf;
#  m++;
#endfor

## Scale the square to keep the new frontal width = 1
#lside = 1/(sin(alpha)+cos(alpha));
#xsurf = xsurf.*lside;
#ysurf = ysurf.*lside;

## Rotate square by alpha
#phi_orig = atan2(ysurf,xsurf);
#phi = atan2(ysurf,xsurf)+alpha;
#for n=(1:nsurf)
#  if(phi(n)>(pi))
#    phi(n) = phi(n) - 2*pi;
#  endif
#  if(phi(n)<=(-pi))
#    phi(n) = phi(n) + 2*pi;
#  endif
#endfor
#r = sqrt(xsurf.^2+ysurf.^2);
#xsurfnew = r.*cos(phi);
#ysurfnew = r.*sin(phi);

## Find the new point closest to the x-axis
#[a,b] = min(abs(phi));

#for n=(1:nsurf)
#  if((n)<=((nsurf-b)+1))
#    xsurf(n) = xsurfnew(n+b-1);
#    ysurf(n) = ysurfnew(n+b-1);
#  else
#    xsurf(n) = xsurfnew(n-(nsurf-b+1));
#    ysurf(n) = ysurfnew(n-(nsurf-b+1));
#  endif
#endfor

## Rectangle
# nsurf = 24;
# xsurf = zeros(nsurf,1);
# ysurf = zeros(nsurf,1);

# aspect = 2; # Ratio of length to thickness

# l_quarter = (aspect + 1)/2;
# delta_l_ideal = l_quarter/(nsurf/4);
# nshort = floor(1/(2*delta_l_ideal));

# for n=(1:nshort)
#   xsurf(n) = aspect/2;
#   ysurf(n) = (n-1)/(2*nshort);
# endfor

# m=1;
# for n=(nshort+1:nsurf/2-nshort)
#   xsurf(n) = aspect/2 - (m-1)*(aspect/2)/((nsurf/4)-nshort);
#   ysurf(n) = 0.5;
#   m++;
# endfor

# m=1;
# for n=(nsurf/2-nshort+1:nsurf/2+nshort)
#   xsurf(n) = -aspect/2;			# x constan negative from 0 axis 				
#   ysurf(n) = 0.5-(m-1)/(2*nshort);            # starts from the point (1,0.5)
#   m++;
# endfor

# m=1;
# for n=(nsurf/2+nshort+1:nsurf-nshort)
#   xsurf(n) = -aspect/2 + (m-1)*(aspect/2)/((nsurf/4)-nshort);
#   ysurf(n) = -0.5;;
#   m++;
# endfor

# m=1;
# for n=(nsurf-nshort+1:nsurf)
#   xsurf(n) = aspect/2;                    # x constant opsitive from 0 axis
#   ysurf(n) = -0.5 + (m-1)/(2*nshort);
#   m++;
# endfor

## Diamond
# nsurf = 24;
# xsurf = zeros(nsurf,1);
# ysurf = zeros(nsurf,1);

# for n=(1:nsurf/4)
#   xsurf(n) = 0.5 - (n-1)*0.5/(nsurf/4);
#   ysurf(n) = -xsurf(n) + 0.5;
# endfor

# m=1;
# for n=(nsurf/4+1:nsurf/2)
#   xsurf(n) = -(m-1)*0.5/(nsurf/4);
#   ysurf(n) = xsurf(n) + 0.5;
#   m++;
# endfor

# m=1;
# for n=(nsurf/2+1:3*nsurf/4)
#   xsurf(n) = -0.5 + (m-1)*0.5/(nsurf/4);
#   ysurf(n) = -xsurf(n) - 0.5;
#   m++;
# endfor

# m=1;
# for n=(3*nsurf/4+1:nsurf)
#   xsurf(n) = (m-1)*0.5/(nsurf/4);
#   ysurf(n) = xsurf(n) - 0.5;
#   m++;
# endfor

## Rear-pointing equilateral triangle, front face at x=0
aoa_degrees=0
alpha = -((aoa_degrees)/(180))*pi
 nsurf = 24;
 xsurf = zeros(nsurf,1);
 ysurf = zeros(nsurf,1);

 for n=(1:nsurf/3)
   xsurf(n) = 1/sqrt(2) - (n-1)*0.5/(nsurf/4);
   ysurf(n) = -xsurf(n)/sqrt(2)+0.5;
 endfor

 m=1;
 for n=(nsurf/3+1:2*nsurf/3)
   xsurf(n) = 0;
   ysurf(n) = 0.5-(m-1)*3/nsurf;
   m++;
 endfor

 m=1;
 for n=(2*nsurf/3+1:nsurf)
   xsurf(n) = (m-1)*0.5/(nsurf/4);
   ysurf(n) = xsurf(n)/sqrt(2)-0.5;
   m++;
 endfor

## Scale the square to keep the new frontal width = 1
lside = 1/(sin(alpha)+cos(alpha));
xsurf = xsurf.*lside;
ysurf = ysurf.*lside;

## Rotate square by alpha
phi_orig = atan2(ysurf,xsurf);
phi = atan2(ysurf,xsurf)+alpha;
for n=(1:nsurf)
  if(phi(n)>(pi))
    phi(n) = phi(n) - 2*pi;
  endif
  if(phi(n)<=(-pi))
    phi(n) = phi(n) + 2*pi;
  endif
endfor
r = sqrt(xsurf.^2+ysurf.^2);
xsurfnew = r.*cos(phi);
ysurfnew = r.*sin(phi);



# Find the new point closest to the x-axis
[a,b] = min(abs(phi));

for n=(1:nsurf)
  if((n)<=((nsurf-b)+1))
    xsurf(n) = xsurfnew(n+b-1);
    ysurf(n) = ysurfnew(n+b-1);
  else
    xsurf(n) = xsurfnew(n-(nsurf-b+1));
    ysurf(n) = ysurfnew(n-(nsurf-b+1));
  endif
endfor

## Circular cylinder
# nsurf = 16;
# theta = linspace(0,(2*pi*(1-(1/nsurf))),nsurf);
# xsurf = cos(theta);
# ysurf = sin(theta);

## Elliptical cylinder
# a = 0.5;
# #b = 0.05;
# nsurf = 24;
# theta = linspace(0,(2*pi*(1-(1/nsurf))),nsurf);
# xsurf = a*cos(theta);
# #ysurf = b*sin(theta);
# ysurf = surffunc(xsurf);
# for n=(nsurf/2+1:nsurf)
#   ysurf(n) = -ysurf(n);
# endfor

# ## Calculate curve deviation for surface elements of ellipse
# dev = zeros(nsurf,1);
# for n=(1:nsurf-1)
#   [r,Q,temp,dev(n)] = getrad(xsurf(n),xsurf(n+1));
# endfor
# n=nsurf
#   [r,Q,temp,dev(n)] = getrad(xsurf(n),xsurf(1));
###########################################################################################################
## KJ trying out to add surface points manually
# nsurf = 16;
# xsurf = [0.25, 0,  -0.125, -0.25, -0.375, -0.5, -0.5, -0.5, -0.5,  -0.5, -0.375, -0.25, -0.125, 0,   0.25, 0]';
# ysurf = [0.25, 0.5, 0.5,    0.5,   0.5,   0.5,  0.25,  0,   -0.25, -0.5, -0.5,   -0.5,  -0.5,   -0.5,-0.25,0.5]';




###########################################################################################################
## END OF SECTION DEFINING SURFACE POINTS AND REQUIRED ELEMENT CURVATURE

## Generate circle surrounding body
theta = linspace(0,(2*pi*(1-(1/nsurf))),nsurf);
xcircle = rcircle.*cos(theta);
ycircle = rcircle.*sin(theta);

## Generate square surrounding circle
nsquareside = nsurf/4;
lsquare_inc = 2*lsquare/nsquareside;
xsquare = zeros(nsurf,1);
ysquare = zeros(nsurf,1);
m=1;
for n=(1:nsurf/8)
  xsquare(m) = lsquare;
  ysquare(m) = (n-1)*lsquare_inc;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = lsquare-(n-1)*lsquare_inc;
  ysquare(m) = lsquare;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = -(n-1)*lsquare_inc;
  ysquare(m) = lsquare;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = -lsquare;
  ysquare(m) = lsquare-(n-1)*lsquare_inc;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = -lsquare;
  ysquare(m) = -(n-1)*lsquare_inc;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = -lsquare+(n-1)*lsquare_inc;
  ysquare(m) = -lsquare;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = (n-1)*lsquare_inc;
  ysquare(m) = -lsquare;
  m++;
endfor
for n=(1:nsurf/8)
  xsquare(m) = lsquare;
  ysquare(m) = -lsquare+(n-1)*lsquare_inc;
  m++;
endfor

## Generate points in between these shells

## Between body and circle. Points ordered starting closest to body, on
## rays from surface to circle, ray closest to centre axis downstream
## first, i.e.; ray 1 from in to out, then ray 2 from in to out, etc.

xbc = zeros(nsurf*(nbclayers-1),1);
ybc = zeros(nsurf*(nbclayers-1),1);

for n=(1:nsurf)
  lray = sqrt((ycircle(n)-ysurf(n))^2+(xcircle(n)-xsurf(n))^2);
  theta = atan2((ycircle(n)-ysurf(n)),(xcircle(n)-xsurf(n)));
  xray = linstretch(bl_inc,nbclayers+1,lray);
  for m=(1:nbclayers-1)
    xbc((n-1)*(nbclayers-1)+m) = xsurf(n)+xray(m+1)*cos(theta);
    ybc((n-1)*(nbclayers-1)+m) = ysurf(n)+xray(m+1)*sin(theta);
  endfor
endfor

## Between circle and square. Points ordered starting closest to body,
## on rays from surface to circle, ray closest to centre axis downstream
## first, i.e.; ray 1 from in to out, then ray 2 from in to out, etc.

xcs = zeros(nsurf*(ncslayers-1),1);
ycs = zeros(nsurf*(ncslayers-1),1);

for n=(1:nsurf)
  lray = sqrt((ysquare(n)-ycircle(n))^2+(xsquare(n)-xcircle(n))^2);
  theta = atan2((ysquare(n)-ycircle(n)),(xsquare(n)-xcircle(n)));
  xray = linspace(0,lray,ncslayers+1);
  for m=(1:ncslayers-1)
    xcs((n-1)*(ncslayers-1)+m) = xcircle(n)+xray(m+1)*cos(theta);
    ycs((n-1)*(ncslayers-1)+m) = ycircle(n)+xray(m+1)*sin(theta);
  endfor
endfor

## Now generate mesh outside of the body box

## Upstream
linlet_inc1 = lsquare-xcs(ncslayers-1);
xinlet = -lsquare-linstretch(linlet_inc1,ninlet+1,linlet);

## Downstream
xoutlet = lsquare + linstretch(linlet_inc1,noutlet+1,loutlet);

## Transverse
#ltrans_inc1 = lsquare-xcs(ncslayers-1);
ltrans_inc1 = lsquare_inc;
if(ntrans!=1)
  ytrans = lsquare + linstretch(ltrans_inc1,ntrans+1,ltrans);
else
  ytrans = lsquare+ltrans;
endif

## Output node locations and connectivity to mesh file. Build elements
## in layers, number points in layers, anticlockwise, for internal mesh.
## Then do downstream, then upstream, then transverse.

fid = fopen("mesh.file","w");

## Node locations
nnodes = length(xsurf) + length(xbc) + length(xcircle) + length(xcs) + \
    length(xsquare) + (nsurf/4+1)*noutlet + (nsurf/4+1)*ninlet \
    +2*(ninlet+(nsurf/4+1)+noutlet)*ntrans;

fprintf(fid,"%d nodes\n", nnodes);

m=1;
## Surface nodes
for n=(1:nsurf)
  fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xsurf(n), ysurf(n), 0.0, m);
  m++;
endfor
## Nodes on layers between surface and circle
for k=(1:(nbclayers-1))
  for n=(k:nbclayers-1:(nbclayers-1)*nsurf)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xbc(n), ybc(n), 0.0, m);
    m++;
  endfor
endfor
## Nodes on circle
for n=(1:nsurf)
  fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xcircle(n), ycircle(n), 0.0, m);
  m++;
endfor
## Nodes on layers between circle and square
for k=(1:(ncslayers-1))
  for n=(k:ncslayers-1:(ncslayers-1)*nsurf)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xcs(n), ycs(n), 0.0, m);
    m++;
  endfor
endfor
## Nodes on square
for n=(1:nsurf)
  fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xsquare(n), ysquare(n), 0.0, m);
  m++;
endfor

## Nodes in downstream mesh. Increasing downstream, then up.
for k=(7*(nsurf/8)+1:nsurf)
  for n=(1:noutlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xoutlet(n+1), ysquare(k), 0.0, m);
    m++;
  endfor
endfor
for k=(1:nsurf/8)
  for n=(1:noutlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xoutlet(n+1), ysquare(k), 0.0, m);
    m++;
  endfor
endfor
k=nsurf/8+1;
for n=(1:noutlet)
  fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xoutlet(n+1), ysquare(k), 0.0, m);
  m++;
endfor

## Nodes in upstream mesh. Increasing upstream, then up.
for k=(7*(nsurf/8)+1:nsurf)
  for n=(1:ninlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xinlet(n+1), ysquare(k), 0.0, m);
    m++;
  endfor
endfor
for k=(1:nsurf/8)
  for n=(1:ninlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xinlet(n+1), ysquare(k), 0.0, m);
    m++;
  endfor
endfor
k=nsurf/8+1;
for n=(1:ninlet)
  fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xinlet(n+1), ysquare(k), 0.0, m);
  m++;
endfor

## Nodes in positive transverse mesh. Increasing downstream, then up.
for k=(2:ntrans+1)
  for n=(1:ninlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xinlet(ninlet+2-n), ytrans(k), 0.0, m);
    m++;
  endfor
  for n=(1:nsurf/4+1)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", -lsquare+(n-1)*lsquare_inc, ytrans(k), 0.0, m);
    m++;
  endfor
  for n=(1:noutlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xoutlet(n+1), ytrans(k), 0.0, m);
    m++;
  endfor
endfor

## Nodes in negative transverse mesh. Increasing downstream, then down.
for k=(2:ntrans+1)
  for n=(1:ninlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xinlet(ninlet+2-n), -ytrans(k), 0.0, m);
    m++;
  endfor
  for n=(1:nsurf/4+1)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", -lsquare+(n-1)*lsquare_inc, -ytrans(k), 0.0, m);
    m++;
  endfor
  for n=(1:noutlet)
    fprintf(fid," %9.4f %9.4f %9.4f %5d\n", xoutlet(n+1), -ytrans(k), 0.0, m);
    m++;
  endfor
endfor

## Element connectivity
nelem = (nbclayers+ncslayers)*nsurf + (nsurf/4)*(noutlet+ninlet) + \
    2*(nsurf/4+noutlet+ninlet)*ntrans;

fprintf(fid,"%d elements\n", nelem);

k=1;
for m=(1:(nbclayers+ncslayers))
  if(m==1)
    ## First row: Add boundary condition to body surface
    for n=(1:nsurf-1)
      fprintf(fid, "%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
  	      k, (m-1)*nsurf+n, m*nsurf+n, m*nsurf+n+1, (m-1)*nsurf+n+1, \
  	      0, 0, 0, 2);
      k++;
    endfor
    n=nsurf;
    fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
  	    k, (m-1)*nsurf+n, m*nsurf+n, m*nsurf+1, (m-1)*nsurf+1, \
  	    0, 0, 0, 2);
    k++;
  else
    for n=(1:nsurf-1)
      fprintf(fid, "%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, (m-1)*nsurf+n, m*nsurf+n, m*nsurf+n+1, (m-1)*nsurf+n+1, \
	      0, 0, 0, 0);
      k++;
    endfor
    n=nsurf;
    fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	    k, (m-1)*nsurf+n, m*nsurf+n, m*nsurf+1, (m-1)*nsurf+1, \
	    0, 0, 0, 0);
    k++;
  endif
endfor

## Downstream elements
nprev = nsurf*(nbclayers+ncslayers+1); #Number of nodes previous to the
				#current section

mmatch = zeros(nsurf/4+1,1); #Node numbers of surface of inner mesh that
				#needs to be stitched to downstream mesh
for n=(1:nsurf/8)
  mmatch(n) = nprev-nsurf/8+n;
endfor
m=1;
for n=(nsurf/8+1:nsurf/4+1)
  mmatch(n) = nprev-nsurf+m;
  m++;
endfor

m=1;
for m=(1:nsurf/4)
  for n=(1:noutlet)
    if(n==1)
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, mmatch(m), nprev+(m-1)*noutlet+1, nprev+m*noutlet+1, mmatch(m+1), 0, 0, 0, 0);
      k++;
    elseif(n==noutlet)
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+(m-1)*noutlet+n-1, nprev+(m-1)*noutlet+n, nprev+m*noutlet+n, nprev+m*noutlet+n-1, 0, 4, 0, 0);
      k++;
    else
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+(m-1)*noutlet+n-1, nprev+(m-1)*noutlet+n, nprev+m*noutlet+n, nprev+m*noutlet+n-1, 0, 0, 0, 0);
      k++;
    endif
  endfor
endfor

# ## Upstream elements
nprev = nsurf*(nbclayers+ncslayers+1) + (nsurf/4+1)*noutlet; #Number of
				#nodes previous to the current section

mmatch = zeros(nsurf/4+1,1); #Node numbers of surface of inner mesh that
				#needs to be stitched to downstream mesh
for n=(1:nsurf/4+1)
  mmatch(n) = nsurf*(nbclayers+ncslayers+1)-3*nsurf/8+2-n;
endfor

m=1;
for m=(1:nsurf/4)
  for n=(1:ninlet)
    if(n==1)
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, mmatch(m), mmatch(m+1), nprev+m*ninlet+1, nprev+(m-1)*ninlet+1, 0, 0, 0, 0);
      k++;
    elseif(n==ninlet)
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+(m-1)*ninlet+n-1, nprev+m*ninlet+n-1, nprev+m*ninlet+n, nprev+(m-1)*ninlet+n, 0, 0, 3, 0);
      k++;
    else
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+(m-1)*ninlet+n-1, nprev+m*ninlet+n-1, nprev+m*ninlet+n, nprev+(m-1)*ninlet+n, 0, 0, 0, 0);
      k++;
    endif
  endfor
endfor

## Positive transverse elements
nprev = nsurf*(nbclayers+ncslayers+1) + \
    (nsurf/4+1)*noutlet + (nsurf/4+1)*ninlet; #Number of nodes previous
				#to the current section

ntransrow = ninlet+nsurf/4+1+noutlet;

mmatch = zeros(ntransrow,1);
for n=(1:ninlet)
  mmatch(n) = nprev+1-n;
endfor
m=1;
for n=(ninlet+1:ninlet+nsurf/4+1)
  mmatch(n) = nsurf*(nbclayers+ncslayers+1)-5*nsurf/8+2-m;
  m++;
endfor
m=1;
for n=(ninlet+nsurf/4+2:ninlet+nsurf/4+1+noutlet)
  mmatch(n) = nsurf*(nbclayers+ncslayers+1)+(nsurf/4)*noutlet+m;
  m++;
endfor

bc = zeros(4,1);

for m=(1:ntrans)
  for n=(1:ntransrow-1)
    ## Set tags for boundary conditions
    bc(1) = bc(2) = bc(3) = bc(4) = 0;
    if(m==ntrans)
      bc(3) = 3;
    endif
    if(n==1)
      bc(4) = 3;
    endif
    if(n==ntransrow-1)
      bc(2) = 4;
    endif
    
    if(m==1)
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, mmatch(n), mmatch(n+1), nprev+n+1, nprev+n, bc(1), bc(2), bc(3), bc(4));
      k++;
    else
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+(m-2)*ntransrow+n, nprev+(m-2)*ntransrow+n+1, nprev+(m-1)*ntransrow+n+1, nprev+(m-1)*ntransrow+n, bc(1), bc(2), bc(3), bc(4));
      k++;
    endif
  endfor
endfor

## Negative transverse elements
nprev = nsurf*(nbclayers+ncslayers+1) + \
    (nsurf/4+1)*noutlet + (nsurf/4+1)*ninlet +\
    ntransrow*ntrans; #Number of nodes previous to the current section

ntransrow = ninlet+nsurf/4+1+noutlet;

mmatch = zeros(ntransrow,1);
for n=(1:ninlet)
  mmatch(n) = nsurf*(nbclayers+ncslayers+1) + \
    (nsurf/4+1)*noutlet + ninlet + 1 -n;
endfor
m=1;
for n=(ninlet+1:ninlet+nsurf/4+1)
  mmatch(n) = nsurf*(nbclayers+ncslayers+1)-3*nsurf/8+m;
  m++;
endfor
m=1;
for n=(ninlet+nsurf/4+2:ninlet+nsurf/4+1+noutlet)
  mmatch(n) = nsurf*(nbclayers+ncslayers+1)+m;
  m++;
endfor

for m=(1:ntrans)
  for n=(1:ntransrow-1)
    ## Set tags for boundary conditions
    bc(1) = bc(2) = bc(3) = bc(4) = 0;
    if(m==ntrans)
      bc(1) = 3;
    endif
    if(n==1)
      bc(4) = 3;
    endif
    if(n==ntransrow-1)
      bc(2) = 4;
    endif

    if(m==1)
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+n,  nprev+n+1, mmatch(n+1), mmatch(n), bc(1), bc(2), bc(3), bc(4));
      k++;
    else
      fprintf(fid,"%5d %5d %5d %5d %5d %5d %5d %5d %5d\n", \
	      k, nprev+(m-1)*ntransrow+n, nprev+(m-1)*ntransrow+n+1, nprev+(m-2)*ntransrow+n+1, nprev+(m-2)*ntransrow+n,  bc(1), bc(2), bc(3), bc(4));
      k++;
    endif
  endfor
endfor

## Write out information for curvature of surface elements
if(curves != 0)
  fprintf(fid,"%d curves\n", nsurf);
  for n=(1:nsurf)
    fprintf(fid,"%5d %5d %10.6f %10.6f\n", n, 2, 0, -dev(n));
  endfor
endif

fclose(fid);
