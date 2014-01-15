function y = surffunc(x)

  ## Function describing surface
  
  ## Ellipse
  asq = 1/4;
  bsq = 1/16;
  
  y = sqrt(bsq-bsq.*x.^2./asq);

  ## r = sqrt((a*b)./(b*cos(theta).^2 + a*sin(theta).^2));
  ## x = r.*cos(theta);
  ## y = r.*sin(theta);
  
end

