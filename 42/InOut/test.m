R=6378145;
h=500000;
ratio = (R+h)/R;
thetaMax = asin(1/ratio);
halfFOV = 140/2*pi/180;
theta = 0;

if (theta > halfFOV)
		theta1 = theta - halfFOV;
		theta2 = Limit_max(theta + halfFOV,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 - alpha1;
elseif (theta < halfFOV)
		theta1 = Limit_max(halfFOV - theta,thetaMax);
		theta2 = Limit_max(halfFOV + theta,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 + alpha1;
else
		theta1 = theta;
		theta2 = Limit_max(halfFOV + theta,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 + alpha1;
end
	swathWidth = alpha/(2*pi)*R;
	swathArea = 2*pi*R^2*(1-cos(alpha/2));
   
function xout=Limit_max(x, max)
if(x > max)
    xout = max;
else
    xout = x;
end
end