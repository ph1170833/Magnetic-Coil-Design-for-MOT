tspan = 0:1e-5:1e-2;
clf
x0 = 0; %Initial Position
Current=1.34; %Current through Coils: Gives 10Gauss/cm Gradient
I = 0.8e3; %Intensity of Laser
for i=7:16
vcap = i;
[t,Y] = ode45(@(t,y) TrapOsc(t,y,I,Current), tspan, [x0; vcap]);
plot(t, Y(:,1))
ylim([-1e-3 3.5e-3])
hold on
end

xlabel('Time (s)')
ylabel('Displacement from Trap Center (m)') 

 Legend=cell(10,1);
 for iter=1:10
   x=iter+6;
   %Legend{iter}=strcat('Intensity(mW/cm^2)=', num2str(x));
   Legend{iter}=strcat(num2str(x), 'm/s');
 end
 title(legend,'Initial Velocity')
 legend(Legend)