tspan = 0:1e-4:1e-2;
clf
vcap=7; %Capture Velocity
x0 = 0; %Initial Position
Current=0.67; %Current through Coils: Gives 5Gauss/cm Gradient
for i=0.5e3:0.25e3:2.5e3
I = i; %Intensity of Laser
[t,Y] = ode45(@(t,y) TrapOsc(t,y,I,Current), tspan, [x0; vcap]);
plot(t, Y(:,1))
ylim([-1e-3 3.5e-3])
hold on
end

xlabel('Time (s)')
ylabel('Displacement from Trap Center (m)') 

 Legend=cell(9,1);
 for iter=1:9
   x=25*iter+25;
   %Legend{iter}=strcat('Intensity(mW/cm^2)=', num2str(x));
   Legend{iter}=strcat(num2str(x), 'mW/cm^2');
 end
 title(legend,'Laser Intensity')
 legend(Legend)