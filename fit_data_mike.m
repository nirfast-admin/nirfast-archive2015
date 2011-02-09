function [mua,mus,lnI_offset,phase_offset,fem_data] = fit_data_mike(mesh,...
                                                               data,...
                                                               frequency,...
                                                               iteration)
                                                           
% [mua,mus,lnI_offset,phase_offset,fem_data] = fit_data(mesh,...
%                                   data,frequency,iteration)
%
% fits data to a given model to give initial guesses for
% reconstruction as well as data offsets
%
% mesh is the input mesh (variable)
% data is the boundary data (variable)
% frequency is the modulation frequency (MHz)
% iteration is the number of iterations for fitting
% outputs are the intial guesses

link_orig = mesh.link;
link_r = importdata('link_reflectance.txt');
mesh.link = link_r;
data_all = data;
a = link_orig==0;
mesh.link(a) = 0;
% calculate the source / detector distance for each measurement
k = 1; r = 1;
datanum = 0;
[ns,junk]=size(mesh.source.coord);
for i = 1 : ns
  for j = 1 : length(mesh.link(i,:))
      datanum = datanum + 1;
      % original fit data loop
      if isnan(data(datanum,1)) || isnan(data(datanum,2))
          mesh.link(i,j) = 0;
          link_orig(i,j) = 0; 
      end
      if mesh.link(i,j) ~= 0
          jj = mesh.link(i,j);
          dist(k,1) = sqrt(sum((mesh.source.coord(i,:) - ...
                    mesh.meas.coord(jj,:)).^2));
          k = k+1;
      else
        data(datanum,1) = NaN;
      end
      % added to find all data in link_original. needed for offset.
      if link_orig(i,j) == 0
          data_all(datanum,1) = NaN;
      end
  end
end

% % convert log amplitude into amplitude
% data(:,1) = exp(data(:,1));

% Set lnrI, lnr and phase!
[j,k] = size(data(:,1));
[j2,k2] = size(dist);

% deal with NaNs
dist_orig = dist;
ind = unique([find(isnan(data(:,1))==1); find(isnan(data(:,2))==1)]);
ind = setdiff(1:size(data,1),ind);
data = data(ind,:);

% find all data for calculating offset later.
ind_all = unique([find(isnan(data_all(:,1))==1); find(isnan(data_all(:,2))==1)]);
ind_all = setdiff(1:size(data_all,1),ind_all);
data_all = data_all(ind_all,:);

lnrI = log(data(:,1).*dist);
lnI = log(data(:,1));
phase = data(:,2);

figure;
subplot(2,2,1);
plot(dist,lnrI,'.')
ylabel('lnrI');
xlabel('Source / Detector distance');
subplot(2,2,2);
plot(dist,phase,'.')
ylabel('Phase');
xlabel('Source / Detector distance');
drawnow
pause(0.001)

% Calculate the coeff of a polynomial fit of distance vs. Phase or lnrI
m0 = polyfit(dist,phase,1); m0 = m0(1);
m1 = polyfit(dist,lnrI,1); m1 = m1(1);

% fit data using an analytical model
% based on Pogue paper
omega = 2*pi*frequency*1e6;
c=(3e11./mean(mesh.ri));

mua = 0.01; mus = 1; kappa = 1/(3*(mua+mus));

for i = 1 : 25
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  
  mua = mua + 0.0001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha1 = -f1*cos(0.5*atan2(omega/c,mua));
  phi1 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  mua = mua - 0.0001;   
  da = (alpha-alpha1)/0.0001;
  ds = (phi-phi1)/0.0001;
  mua = mua - (m1-alpha)/da;
  mus = mus - (m0-phi)/ds*0.1;
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  
  mus = mus + 0.001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha2 = -f1*cos(0.5*atan2(omega/c,mua));
  phi2 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  mus = mus - 0.001; 
  da = (alpha-alpha2)/0.001;
  ds = (phi-phi2)/0.001;
  mua = mua - (m1-alpha)/da*0.01;
  mus = mus - (m0-phi)/ds*0.2;
  kappa = 1/(3*(mua+mus));
end

disp('Global values calculated from Analytical fit');
disp(['Absorption = ' num2str(mua) ' mm-1']);
disp(['Scatter    = ' num2str(mus) ' mm-1']);
disp('============================================');

% Set the global values onto mesh.
mesh.mua(:) = mua;
mesh.mus(:) = mus;
mesh.kappa(:) = kappa;

% Fit for mua and mus using FEM
dist = dist_orig;
jj = 0;
while jj ~= iteration
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
    
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);
  
  phi0 = polyfit(dist,femphase,1); phi0 = phi0(1);
  alpha0 = polyfit(dist,femlnrI,1); alpha0 = alpha0(1);
  
  mesh.mua(:) = mesh.mua(:)+0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);
  
  phi1 = polyfit(dist,femphase,1); phi1 = phi1(1);
  alpha1 = polyfit(dist,femlnrI,1); alpha1 = alpha1(1);
  
  mesh.mua(:) = mesh.mua(:)-0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  da = (alpha0-alpha1)/0.0001;
  ds = (phi0-phi1)/0.0001;
  
  mesh.mua(:) = mesh.mua(:) - (m1-alpha0)/da;
  mesh.mus(:) = mesh.mus(:) - (m0-phi0)/ds.*0.1;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
 [mesh.mua(1) mesh.mus(1)];
 
 %
 
 [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);

  phi0 = polyfit(dist,femphase,1); phi0 = phi0(1);
  alpha0 = polyfit(dist,femlnrI,1); alpha0 = alpha0(1);
  
  mesh.mus(:) = mesh.mus(:)+0.001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI = log(fem_data(:,1).*dist);
  femphase = fem_data(:,2);
  
  phi2 = polyfit(dist,femphase,1); phi2 = phi2(1);
  alpha2 = polyfit(dist,femlnrI,1); alpha2 = alpha2(1);
  
  mesh.mus(:) = mesh.mus(:)-0.001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  da = (alpha0-alpha2)/0.001;
  ds = (phi0-phi2)/0.001;
  
  mesh.mua(:) = mesh.mua(:) - abs((m1-alpha0)/da*0.002);
  mesh.mua(find(mesh.mua(:)<0)) = mesh.mua(find(mesh.mua(:)<0)) + da;
  mesh.mus(:) = mesh.mus(:) - (m0-phi0)/ds.*0.2;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  err_a = abs(mean(mesh.mua)-mua)./mua;
  err_s = abs(mean(mesh.mus)-mus)./mus;
  
  if ((err_a < 0.001) & (err_s < 0.001))
    jj = iteration;
  else
    jj = jj + 1;
  end
  
  mua = mean(mesh.mua);
  mus = mean(mesh.mus);
  disp('Global values calculated from Numerical fit');
  disp(['Iteration ' num2str(jj) ' of ' num2str(iteration)]);
  disp(['Absorption = ' num2str(mua) ' mm-1 with error of ' num2str(err_a)]);
  disp(['Scatter    = ' num2str(mus) ' mm-1 with error of ' num2str(err_s)]);
  disp('-------------------------------------------------');
end
%%%%%% set back to original
mesh.link = link_orig;
% Calculate data based on these global values
[fem_data]=femdata(mesh,frequency);
fem_data = fem_data.paa;


% Arrange data to calculate offset
femlnI = log(fem_data(:,1));
femphase = fem_data(:,2);

% Find offset
% make sure any original NaNs are present
lnI = log(data_all(:,1));
phase = data_all(:,2);

% Set offset based on particular source / detector
% we do this because data is not symmetrical!
[n,m]=size(mesh.link);
spot = 1;
lnI_offset = zeros(size(lnI,1),1);
phase_offset = zeros(size(phase,1),1);
for i=1:n
    num_non = sum(mesh.link(i,:)~=0);
    lnI_offset(spot:spot+num_non-1) = mean(lnI(spot:spot+num_non-1)-femlnI(spot:spot+num_non-1));
    phase_offset(spot:spot+num_non-1) = mean(phase(spot:spot+num_non-1)-femphase(spot:spot+num_non-1));
    spot = spot + num_non;
end

subplot(2,2,3);
plot(lnI,'k');
hold on
plot(femlnI+lnI_offset,'r--');
axis tight;
xlabel('log Amplitude');
legend('original','Calibrated');
subplot(2,2,4);
plot(phase,'k');
hold on
plot(femphase+phase_offset,'r--');
axis tight;
xlabel('Phase');
legend('original','Calibrated');

% restore NaNs to data
fem_datatmp = fem_data;
lnI_offsettmp = lnI_offset;
phase_offsettmp = phase_offset;
fem_data = [];
lnI_offset = [];
phase_offset = [];
datanum = 0;
for i = 1 : n
  for j = 1 : m
      if mesh.link(i,j) ~= 0
          datanum = datanum + 1;
          fem_data = [fem_data; fem_datatmp(datanum,:)];
          lnI_offset = [lnI_offset; lnI_offsettmp(datanum)];
          phase_offset = [phase_offset; phase_offsettmp(datanum)];
      else
          fem_data = [fem_data; NaN NaN];
          lnI_offset = [lnI_offset; NaN];
          phase_offset = [phase_offset; NaN];
      end
  end
end
