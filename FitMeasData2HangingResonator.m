% Copyright Chris Wilen 2018
% 
function fitMeasData2HangingResonator(data_variable1)

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    data_variable1 = selectDepDataVars(data, true);
    if isempty(data_variable1)
        return
    end
    data_variable1 = data_variable1{1};
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Check that the data variable exists (compute it if necessary).
data_is_magphase = 0;
if strcmp(data_variable1, 'S21') || strcmp(data_variable1, 'S43')
    data_variable2 = [data_variable1,'_Phase'];
    data_is_magphase = 1;
elseif strcmp(data_variable1,'I')
    data_variable2 = 'Q'
elseif strcmp(data_variable1,'Q')
    data_variable1 = 'I'
    data_variable2 = 'Q'
else
    error(['Dependent variables need to be (S21, S21_Phase),', ...
           '(S43, S43_Phase), or (I,Q).'])
end
    
[data, data_variable1] = checkDataVar(data, data_variable1);
[data, data_variable2] = checkDataVar(data, data_variable2);

for data_variable = {data_variable1, data_variable2}
    data_var = data_variable{1};
    fitted_name = ['Fitted_', data_var];
    if strcmp(data.units.(data_var), 'dB')
        data.(data_var) = 10.^(data.(data_var) / 20);
        data.units.(data_var) = 'Arb. Units';
        disp('The values in dB are converted to arbitrary voltage units.')
    elseif strcmp(data.units.(data_var), 'V^2')
        data.(data_var) = sqrt(data.(data_var));
        data.units.(data_var) = 'V';
        disp(['The data values in V are square rooted to obtain a quantity',...
             ' proportional to the voltage.'])
    end
end

dep_vals1 = data.(data_variable1);
dep_rels1 = data.rels.(data_variable1)
dep_vals2 = data.(data_variable2);
dep_rels2 = data.rels.(data_variable2);


if isempty(dep_rels1)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable1, '_', ' '), ''' are not specified.'])
end

% Plot 1D data.  
if length(dep_rels1) == 1
    indep_name = dep_rels1{1};
    indep_vals = data.(indep_name);
    
    f = indep_vals;
    if data_is_magphase
        s_mag = dep_vals1;
        p = pi/180*dep_vals2;
    else
        s_mag = sqrt(dep_vals1.^2 + dep_vals2.^2);
        p = atan(dep_vals2./dep_vals1);
    end
    s_dB = 20*log10(s_mag);
    s_dB = s_dB - mean(s_dB);
    s_mag = 10.^(s_dB/20);
    p = p - mean(p);
    i = s_mag.*cos(p); 
    q = 1j*s_mag.*sin(p);
    s = (i+q);
    
    figure()
    
    subplot(2,2,1);
    plot(f, s_dB)
    title('S21 in dB');
    
    subplot(2,2,3);
    plot(f, 180/pi*p)
    title('S21 Phase');
    
    subplot(2,2,2);
    plot(i, -1j*q);
    axis equal;
    title('IQ')
    
    subplot(2,2,4);
    plot(real(1./s),imag(1./s));
    % axis equal;
    title('1/S21 IQ');
    
    
    opts = optimoptions(@lsqcurvefit,'Display','off','MaxIterations',10000,...
                        'MaxFunEvals',10000, 'TolX', 1e-9,'TolFun',1e-9, ...
                        'StepTolerance',1e-9);%, 'PlotFcn',@optimplotresnorm);
    
    objfcn = @(v,x)[real(v(10)*( 1 - (v(1)/v(2)-2j*v(1)*v(3)/v(4)) ./ (1+2j*v(1).*(x-v(4))/v(4)) ).* exp(1j*(2*pi*v(5)*(x-v(4))+v(6))) ) + 1*v(7), ...
                    imag(v(10)*( 1 - (v(1)/v(2)-2j*v(1)*v(3)/v(4)) ./ (1+2j*v(1).*(x-v(4))/v(4)) ).* exp(1j*(2*pi*v(5)*(x-v(4))+v(6))) ) + 1*v(8)];
    % Q0, Qc, df, f0, tau, phi, b_re, b_im, m, scale factor
     v0 = [350000; 450000; 0.05; 5.53446; -80e-9; 0; mean(i); mean(q); 0.5; 1];
    % v0 = [2000; 3000; 0.05; 7.5466; -80e-9; 0; mean(i); mean(q); 0.5; 1];
    [vestimated,resnorm] = lsqcurvefit(objfcn,v0,f.',[i.',-1j*q.'],[],[],opts)
    1/(1/vestimated(1) - 1/vestimated(2))
    temp = objfcn(vestimated,f.');
    ft = temp(:,1) + 1j*temp(:,2);
    
    subplot(2,2,1); hold on;
    plot(f, 20*log10(abs(ft)))
    subplot(2,2,3); hold on;
    plot(f, 180/pi*atan(imag(ft)./real(ft)))
    subplot(2,2,2); hold on;
    plot(real(ft),imag(ft))
    subplot(2,2,4); hold on;
    plot(real(1./ft),imag(1./ft))

elseif length(dep_rels1) == 2 % Plot 2D data.
    error(['2D data not yet implemented.'])

end