classdef GIRF_fft < handle
    
    properties
        gstf
        girf
        t_axis
        dt
        f_axis
        df
        dwelltime_compensated
    end % properties
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GIRF_fft(dt)
        %%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
            obj.dt = dt;
            obj.dwelltime_compensated = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_fft(obj, input, output, skipTriangles)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % input: [ROPs, triangles]
            % output: [channels, ROPs, triangles]
            % skipTriangles is just a 1D array
            
            obj.gstf = zeros(size(output,2),size(output,1));
            
            for channel=1:1:size(output,1)
                input_signal = input;
                output_signal = squeeze(output(channel,:,:));
                
                % Delete unwanted triangles
                for triang=size(input_signal,2):-1:1
                    if ismember(triang, skipTriangles)
                        input_signal(:,triang) = [];
                        output_signal(:,triang) = [];
                    end
                end % triang
                
                out_fft = fft1d(output_signal,1);
                in_fft = fft1d(input_signal,1);
                Norm = sum(abs(in_fft).^2,2);
                
                obj.gstf(:,channel) = sum((conj(in_fft)).*(out_fft),2)./Norm;
            end % channel
            obj.girf = real(fft1d(obj.gstf,1));
            
            % Calculate the time and frequency axes
            F = 1/obj.dt;
            obj.df = F/(size(obj.gstf,1)-1);
            obj.f_axis = (-F/2:obj.df:F/2);
            T = 1/obj.df;
%             obj.t_axis = (0:obj.dt:T);
            obj.t_axis = (0:obj.dt:T) + obj.dt/2; % 2023-10-09
        end % calcH_fft
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dwelltime_compensation(obj, dwelltime)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~obj.dwelltime_compensated
                obj.gstf = obj.gstf./sinc(dwelltime*obj.f_axis).';
                obj.girf = real(fft1d(obj.gstf,1));
                obj.dwelltime_compensated = 1;
            end
        end % dwelltime_compensation
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [outputMatrix_calc, residualMatrix] = forwardCalculation_fft(obj, input, outputMatrix, i_shift, f_axis_highRes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % input{i}: [timepoints, numTriang] (timepoints from start of TR until end of last ADC of the ith measurement)
            % outputMatrix: [lengthADC, numADC, numTriang, numMeasurements, channels]
            lengthADC = size(outputMatrix,1);
            numADC = size(outputMatrix,2);
            numTriang = size(outputMatrix,3);
            numMeas = size(outputMatrix,4);
            numChannels = size(outputMatrix,5);
            
            gstf_interp = interp1(obj.f_axis, obj.gstf, f_axis_highRes);
            lengthH_interp = size(f_axis_highRes,2);
            
            outputMatrix_calc = zeros(size(outputMatrix));
            for meas=1:1:numMeas
                input_meas = input{meas}; % [timepoints, triangles] (timepoints from start of TR until end of last ADC)
                if size(input_meas,1)<lengthH_interp
                    input_meas = vertcat(input_meas, zeros(lengthH_interp-size(input_meas,1),size(input_meas,2)));
                else
                    input_meas = input_meas(1:lengthH_interp,:);
                end
                input_meas = fft1d(input_meas,1); % [lengthH_interp, triangles]
                output_calc_meas = zeros(numChannels, size(input_meas,1), size(input_meas,2));
                for channel=1:1:numChannels
                    output_calc_meas(channel,:,:) = input_meas.*gstf_interp(:,channel);
                end
                clearvars input_meas;
                output_calc_meas = real(ifft1d(output_calc_meas,2)); % [channels, lengthH_interp, triangles]
                
                for channel=1:1:numChannels
                    for triang=1:1:numTriang
                        for adc=1:1:numADC
                            idx = round(i_shift{meas}(adc,triang)+1);
                            if idx+lengthADC > lengthH_interp
                                outputMatrix_calc(1:lengthH_interp-idx+1,adc,triang,meas,channel) = squeeze(output_calc_meas(channel,idx:end,triang));
                            else
                                outputMatrix_calc(:,adc,triang,meas,channel) = squeeze(output_calc_meas(channel,idx:idx+lengthADC-1,triang));
                            end
                        end
                    end
                end
                clearvars output_calc_meas;
            end % meas
            
            % Calculate residuals between measured and calculated output gradient signal
            residualMatrix = outputMatrix - outputMatrix_calc;
        end % forwardCalculation_fft
        
    end % methods
    
end % classdef







