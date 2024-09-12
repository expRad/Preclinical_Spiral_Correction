classdef GIRF_matrix < handle

    % Copyright (c) 2024 Hannah Scholten
    
    properties
        girf
        gstf
        fieldOffsets
        t_axis
        dt
        f_axis
        df
        dwelltime_compensated
        lengthH
    end % properties
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GIRF_matrix(dt, lengthH)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
            obj.dt = dt;
            obj.dwelltime_compensated = 0;
            obj.lengthH = lengthH;
            % Calculate the time and frequency axes
            F = 1/obj.dt;
            obj.df = F/(lengthH-1);
            obj.f_axis = (-F/2:obj.df:F/2);
            T = 1/obj.df;
%             obj.t_axis = (0:dt:T);
            obj.t_axis = (0:obj.dt:T) + obj.dt/2; % 2023-10-09
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_matrix(obj, inputMatrix, outputMatrix)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % inputMatrix: [timepoints, lengthH + 1]
            % output: [timepoints, channels]
            
            H = sparse(inputMatrix)\outputMatrix; % [lengthH + 1, channels]
            H(end,:)
            obj.girf = H(1:end-1,:); % GIRF is contained in all but the last row
            obj.fieldOffsets = H(end,:); % last row contains the field offsets
            
            obj.gstf = fft1d(obj.girf,1); % [lengthH, channels]
        end % calcH_matrix
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_matrix_Tikhonov(obj, inputMatrix, outputMatrix, lambda, alpha)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % inputMatrix: [timepoints, lengthH + 1]
            % output: [timepoints, channels]
            obj.dwelltime_compensated = 0;
            
            % Prepare regularization matrix
            L = eye(size(inputMatrix,2)); % identity matrix
            L(1:obj.lengthH,1:obj.lengthH) = exp(alpha*obj.t_axis.') .* L(1:obj.lengthH,1:obj.lengthH); % introduce an exponential weighting into the regularization
%             figure;
%             plot(L.');
            b = zeros(size(inputMatrix,2), size(outputMatrix,2));
            % Using the identity matrix for the Tikhonov regularization favors the solution with the smalles norm.

            H = sparse([inputMatrix; lambda*L])\[outputMatrix; b]; % [lengthH + 1, channels]
            % https://www.math.uni-frankfurt.de/~harrach/talks/2014Bangalore_Lecture2.pdf
            
            obj.girf = H(1:end-1,:); % GIRF is contained in all but the last row
            obj.fieldOffsets = H(end,:); % last row contains the field offsets
            
            obj.gstf = fft1d(obj.girf,1); % [lengthH, channels]
            
        end % calcH_matrix_Tikhonov
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_matrix_Tikhonov_freqWeight(obj, inputMatrix, outputMatrix, lambda, alpha_array)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % inputMatrix: [timepoints, lengthH + 1]
            % output: [timepoints, channels]
            obj.dwelltime_compensated = 0;
            
            % Prepare regularization matrix
            L = eye(size(inputMatrix,2)); % identity matrix, size=lengthH+1
            FourierMatrix = zeros(obj.lengthH);
            for column=1:1:obj.lengthH
                omega = 2*pi*obj.f_axis(column);
                FourierMatrix(:,column) = squeeze(exp(1i*omega*obj.t_axis));
            end % column
            disp('        FourierMatrix created.')
%             FourierWithWeighting = (exp(alpha*obj.t_axis.').*eye(obj.lengthH)) * FourierMatrix;
            FourierWithWeighting = FourierMatrix;
            for row=1:1:obj.lengthH
                weights = exp(alpha_array*obj.t_axis(row));
%                 weights(weights>5e8) = 5e8;
                FourierWithWeighting(row,:) = FourierWithWeighting(row,:).*weights;
            end
            
%             figure;
%             plot(real(FourierWithWeighting(:,1)));
%             hold on;
%             plot(real(FourierWithWeighting(:,ceil(end/2))));
%             plot(real(FourierWithWeighting(:,floor(end/2))));
            
            L(1:obj.lengthH,1:obj.lengthH) = (FourierWithWeighting/FourierMatrix);
            clearvars FourierMatrix FourierWithWeighting;
%             find(isnan(L)==1)
%             figure;
%             plot(abs(L.'));
            disp('        RegularizationMatrix created.')
            disp(['issparse(L) = ',num2str(issparse(L))])
            disp(['max(L) = ',num2str(max(L,[],'all'))]);
            b = zeros(size(inputMatrix,2), size(outputMatrix,2));

%             H = real([inputMatrix; lambda*L]\[outputMatrix; b]); % [lengthH + 1, channels]
            H = real(pinv([inputMatrix; lambda*L])*[outputMatrix; b]); % [lengthH + 1, channels]
            % https://www.math.uni-frankfurt.de/~harrach/talks/2014Bangalore_Lecture2.pdf
            
            obj.girf = H(1:end-1,:); % GIRF is contained in all but the last row
            obj.fieldOffsets = H(end,:); % last row contains the field offsets
            
            obj.gstf = fft1d(obj.girf,1); % [lengthH, channels]
            
        end % calcH_matrix_Tikhonov_freqWeight
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dwelltime_compensation(obj, dwelltime)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~obj.dwelltime_compensated
                obj.gstf = obj.gstf./sinc(dwelltime*obj.f_axis).';
                gstf2 = fft1d(obj.girf,1);
                gstf2 = gstf2./sinc(dwelltime*obj.f_axis).';
                obj.girf = real(ifft1d(gstf2,1));
                obj.dwelltime_compensated = 1;
            end
        end % dwelltime_compensation
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function correct_GSTF_phase(obj, phaseAtZero, corr_delay)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            linear_array = (1:1:size(obj.gstf,1))*pi;
            corrected_phase = angle(obj.gstf.*exp(1i*(obj.f_axis*2*pi*corr_delay+linear_array)).');
            corrected_phase_atZero = corrected_phase(floor(size(corrected_phase,1)/2),2);
            disp(['phaseDifferenceAtZero = ',num2str(abs(corrected_phase_atZero - phaseAtZero))])
            if abs(corrected_phase_atZero - phaseAtZero) > 1
                linear_array = (0:1:size(obj.gstf,1)-1)*pi;
                disp('Changed linear_array.')
            end
            obj.gstf = obj.gstf.*exp(1i*(obj.f_axis*2*pi*corr_delay+linear_array)).';
        end % correct_GSTF_phase
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [outputMatrix_calc, residualMatrix] = forwardCalculation_matrix(obj, inputMatrix, outputMatrix)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % inputMatrix: [measurement points, lengthH]
            % outputMatrix: [measurement points, channels]
            
            % Perform forward calculation with the input matrix and the calculated GIRF
            outputMatrix_calc = inputMatrix(:,1:end-1) * obj.girf + obj.fieldOffsets;
            % Calculate residuals between measured and calculated output gradient signal
            residualMatrix = outputMatrix - outputMatrix_calc;
        end % forwardCalculation_matrix
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [outputMatrix_calc, residualMatrix] = forwardCalculation_fft(obj, input, outputMatrix, i_shift)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Dimensions:
            % input{i}: [timepoints, numTriang] (timepoints from start of TR until end of last ADC of the ith measurement)
            % outputMatrix: [lengthADC, numADC, numTriang, numMeasurements, channels]
            lengthADC = size(outputMatrix,1);
            numADC = size(outputMatrix,2);
            numTriang = size(outputMatrix,3);
            numMeas = size(outputMatrix,4);
            numChannels = size(outputMatrix,5);
            
            outputMatrix_calc = zeros(size(outputMatrix));
            for meas=1:1:numMeas
                input_meas = input{meas}; % [timepoints, triangles] (timepoints from start of TR until end of last ADC)
                if size(input_meas,1)<obj.lengthH
                    input_meas = vertcat(input_meas, zeros(obj.lengthH-size(input_meas,1),size(input_meas,2)));
                else
                    input_meas = input_meas(1:obj.lengthH,:);
                end
                input_meas = fft1d(input_meas,1); % [lengthH, triangles]
                output_calc_meas = zeros(numChannels, size(input_meas,1), size(input_meas,2));
                for channel=1:1:numChannels
                    output_calc_meas(channel,:,:) = input_meas.*obj.gstf(:,channel);
                end
                clearvars input_meas;
                output_calc_meas = real(ifft1d(output_calc_meas,2)); % [channels, lengthH, triangles]
                
                for channel=1:1:numChannels
                    for triang=1:1:numTriang
                        for adc=1:1:numADC
                            idx = round(i_shift{meas}(adc,triang)+1);
                            if idx+lengthADC > obj.lengthH
                                outputMatrix_calc(1:obj.lengthH-idx+1,adc,triang,meas,channel) = squeeze(output_calc_meas(channel,idx:end,triang))  + obj.fieldOffsets(channel);
                            else
                                outputMatrix_calc(:,adc,triang,meas,channel) = squeeze(output_calc_meas(channel,idx:idx+lengthADC-1,triang)) + obj.fieldOffsets(channel);
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








