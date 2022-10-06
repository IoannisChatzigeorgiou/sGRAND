% --------------------------------------------------------------------
% This MATLAB script simulates a system that consists of a source, 
% and a destination. The source uses a k/n random linear code (RLC) 
% and a 16-QAM modulator. The modulated symbols are transmitted over
% a Rayleigh fading channel. The destination attempts to recover the
% the sequence of k bits using symbol-level GRAND.
%
% The system is described in 
% "Symbol-Level GRAND for High-Order Modulation over Flat Fading Channels",
% which has been uploaded to arXiv (https://arxiv.org/abs/2207.07748)
% 
% Authors: 
% Ioannis Chatzigeorgiou, Lancaster University, United Kingdom
% Francisco Monteiro, Instituto de Telecommunicacoes 
% and ISCTE-Instituto Universitario de Lisboa, Portugal
%
% If you use this code or the simulation results, 
% please cite the paper above.
%
% If you modify this code, please cite the paper above and 
% share your code with the research community.
% --------------------------------------------------------------------

% INITIALISE PARAMETERS 
% Note: This is the Eb/N0 in dB, it is not the signal-to-noise ratio (SNR).
% The SNR is defined later.
Eb_N0dB = 40:-1:10;
realisations = [4e6 ...
                2e6*ones(1,5) ...
                6e5*ones(1,5) ...
                4e5*ones(1,5) ...
                2e5*ones(1,5) ...
                7e4*ones(1,5) ...
                5e4*ones(1,5)];                  
            
% Length of input word and output codeword
k = 103;
n = 128;

% Variables for symbol-level GRAND
% File name of lookup table
sim_name = ['lookup_k',num2str(k),'n',num2str(n),'M16.mat'];
% Maximum weight of error patterns considered,
% error structures that generate error patterns with weight higher
% than max_weight will not be considered.
max_weight = 3;
% Number of most likely error structures to be considered 
% (all of them will have weight less than or equal to max_weight)  
top_rows   = 5;

% -- END OF INITIALISATION --

% Obtain truncated lookup table based on the full lookup table
% and parameters max_weight and top_rows
load(sim_name, 'Eb_N0', 'analytical_ordered_error_types')

% Load the full lookup table, not just the entries that
% correspond to the values of Eb_N0dB
trunc_analytical_ordered_error_types = zeros(top_rows, 3, length(Eb_N0));
% Store the number of error types for each value of Eb/N0
% because they might be fewer than top_rows
trun_error_types_per_entry           = zeros(length(Eb_N0));
% Store the Eb/N0 step
lookup_size                          = length(Eb_N0);
lookup_Eb_N0_step                    = Eb_N0(2) - Eb_N0(1);

for Eb_N0_idx = 1:length(Eb_N0)

    temp = analytical_ordered_error_types(:, :, Eb_N0_idx);
    rows_recorded = 0;
    
    for i=1:size(temp, 1)
       
        % Check the Hamming weight of each error structure
        if (temp(i, 1) > 0) && (temp(i, 2) + 2*temp(i, 3) <= max_weight) && (rows_recorded < top_rows)
            rows_recorded = rows_recorded + 1;
            trunc_analytical_ordered_error_types(rows_recorded, :, Eb_N0_idx) = temp(i, :);            
        end
        
    end    
    
    if rows_recorded == 0
        lookup_size = Eb_N0_idx-1;
        Eb_N0 = Eb_N0(1:lookup_size);
        trun_error_types_per_entry = trun_error_types_per_entry(1:lookup_size);
        trunc_analytical_ordered_error_types = trunc_analytical_ordered_error_types(1:top_rows, 1:3, 1:lookup_size);
        break        
    end
    
    trun_error_types_per_entry(Eb_N0_idx) = rows_recorded;
    
end

clear Eb_N0 analytical_ordered_error_types

code_rate                     = k / n;
SNR                           = 4.*code_rate.*10.^(Eb_N0dB./10);
stddev_awgn_per_dimension     = sqrt(2./SNR);
stddev_Rayleigh_per_dimension = sqrt(1/2);

% Gray-coded 16-QAM constellation used:
%
% 0111  0011 | 1011  1111
% 0110  0010 | 1010  1110
% -----------------------
% 0100  0000 | 1000  1100
% 0101  0001 | 1001  1101

QAM_array = build_16QAM;
[adj_error_patterns, num_of_adj_err_patterns] = build_lookup_tables(QAM_array);

for Eb_N0_idx = 1:length(Eb_N0dB)

    disp(['Eb_N0 = ', num2str(Eb_N0dB(Eb_N0_idx)), ' dB']);   
        
    qam16_BLER          = 0;
    avg_patterns_tested = 0;   
    
    for rlz_idx = 1:realisations(Eb_N0_idx)

        if mod(rlz_idx, 100)==0
            disp(['-- Realisation: ', num2str(rlz_idx)])
        end
        
        % Random code generation (a new code for each iteration)
        P = randi([0 1], k, n-k);
        G = [eye(k) P];    
        H = [P' eye(n-k)];  
        % G*H' = 0 (mod 2)
           
        % Random input sequence
        u = randi([0 1], 1, k);        
       
        % Output codeword of RLC
        x = mod(u*G, 2);
        
        % Group bits into 4-bit words
        x_strings = reshape(x, 4, (1/4)*n);

        % 16-QAM (GRAY) mapping (lame but speedy implementation)
        % Ensures that if Eb=1 (energy per bit),
        % then Es=4 (energy per symbol), given that 1 symbol carries 4 bits
        s = zeros(1,(1/4)*n);
        s(x_strings(1,:)==1 & x_strings(2,:)==1 & x_strings(3,:)==1 & x_strings(4,:)==1) =  (3+1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==1 & x_strings(3,:)==1 & x_strings(4,:)==0) =  (3+1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==1 & x_strings(3,:)==0 & x_strings(4,:)==0) =  (3-1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==1 & x_strings(3,:)==0 & x_strings(4,:)==1) =  (3-1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==0 & x_strings(3,:)==1 & x_strings(4,:)==1) =  (1+1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==0 & x_strings(3,:)==1 & x_strings(4,:)==0) =  (1+1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==0 & x_strings(3,:)==0 & x_strings(4,:)==0) =  (1-1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==1 & x_strings(2,:)==0 & x_strings(3,:)==0 & x_strings(4,:)==1) =  (1-1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==0 & x_strings(3,:)==1 & x_strings(4,:)==1) = (-1+1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==0 & x_strings(3,:)==1 & x_strings(4,:)==0) = (-1+1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==0 & x_strings(3,:)==0 & x_strings(4,:)==0) = (-1-1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==0 & x_strings(3,:)==0 & x_strings(4,:)==1) = (-1-1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==1 & x_strings(3,:)==1 & x_strings(4,:)==1) = (-3+1i*3)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==1 & x_strings(3,:)==1 & x_strings(4,:)==0) = (-3+1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==1 & x_strings(3,:)==0 & x_strings(4,:)==0) = (-3-1i*1)/sqrt(2.5) ;
        s(x_strings(1,:)==0 & x_strings(2,:)==1 & x_strings(3,:)==0 & x_strings(4,:)==1) = (-3-1i*3)/sqrt(2.5) ;

        % Additive White Gaussian Noise
        % Both the REAL and the IMAG part of the transmitted_symbols should
        % experience the additive noise
        z = stddev_awgn_per_dimension(Eb_N0_idx)*randn(1,(1/4)*n) + 1i*stddev_awgn_per_dimension(Eb_N0_idx).*randn(1,(1/4)*n);                
        % Quasi-static Rayleigh fading
        h = stddev_Rayleigh_per_dimension*randn(1,1) + 1i*stddev_Rayleigh_per_dimension*randn(1,1);
        r = h .* s + z;
        % Coherent detection
        r = r .* (conj(h)/abs(h)^2);             
                        
        % Hard decoding of received symbols
        y_strings = zeros(4,(1/4)*n);       
        y_strings(1,real(r)>=0) = 1 ;
        y_strings(1,real(r)<0)  = 0 ;
        y_strings(2,(abs(real(r))-2/sqrt(2.5))>=0) = 1 ;
        y_strings(2,(abs(real(r))-2/sqrt(2.5))<0)  = 0 ;
        y_strings(3,imag(r)>=0) = 1 ;
        y_strings(3,imag(r)<0)  = 0 ;
        y_strings(4,(abs(imag(r))-2/sqrt(2.5))>=0) = 1 ;
        y_strings(4,(abs(imag(r))-2/sqrt(2.5))<0)  = 0 ;
               
        % Sequence of hard-detected bits
        y = reshape(y_strings,1,n);
         
        % Find the appropriate error type / structure for the
        % current value of the fading coefficient
        Eb_N0_lut_idx = 1 + round((Eb_N0dB(Eb_N0_idx) + 20*log10(abs(h)))* 1/lookup_Eb_N0_step);        
        if Eb_N0_lut_idx < 1, Eb_N0_lut_idx = 1; end
        if Eb_N0_lut_idx > lookup_size, Eb_N0_lut_idx = lookup_size; end     
        fading_specific_ordered_error_types = ...
            trunc_analytical_ordered_error_types(1:trun_error_types_per_entry(Eb_N0_lut_idx), :, Eb_N0_lut_idx);      

        % Symbol-level GRAND
        [est_e, tests, stopped] = ...
            symbol_level_GRAND(y_strings, H, ...
            adj_error_patterns, ...
            num_of_adj_err_patterns, ...
            fading_specific_ordered_error_types);
                
        if stopped == 0
            % Estimated codeword
            est_x = mod(y + est_e, 2);  
            
            % Error measurement
            errors = (est_x ~= x);  
            
        else         
            errors = (y ~= x); 
        end
                     
        % Block error rate measurement   
        if sum(errors)>0
            qam16_BLER = qam16_BLER + 1;
        end
        
        % Complexity measurement
        avg_patterns_tested = avg_patterns_tested + tests;
        
    end
    
    qam16_BLER          = qam16_BLER / realisations(Eb_N0_idx);
    avg_patterns_tested = avg_patterns_tested / realisations(Eb_N0_idx);
    
    % The values of max_weight and top_rows are included in the name
    % of the data file
    fname = ['BLER_Rayleigh_RLC_16QAM_symbol_GRAND_mw',num2str(max_weight),'tr',num2str(top_rows),'_', num2str(Eb_N0dB(Eb_N0_idx)), 'dB.mat'];
    save(fname, 'qam16_BLER', 'avg_patterns_tested', 'k', 'n')
    
end


% ---------------- Definition of called functions ----------------

% --- FUNCTION build_16QAM

function get_QAM_array = build_16QAM
% The function returns a look-up table get_QAM_array{bit}{x}{y}.
% get_QAM_array(:, x, y) holds the QAM symbol in the x-th column (x axis) 
% and y-th row (y axis). For example, the positions of the symbols for 
% 16-QAM will be:
% (:,1,1) (:,2,1) (:,3,1) (:,4,1)
% (:,1,2) (:,2,2) (:,3,2) (:,4,2)
% (:,1,3) (:,2,3) (:,3,3) (:,4,3)
% (:,1,4) (:,2,4) (:,3,4) (:,4,4)
%
% THIS FUNCTION SHOULD BE CALLED ONLY ONCE, WHEN THE CARDINALITY OF THE
% MODULATION SCHEME IS SELECTED. DO NOT INCLUDE IT IN A LOOP!

cardinality           = 16;
number_of_PAM_symbols = sqrt(cardinality);
bits_per_PAM_symbol   = log2(cardinality)/2;
PAM_vector            = [1 1; 1 0; 0 0; 0 1];
get_QAM_array         = zeros(2*bits_per_PAM_symbol, number_of_PAM_symbols, number_of_PAM_symbols);

% Gray-coded M-QAM constellation
for col=1:number_of_PAM_symbols
    for row=1:number_of_PAM_symbols
        get_QAM_array(1:bits_per_PAM_symbol,col,row) = PAM_vector(number_of_PAM_symbols-col+1, :);
        get_QAM_array(bits_per_PAM_symbol+1:2*bits_per_PAM_symbol,col,row) = PAM_vector(row, :);
    end
end

% % Display constellation (use only when debugging)
% for row=1:number_of_PAM_symbols
%     constellation_row = [];
%     for col=1:number_of_PAM_symbols
%         constellation_row = [constellation_row, char(get_QAM_array(:,col,row)' + '0'), ' ']; 
%     end
%     disp(constellation_row)
% end

end % FUNCTION build_16QAM


% --- FUNCTION build_lookup_tables

function [get_adj_errors, get_num_of_adj_errors] = build_lookup_tables(input_QAM_array)
% The function accepts one input:
% - input_QAM_array: the array obtained from function build_square_QAM.
%
% The function returns two look-up tables:
% - get_adj_errors{i}{j}{k}: this is the k-th error pattern of the symbol 
% that is adjacent to the i-th received symbol and can be found in the 
% j-th neighbourhood around the received symbol.
% - get_num_of_adj_errors(i, j): this is size of the j-th neighbourhood 
% of error patterns around the i-th received symbol.
%
% THIS FUNCTION SHOULD BE CALLED ONLY ONCE, WHEN THE CARDINALITY OF THE
% MODULATION SCHEME IS SELECTED. DO NOT INCLUDE IT IN A LOOP!

cardinality           = size(input_QAM_array, 1)*size(input_QAM_array, 2);
number_of_PAM_symbols = sqrt(cardinality);
bits_per_PAM_symbol   = log2(cardinality)/2;

% Consider three different neighbourhoods around the received symbol.
% All of of them contain four members or less, depending on the position 
% of the received symbol in the M-QAM constellation)
neighbourhoods{1} = {[0 1], [1 0],  [0 -1],  [-1 0]}; % Most likely
neighbourhoods{2} = {[1 1], [-1 1], [-1 -1], [1 -1]}; % Second most likely

get_num_of_adj_errors = zeros(cardinality, 2);

for i = 1:number_of_PAM_symbols
    
    for j = 1:number_of_PAM_symbols
        
        % convert the value of the binary QAM symbol into an integer
        % and ADD 1 (to avoid having value 0 as an index)
        idx = 0;
        for bit = 1:2*bits_per_PAM_symbol
            idx = idx + input_QAM_array(bit,i,j)*(2^(2*bits_per_PAM_symbol-bit));
        end
        idx = idx + 1;
        
        % We consider two neighbourhoods
        for hood = 1:2  
            
            % Count how many neighbours (adjacent symbols)
            % in the neighbourhood under consideration
            counter = 0; 

            % Each neighbourhood should contain up to four members
            for member = 1:4 
             
                new_i = i + neighbourhoods{hood}{member}(1);
                new_j = j + neighbourhoods{hood}{member}(2);
                
                if (new_i >= 1) && (new_i <= number_of_PAM_symbols) && (new_j >= 1) && (new_j <= number_of_PAM_symbols)

                    counter = counter + 1;
                    
                    % EITHER activate the line below to store the symbols
                    % that are adjacent to the received symbol
                    % adj_symbols{idx}{hood}{counter} = input_QAM_array(:,new_i,new_j);
                    
                    % OR ativate the line below to store the error patterns
                    % of the symbols that are adjacent to the received symbol
                    get_adj_errors{idx}{hood}{counter} = mod(input_QAM_array(:,i,j)+input_QAM_array(:,new_i,new_j), 2);         
                    
               end
                               
            end % member
            
            get_num_of_adj_errors(idx, hood) = counter;
            
        end % hood
        
    end % j
    
end % i

end


% --- FUNCTION symbol_level_GRAND

function [estimated_noise, num_of_queries, terminated] = symbol_level_GRAND(y_rec_bin_symbs_matrix, parity_check_matrix, adjacent_patterns, num_of_patterns, seq_error_types)

% Number of error types to be considered in the algorithm.
% Change the size of the input argument 'seq_error_types', i.e., 
% increase or decrease the number of rows in 'seq_error_types', 
% to increase or decrease the number of error types.
max_num_error_types = size(seq_error_types, 1);

found      = 0;
terminated = 0;

bits_in_QAM_symbol       = size(y_rec_bin_symbs_matrix, 1); % Bits carried by each QAM symbol
received_QAM_symbols     = size(y_rec_bin_symbs_matrix, 2); % Number of QAM symbols
codeword_length_in_bits  = received_QAM_symbols * bits_in_QAM_symbol; % Length of codeword
estimated_noise          = zeros(codeword_length_in_bits, 1);
y_rec_dec_symbs_matrix   = bin2dec(char(y_rec_bin_symbs_matrix' + '0'))'; % QAM symbols in decimal format
y_decoded_bits           = reshape(y_rec_bin_symbs_matrix, 1, codeword_length_in_bits); % change from binary symbols to binary sequence
syndrome                 = mod(y_decoded_bits * parity_check_matrix', 2);
% Syndrome calculation:
% s = y*H' => s = (x+e)*H' => s = u*G*H' + e*H'=> s = e*H'

if sum(syndrome)==0 
% The received codeword appears to be a correct codeword,
% no need to search for an error pattern.
    found=1;
end

current_error_type = 1;
num_of_queries     = 1; % The first query was for the all-zero error pattern

%Testing all the error patterns with a number of errors 1, 2, 3,... up to
%the defined maximum number of errors :
while (~found) && (current_error_type <= max_num_error_types)
    
    outer_test = 1;
    
    % Compute all combinations of errors (the number of which is stored 
    % in the first column of the appropriate row of 'seq_error_types')
    % in a sequence of length 'received_QAM_symbols'.
    symb_errors_in_hood    = zeros(2,1);
    symb_errors_in_hood(1) = seq_error_types(current_error_type, 2);
    symb_errors_in_hood(2) = seq_error_types(current_error_type, 3);   
    total_symb_errors      = seq_error_types(current_error_type, 1);
        
    outer_positions = nchoosek(1:received_QAM_symbols, total_symb_errors);
    inner_positions = nchoosek(1:total_symb_errors, symb_errors_in_hood(1));
    
    if (symb_errors_in_hood(1) == 0) || (symb_errors_in_hood(2) == 0)
        num_inner_patterns = 1;
    else
        num_inner_patterns = size(inner_positions, 1);
    end

    % Store all possible combinations for
    % a given value of 'total_symb_errors'.
    num_outer_patterns = size(outer_positions,1);    
        
    while (~found) && (outer_test <= num_outer_patterns)
        
        inner_test = 1;
        
        % Identify the received symbols that their positions have been 
        % flagged as containing errors
        symbols_in_error = y_rec_dec_symbs_matrix(outer_positions(outer_test, :));   

        while (~found) && (inner_test <= num_inner_patterns)
        
            patterns_per_symbol = zeros(1, total_symb_errors);
                
            % Go over the indices that correspond to symbols that 
            % will transition to neighbourhood #1
            if symb_errors_in_hood(1) > 0
                
                for symb_idx = inner_positions(inner_test, :)
                    
                    lookup_idx = symbols_in_error(symb_idx) + 1;
                    
                    for err_pattern=1:num_of_patterns(lookup_idx, 1)
                        store_adj_errors{symb_idx}{err_pattern} = ...
                            adjacent_patterns{lookup_idx}{1}{err_pattern}';
                    end
                    
                    patterns_per_symbol(symb_idx) = ...
                        patterns_per_symbol(symb_idx) + num_of_patterns(lookup_idx, 1);
                    
                end
                
            end
        
            % Go over the remaining indices that correspond to symbols that 
            % will transition to neighbourhood #2                      
            if symb_errors_in_hood(2) > 0
                
                if symb_errors_in_hood(1) > 0
                    % 'inner_positions' gives the positions of symbols for which
                    % neighbourhood #1 will be considered.
                    % 'inner_positions_complement' finds the positions of the
                    % remaining symbols for which neighbourhood #2 will be
                    % considered.
                    inner_positions_complement = setdiff(1:total_symb_errors, inner_positions(inner_test, :));
                else
                    inner_positions_complement = 1:total_symb_errors;
                end
                                
                for symb_idx = inner_positions_complement
                    
                    lookup_idx = symbols_in_error(symb_idx) + 1;
                    
                    for err_pattern=1:num_of_patterns(lookup_idx, 2)
                        store_adj_errors{symb_idx}{err_pattern} = ...
                            adjacent_patterns{lookup_idx}{2}{err_pattern}';
                    end
                    
                    patterns_per_symbol(symb_idx) = ...
                        patterns_per_symbol(symb_idx) + num_of_patterns(lookup_idx, 2);
                    
                end
                
            end            
            
            % Try out all bit error patterns in the consdired neighbourhood
            % of the received symbols for a specific symbol error pattern            
            tested_patterns_per_symbol = ones(1, total_symb_errors);
            % The i-th element in tested_patterns_per_symbol takes values
            % in the range 1, ..., patterns_per_symbol(i).

            point_to_position  = total_symb_errors;

            y_QAM_error_matrix = zeros(received_QAM_symbols, bits_in_QAM_symbol);
           
            while (~found) && (tested_patterns_per_symbol(1) <= patterns_per_symbol(1))
                        
                for symb_idx = 1:total_symb_errors                    
                    
                    y_QAM_error_matrix(outer_positions(outer_test, symb_idx), :) = ...
                        store_adj_errors{symb_idx}{tested_patterns_per_symbol(symb_idx)};
                
                end

                estimated_noise = reshape(y_QAM_error_matrix', 1, codeword_length_in_bits);                      
                est_syndrome    = mod(estimated_noise * parity_check_matrix', 2);
                num_of_queries  = num_of_queries + 1;
                                                                        
                if est_syndrome == syndrome
                
                    found = 1; 
                                
                else
                
                    while (tested_patterns_per_symbol(point_to_position) + 1 > patterns_per_symbol(point_to_position)) && (point_to_position > 1)
                        tested_patterns_per_symbol(point_to_position) = 1;
                        point_to_position = point_to_position - 1;
                    end
                
                    tested_patterns_per_symbol(point_to_position) = tested_patterns_per_symbol(point_to_position) + 1;
                    point_to_position = total_symb_errors;
                    
                end
                   
            end
                
            inner_test = inner_test + 1;                        
                    
        end
        
        outer_test = outer_test + 1;
        
    end
    
    current_error_type = current_error_type + 1;
    
end

%In case guessing up to the defined maximum number of errors did not find
%any valid codeword:
if (current_error_type == max_num_error_types) && (~found)
    terminated = 1;
end

end
