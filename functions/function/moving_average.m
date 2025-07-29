function smoothed_signal = moving_average(signal, window_size, type)
    % Function to apply different types of moving averages to a multi-channel signal
    %
    % Parameters:
    % signal      - Input signal to be smoothed (samples, number of channels)
    % window_size - Size of the moving window
    % type        - Type of moving average ('simple', 'weighted', 'exponential')
    
    if nargin < 3
        type = 'simple'; % Default to simple moving average
    end
    
    [num_samples, num_channels] = size(signal);
    smoothed_signal = zeros(num_samples, num_channels);
    
    for channel = 1:num_channels
        switch type
            case 'simple'
                smoothed_signal(:, channel) = simple_moving_average(signal(:, channel), window_size);
                
            case 'weighted'
                smoothed_signal(:, channel) = weighted_moving_average(signal(:, channel), window_size);
                
            case 'exponential'
                smoothed_signal(:, channel) = exponential_moving_average(signal(:, channel), window_size);
                
            otherwise
                error('Unknown type. Use ''simple'', ''weighted'', or ''exponential''.');
        end
    end
end

function sma = simple_moving_average(signal, window_size)
    % Simple moving average
    sma = filter(ones(1, window_size)/window_size, 1, signal);
end

function wma = weighted_moving_average(signal, window_size)
    % Weighted moving average
    weights = (1:window_size) / sum(1:window_size);
    wma = filter(weights, 1, signal);
end

function ema = exponential_moving_average(signal, window_size)
    % Exponential moving average
    alpha = 2 / (window_size + 1);
    ema = filter(alpha, [1 alpha-1], signal, signal(1) * (1-alpha));
end
