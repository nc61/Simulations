function [xdata, signal] = center_data_to_peak_fcn(xdata, signal)
    [~, index_of_peak] = max(abs(signal));
    shift_value = xdata(index_of_peak);
    xdata = xdata - shift_value;
end

