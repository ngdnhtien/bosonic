
from scipy.optimize import curve_fit

def fit_function(x_values, y_values, function, init_params):
    """fit_function _summary_

    _extended_summary_

    Args:
        x_values (_type_): _description_
        y_values (_type_): _description_
        function (_type_): _description_
        init_params (_type_): _description_

    Returns:
        _type_: _description_
    """
    popt, pcov = curve_fit(function, x_values, y_values, init_params)
    y_fit = function(x_values, *popt)

    return popt, pcov, y_fit

def are_elements_unique(arr):
    for i in range(len(arr)):
        for j in range(i + 1, len(arr)):
            if arr[i] == arr[j]:
                print(i, j)
                return False
    return True