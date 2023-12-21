from bigfish.stack import check_parameter
from ..utils import _get_varname
from CustomPandasFramework.computer_interface import get_datetime

def print_parameters(path_out, *parameters, printDateTime= True):
    print("Warning : depreciation. This get_datetime function should be called from CustomPandasFramework")

    """
    Print parameters into a .txt file
    
    Parameters
    ----------
        path_out : str
            full_path to saving location.
        *parameters : args(int, float, str)
            parameters to print into text file.
    """

    check_parameter(path_out = (str))
    if path_out[len(path_out)-1] == '/' : path_out += 'parameters.txt'
    elif path_out[len(path_out)-4 : len(path_out)] != '.txt' : path_out += '.txt'

    parameter_file = open(path_out, "w")

    #Header
    parameter_file.write("PARAMETERS\n")
    parameter_file.write("\n############\n")

    if printDateTime:
        datetime = get_datetime()
        parameter_file.write(datetime)
        parameter_file.write("############\n\n")
    
    lines= []
    for parameter in parameters :
        name = _get_varname(parameter)
        lines += ["{0} : {1}\n".format(name, parameter)]
    parameter_file.writelines(lines)

    parameter_file.close()


def print_dict(dic: dict, path_out):
    print("Warning : depreciation. This get_datetime function should be called from CustomPandasFramework")

    check_parameter(path_out = (str), dic = (dict))
    if path_out[len(path_out) -1] == '/' : path_out += 'dic.txt'
    elif path_out[len(path_out)-4 : len(path_out)] != '.txt' : path_out += '.txt'

    lines = []
    for key, value in dic.items():
        lines += "{0} : {1}\n".format(key, value)


    dict_file = open(path_out, "w")
    dict_file.writelines(lines)
    dict_file.close()

def dict_to_lines(di : dict) :
    print("Warning : depreciation. This get_datetime function should be called from CustomPandasFramework")

    for key, value in di.items() : 
        yield "{0} : {1}\n".format(key, value)