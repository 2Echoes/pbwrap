from .getdata import get_datetime, _get_varname
import pandas as pd

class log :
    def __init__(self, name, path) -> None:
        self.name = name
        self.path = path
        self.creationdate = get_datetime()
        if not self.path.endswith('/') : self.path += '/'
        if not self.name.endswith('.txt') : self.name += '.txt'

    def __str__(self) -> str:
        return "{0} : {1}".format(self.name, self.path)
    


class error_log(log) :
    def __init__(self, name, path) -> None:
        super().__init__(name, path)
        self.errors = []
        self.failingfiles = []
        with open(self.path + self.name, 'w') as logfile :
            logfile.write('Error log : {0}\nCreated on : {1}\n\n'.format(self.name, self.creationdate))

    def write_error(self,filename: str, error: str) :
        with open(self.path + self.name, 'a') as logfile :
            logfile.write("{0} ERROR {1} : {2}\n".format(get_datetime(),filename, error))

    def add_error(self,filename: str, error: str, msg: str) :
        self.failingfiles += [filename]
        self.errors += [error]
        self.write_error(filename= filename, error= msg)

    def output_errors(self) :
        Errors = pd.DataFrame(columns= ['rootfilename'], data= self.errors)
        Errors.to_feather(self.path + 'Errors')


class parameter_log(log) :
    def __init__(self, name, path,) -> None:
        super().__init__(name,path)
        self.parameters = {}
        with open(self.path + self.name, 'w') as logfile :
            logfile.write('Parameter log : {0}\nCreated on : {1}\n\n'.format(self.name, self.creationdate))

    def add_parameters(self, *parameters) :
        for parameter in parameters :
            self.parameters[_get_varname(parameter)] = parameter

    def write(self) :
        with open(self.path + self.name, 'a') as logfile : 
            for parameter in self.parameters :
                logfile.write("{0} : {1}\n".format(parameter, self.parameters[parameter]))
        del self.parameters
        self.parameters = {}

class run_log(log) :
    def __init__(self, name, path) -> None:
        super().__init__(name, path)