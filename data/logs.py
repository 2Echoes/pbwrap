from .getdata import get_datetime, _get_varname
from .output import dict_to_lines
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


    def get_error_number(self) :
        return len(self.errors)

    def write_error(self,filename: str, error: str) :
        with open(self.path + self.name, 'a') as logfile :
            logfile.write("{0} ERROR {1} : {2}\n".format(get_datetime(),filename, error))

    def add_error(self,filename: str, error: str, msg: str) :
        self.failingfiles += [filename]
        self.errors += [error]
        self.write_error(filename= filename, error= msg)

    def output_errors(self) :
        if len(self.errors) != 0 :
            Errors = pd.DataFrame(columns= ['rootfilename', 'error'], data= zip(self.failingfiles, self.errors))
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
        self.sucess_count = 0
        self.error_count = 0
        with open(self.path + self.name, 'w') as logfile :
            logfile.write("Log : {0}\nCreated on : {1}\n\n".format(self.name, self.creationdate))

    def sucess(self):
        self.sucess_count += 1
    
    def failure(self):
        self.error_count += 1

    def update(self, filename: str, rna_computed: list) :
        """Updates the log file during analysis pipeline."""
        with open(self.path + self.name, 'w') as logfile :
            logfile.write("Log : {0}\nCreated on : {1}\n\n".format(self.name, self.creationdate))
            logfile.write("Updated on : {0}\n".format(get_datetime()))
            logfile.write("Number of acquisition that resulted in an error : {0}.\n".format(self.error_count))
            logfile.write("Number of acquisition processed successfully : {0}.\n".format(self.sucess_count))
            logfile.write("Current acquistion : - {0}.\n".format(filename))
            logfile.write("So far {0} gene(s) have been analysed.\nList of analysed gens : {1}\n".format(len(rna_computed), rna_computed))

    def endrun(self, log_report: dict) :
        with open(self.path + self.name, 'w') as logfile :
            lines = ["Log : {0}\nCreated on : {1}\n\n".format(self.name, self.creationdate),
                     "Log finished on : {0} after a process time of {1}s.\n\n".format(get_datetime(), log_report['run time']),
                     "Total acquisition number : {0}\n".format(self.sucess_count + self.error_count),
                     "Success : {0}\n".format(self.sucess_count),
                     "Error : {0}\n".format(self.error_count),
                     "Total cell detected : {0}\n".format(log_report['cell number']),
                     "\n### Integrity Checks ###\n",]
            
            del log_report['run time'], log_report['cell number']
            lines += list(dict_to_lines(log_report))

            logfile.writelines(lines)