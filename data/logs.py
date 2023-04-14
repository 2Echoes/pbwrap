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
        Errors = pd.DataFrame(columns= ['rootfilename', 'error'], data= [self.failingfiles, self.errors])
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
        with open(self.path + self.name, 'w') as logfile :
            logfile.write("Log : {0}\nCreated on : {1}\n\n".format(self.name, self.creationdate))

    def update(self, error_log: error_log, current_acquisition: int, filename: str, acquisition_start= 0) :
        """Updates the log file during analysis pipeline."""
        with open(self.path + self.name, 'w') as logfile :
            error_number = error_log.get_error_number()
            acquisition_number = current_acquisition - error_number - acquisition_start
            logfile.write("Log : {0}\nCreated on : {1}\n\n".format(self.name, self.creationdate))
            logfile.write("Updated on : {0}\n".format(get_datetime()))
            logfile.write("Number of acquisition that resulted in an error : {0}.\n".format(error_number))
            logfile.write("Number of acquisition processed successfully : {0}.\n".format(acquisition_number))
            logfile.write("Current acquistion : {0} - {1}.".format(acquisition_number, filename))

    def endrun(self, log_report: dict) :
        with open(self.path + self.name, 'w') as logfile :
            lines = ["Log : {0}\nCreated on : {1}\n\n".format(self.name, self.creationdate),
                     "Log finished on : {0} after a process time of {1}s.\n\n".format(get_datetime(), log_report['run time']),
                     "Total acquisition number : {0}\n".format(log_report['acquisition number'] + log_report['error number']),
                     "Success : {0}".format(log_report['acquisition number']),
                     "Error : {0}".format(log_report['error number']),
                     "Total cell detected : {0}".format(log_report['cell number']),
                     "\n### Integrity Checks ###\n",
                     "Acquisition DataFrame is empty : {0}.".format(log_report['Acquisition is_empty']),
                     "Acquisition DataFrame has valid id column : {0}.".format(log_report['AcquisitionId is_primarykey']),
                     "Cell DataFrame is empty : {0}.".format(log_report['Cell is_empty']),
                     "Cell DataFrame has valid id column : {0}.".format(log_report['CellId is_primarykey']),
                     "Cell defines (N,1) relation with Acquisition : {0}".format(log_report['Cell defines (N,1) relation with Acquisition'])
                     ]
            logfile.writelines(lines)