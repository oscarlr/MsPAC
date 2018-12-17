#!/bin/env python
from MsPAC import Pipeline

class CleanAssembly(Pipeline):
    def __init__(self,configfile):
        Pipeline.__init__(self,configfile,"clean-assembly")
    
    def map_raw_reads(self):
        template_bash = "%s/map_reads_for_quiver.sh" % self.package_bash_directory
        bashfile = "%s/map_reads_for_quiver.sh" % self.assembly_directory
        params = {
            'raw_reads_in_bam': self.raw_reads_in_bam_format,
            'contigs': self.assembly,
            'output': self.assembly_directory,
            'threads': self.job_threads["high"]
            }
        self.write_to_bashfile(template_bash,bashfile,params)
        self.jobs.append((bashfile,self.job_threads["high"]))
        self.submitjobs()

    def run(self):
        self.map_raw_reads()
        self.split_mappings()
        self.run_quiver()
