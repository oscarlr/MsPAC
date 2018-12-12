#!/bin/env python
from MsPAC import Pipeline

class HaplotypeAssignment(Pipeline):
    def __init__(self,configfile):
        Pipeline.__init__(self,configfile,"phase-bam")
        self.configure()
        
    def assign_reads_to_haplotype(self):
        bashfile = "%s_assign_reads_to_haplotype.sh" % self.bamfile[0:-4]
        template_bash = "%s/assign_reads_to_haplotype.sh" % self.package_bash_directory
        params = {
            'python_scripts': self.package_python_directory,
            'vcffile': self.vcffile,
            'bamfile': self.bamfile,
            'phased_bamfile': self.phased_bamfile,
            'vcf_sample_name': self.vcf_sample_name
            }
        self.write_to_bashfile(template_bash,bashfile,params)
        self.jobs.append((bashfile,"low"))
        self.run_locally()
        
    def run(self):
        self.assign_reads_to_haplotype()

