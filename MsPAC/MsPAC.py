#!/bin/env python
import os
import sys
import ConfigParser
from string import Template
from lsf.lsf import Lsf

class Pipeline(object):
    def __init__(self,configfile,step):
        self.configfile = configfile
        self.step = step
        self.jobs = []
        
        #### All steps params
        ## User mandatory input
        self.directory = None
        self.package_bash_directory = None
        self.package_python_directory = None
        self.cluster = None
        self.job_threads = None
        self.job_walltime = None
        self.job_memory = None
        self.job_queue = None

        #### Multi steps params
        ## User optional input
        self.phased_bamfile = None
        self.phased_bedfile = None

        #### phase-bam only params
        ## User mandatory input
        self.vcffile = None
        self.bamfile = None
        self.vcf_sample_name = None
        ## User optional input
        self.generate_phased_bedfile = None
        self.max_phased_window = None
        self.padding_size = None

        #### prep reads only params
        self.raw_reads_in_bam_format = None
        self.raw_reads_directory = None

        #### assembly only params
        ## User optional input
        self.min_phased_block = None
        self.windows_to_assemble = None
        self.haps_to_assemble = None
        self.assembly_directory = None
        self.flanking_length = None
        self.max_block_length = None
        ## Non user options
        self.hap0_assembly_fa = None
        self.hap0_assembly_fq = None
        self.hap1_assembly_fa = None
        self.hap1_assembly_fq = None
        self.hap2_assembly_fa = None
        self.hap2_assembly_fq = None

        #### sv calling only params 
        self.sv_calling_directory = None
        self.reference = None
        ## Non user options
        self.hap1_assembly_split_fq = None
        self.hap2_assembly_split_fq = None

    def set_options(self,pipeline_options):
        for pipeline_option, option_input in pipeline_options:
            setattr(self,pipeline_option,option_input)        

    def get_bash_scripts_path(self):
        return "%s/MsPAC/bash" % "/".join(os.path.dirname(__file__).split("/")[:-1])

    def get_python_scripts_path(self):
        return "%s/MsPAC/python" % "/".join(os.path.dirname(__file__).split("/")[:-1])

    def set_optional_options(self,config):
        if config.has_option("Params","chromosome"):
            self.chromosome = config.get("Params","chromosome")

    def all_steps_configure(self,config):
        pipeline_options  = \
            [("directory",os.path.abspath(config.get('Input','directory'))),
             ("package_python_directory",self.get_python_scripts_path()),
             ("package_bash_directory",self.get_bash_scripts_path()),
             ("cluster",config.get("Other params",'cluster')),
             ("job_threads",{"high":config.get("HIGH INTENSITY JOB","threads"),
                             "low":config.get("LOW INTENSITY JOB","threads")}),
             ("job_memory",{"high":config.get("HIGH INTENSITY JOB","memory"),
                            "low":config.get("LOW INTENSITY JOB","memory")}),
             ("job_walltime",{"high":config.get("HIGH INTENSITY JOB","walltime"),
                              "low":config.get("LOW INTENSITY JOB","walltime")}),
             ("job_queue",{"high":config.get("HIGH INTENSITY JOB","queue"),
                           "low":config.get("LOW INTENSITY JOB","queue")})]
        self.set_options(pipeline_options)
        self.create_directory(self.directory)

    def phase_bam_dependent_options(self):
        if self.phased_bedfile == "None":
            self.phased_bedfile = "%s.hap_blocks.bed" %  self.phased_bamfile[:-4]

    def phase_bam_configure(self,config):
        pipeline_options  = \
            [("vcffile",config.get("Phase-bam input files",'phased vcf')),
             ("bamfile",config.get("Phase-bam input files",'reads aligned')),
             ("vcf_sample_name",config.get("Phase-bam params",'sample name in VCF')),
             ("phased_bamfile",config.get("Phase-bam params",'output phased bamfile'))]
        self.set_options(pipeline_options)
        #self.phase_bam_dependent_options()

    def assembly_dependent_options(self):
        self.hap0_assembly_fa = "%s/hap0_assembly.fasta" % self.assembly_directory
        self.hap0_assembly_fq = "%s/hap0_assembly.fastq" % self.assembly_directory
        self.hap1_assembly_fa = "%s/hap1_assembly.fasta" % self.assembly_directory
        self.hap1_assembly_fq = "%s/hap1_assembly.fastq" % self.assembly_directory
        self.hap2_assembly_fa = "%s/hap2_assembly.fasta" % self.assembly_directory
        self.hap2_assembly_fq = "%s/hap2_assembly.fastq" % self.assembly_directory

    def assembly_configure(self,config):
        pipeline_options  = \
            [("min_phased_block",config.get("Assembly params",'Minimum phased block length')),
             ("phased_bedfile",config.get("Assembly params",'Phased bedfile')),
             ("windows_to_assemble",[]),
             ("haps_to_assemble",config.get("Assembly params",'Comma-seperated list of haplotypes')),
             ("phased_bamfile",os.path.abspath(config.get("Phase-bam params",'output phased bamfile'))),
             ("flanking_length",int(config.get("Assembly params",'Flanking length'))),
             ("max_block_length",int(config.get("Assembly params",'Max block length'))),
             ("raw_reads_directory",os.path.abspath(config.get("Prep reads params",'Raw reads directory'))), #ugh
             ("assembly_directory",os.path.abspath(config.get("Assembly params",'Assembly directory')))]
        self.set_options(pipeline_options)
        self.phase_bam_dependent_options()
        self.assembly_dependent_options()

    def prep_reads_configure(self,config):
        pipeline_options  = \
            [("raw_reads_in_bam_format",os.path.abspath(config.get("Prep reads params",'BAM fofn'))),
             ("raw_reads_directory",os.path.abspath(config.get("Prep reads params",'Raw reads directory'))),
             ("phased_bamfile",config.get("Phase-bam params",'output phased bamfile'))]
        self.set_options(pipeline_options)

    def sv_calling_dependent_options(self):
        self.hap1_assembly_split_fq = "%s/hap1_assembly_split.fastq" % self.sv_calling_directory
        self.hap2_assembly_split_fq = "%s/hap2_assembly_split.fastq" % self.sv_calling_directory

    def sv_calling_configure(self,config):
        pipeline_options  = \
            [("sv_calling_directory",os.path.abspath(config.get("SV calling params",'SV calling directory'))),
             ("assembly_directory",os.path.abspath(config.get("Assembly params",'Assembly directory'))),
             ("reference",os.path.abspath(config.get("SV calling params",'reference')))]
        self.set_options(pipeline_options)
        self.sv_calling_dependent_options()
        self.assembly_dependent_options()
        self.create_directory(self.sv_calling_directory)

    def configure(self):
        config = ConfigParser.RawConfigParser()
        config.read(self.configfile)
        self.all_steps_configure(config)        
        if self.step == "phase-bam":
            self.phase_bam_configure(config)
        if self.step == "assembly":
            self.assembly_configure(config)
        if self.step == "prep-reads":
            self.prep_reads_configure(config)
        if self.step == "sv-calling":
            self.sv_calling_configure(config)

    def create_directory(self,directory):
        if not os.path.exists(directory):
            os.makedirs(directory)

    def write_to_bashfile(self,template_bash,bashfile,params):
        filein = open(template_bash)
        src = Template(filein.read())
        output_lines = src.safe_substitute(params)
        bashfh = open(bashfile,'w')
        bashfh.write(output_lines)
        filein.close()
        bashfh.close()

    def non_emptyfile(self,checkfile):
        return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0

    def run_locally(self):
        use_cluster = self.cluster
        self.cluster = "no"
        self.submitjobs()
        self.cluster = use_cluster

    def submitjobs(self,wait=True):
        if len(self.jobs) == 0:
            return
        if self.cluster != "Yes":
            for job,intensity in self.jobs:
                os.system("sh %s" % job)
            self.jobs = []
            return
        hpc = Lsf()
        for job,intensity in self.jobs:
            hpc.config(cpu=self.job_threads[intensity],
                       walltime=self.job_walltime[intensity],
                       memory=int(self.job_memory[intensity]) * int(self.job_threads[intensity]),
                       queue=self.job_queue[intensity])
            hpc.submit("%s" % job)
        if wait:
            hpc.wait()
        else:
            dummy=1
            # Bug gets overwritten if told not to wait twice                                         
            #job_id_log = "%s/job.ids" % self.log_directory
            #hpc.write_ids(job_id_log)
        self.jobs = []

    def map_reads(self,reads,reference,directory,name,intensity="low"):
        bashfile = "%s/%s.sh" % (directory,name)
        template_bash = "%s/map_reads.sh" % self.package_bash_directory
        params = {
            'output': "%s/%s" % (directory,name),
            'threads': self.job_threads[intensity],
            'reads': reads,
            'ref': reference,
            'BLASR': self.fast_blasr
            }
        self.write_to_bashfile(template_bash,bashfile,params)
        return (bashfile,intensity)

    def __call__(self):
        if self.step == "phase-bam":
            from haplotype_assignment import HaplotypeAssignment
            print "Assigning reads to haplotypes..."
            assign_reads_to_haplotype = HaplotypeAssignment(self.configfile)
            assign_reads_to_haplotype.run()
        elif self.step == "prep-reads":
            from prepping_reads import PrepReads
            print "Prepping reads..."
            prep_reads = PrepReads(self.configfile)
            prep_reads.run()
        elif self.step == "assembly":
            from assemble_haplotype import HaplotypeAssembly
            print "Assembling haplotypes..."
            assemble_haps = HaplotypeAssembly(self.configfile)
            assemble_haps.run()
        elif self.step == "sv-calling":
            from sv_calling import SVCaller
            print "Calling SVs..."
            call_svs = SVCaller(self.configfile)
            call_svs.run()
        else:
            sys.exit("Choose one of the following steps: phase-bam, prep-reads, assembly, sv-calling")
        #print "Detecting structural variants..."
        #from calling_structural_variants import StructuralVariationDetection
        #calling_structural_variants = StructuralVariationDetection(self.configfile)
        #calling_structural_variants.run()

def run_pipeline(step,configfile):
    pipepine = Pipeline(configfile,step)
    return pipepine()

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: MsPAC <step> <config file>")
    return run_pipeline(sys.argv[1],sys.argv[2])
