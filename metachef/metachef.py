#!/usr/bin/python

def msg():
### splash screen
    return """
    
                                                                
                     #####   ##    ##                           
                  ######  /#### #####                           
                 /#   /  /  ##### #####           #             
                /    /  /   # ##  # ##           ##             
                    /  /    #     #              ##             
                   ## ##    #     #      /##   ######## /###    
                   ## ##    #     #     / ### ######## / ###  / 
                   ## ##    #     #    /   ###   ##   /   ###/  
                   ## ##    #     #   ##    ###  ##  ##    ##   
                   ## ##    #     ##  ########   ##  ##    ##   
                   #  ##    #     ##  #######    ##  ##    ##   
                      /     #      ## ##         ##  ##    ##   
                  /##/      #      ## ####    /  ##  ##    /#   
                 /  #####           ## ######/   ##   ####/ ##  
                /     ##                #####     ##   ###   ## 
                #                                               
                 ##                                             
                                                                
                                                                
                                                              
                        # ###      /                    /##   
                      /  /###  / #/                   #/ ###  
                     /  /  ###/  ##                  ##   ### 
                    /  ##   ##   ##                  ##       
                   /  ###        ##                  ##       
                  ##   ##        ##  /##      /##    ######   
                  ##   ##        ## / ###    / ###   #####    
                  ##   ##        ##/   ###  /   ###  ##       
                  ##   ##        ##     ## ##    ### ##       
                  ##   ##        ##     ## ########  ##       
                   ##  ##        ##     ## #######   ##       
                    ## #      /  ##     ## ##        ##       
                     ###     /   ##     ## ####    / ##       
                      ######/    ##     ##  ######/  ##       
                        ###       ##    ##   #####    ##      
                                        /                     
                                       /                      
                                      /                       
                                     /                        

    MetaChef -- simulate metatranscriptomic datasets from real or simulated reads version 0.2

    """



### simple class for holding some values for each sample
class Samples:
    def __init__(self, metaline):
        self.species, self.abu, self.fw, self.rv, self.gff, self.fna = metaline
        self.abu = float(self.abu)
        self.depth = 0.0        
        self.nabu = 0.0
        self.ndepth = 0.0
        self.frac = 0.0

### input file parser
def parse_meta_input(meta, run):
    samples = []
    #"species","abundance","reads","gff"]
    with open(meta,"r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            pline = line.strip().split("\t")
            if run != "all":
                s = Samples(pline[:-1])
                s.depth = float(pline[-1])
            else:
                s = Samples(pline)
            samples.append(s)
    return samples

### computes total lenght for one fastq or fastq.gz file
def fastq_len(fastq):
    if fastq.endswith(".gz"):
        command = "zcat " + fastq + " | awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}'"
    else:
        command = "awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' " + fastq
    read_len = check_output(command, shell = True)
    return int(read_len)

### computes total lenght of selected features from gff file
def transcripts_len(gff):
    tot_len = 0
    if gff.endswith(".gz"):
        fun = gzip.open
    else:
        fun = open
    for line in fun(gff,"r"):
        if line.startswith("#"):
            continue
        l = line.strip().split("\t")
        if l[2] in ["CDS","ncRNA","rRNA","tRNA","tmRNA"]:
            tot_len += abs(int(l[3]) - int(l[4]))
    return tot_len 

### wraps fastq_len & transcripts_len to estimate coverage for a RNA-Seq experiment
def cov_est(reads, gff):
    fql = 0    
    for f in reads:
        fql += fastq_len(f)
    trl = transcripts_len(gff)
    cov = float(fql)/trl
    return cov

def main():
### splash screen just because
    print msg()
    # check mandatory argument
    if args.meta == None:
        
        print """    
    --meta argument missing

    run python metachef.py -h for help
              """  
        quit()
### setqk path
    seqtk = "/home/ldemarco/MetaTrans/tools/seqtk-master/seqtk"
### get input
### we didn't compute depth beforehand
    if args.run == "all":
        samples = parse_meta_input(args.meta, args.run)
        outfile = open(args.out + "metachef_metadepth.tsv", "w")
        ### for each sample
        for s in samples:
            ### get reads file paths
            reads = [s.fw, s.rv]
            ### get depth and update Samples class object 
            s.depth = cov_est(reads, s.gff)
            outfile.write("\t".join([s.species, str(s.abu), s.fw, s.rv, s.gff, str(s.depth)]) + "\n")      
    else:
### we already have computed depth
        samples = parse_meta_input(args.meta, args.run)
    # normalize depth and abundance into 0 - 1 range (sum --> 1.0)
    norm_a = sum([s.abu for s in samples])
    norm_d = sum([s.depth for s in samples])       
    for s in samples:
        s.nabu = s.abu/norm_a
        s.ndepth = s.depth/norm_d         
    # determine sample with highest ratio between normalized abundance and depth; get ratio between its depth and its normalized abundance --> factor     
    factor = max([(s.nabu/s.ndepth, s.depth/s.nabu) for s in samples])[1]
    # use factor to determine fraction of reads for downsampling; write down log file          
    with open(args.out + "metachef.log","w") as outfile:
        outfile.write("\t".join(["species", "abundance","normalized abundance","estimated depth", "normalized depth"]) + "\n")
        for s in samples:
            s.frac = (factor * s.nabu) / s.depth
            outfile.write("\t".join([str(x) for x in [s.species, s.abu, s.nabu, s.depth, s.frac * s.depth]]) + "\n")
    # compute normalization factor --> total coverage if its abundance was 100%
    for s in samples:
        reads = [s.fw, s.rv]
        for r in reads:
            n = reads.index(r) + 1
            o = args.out + "metachef_" + str(n) + ".fastq.gz"
            os.system("{} sample {} {} | gzip >> {}".format(seqtk,r,str(s.frac),o))    


        
        

        
        
    



if __name__ == "__main__":
    import gzip
    import argparse
    #import argcomplete
    import os
    from glob import glob
    from subprocess import check_output
    parser = argparse.ArgumentParser(description = "MetaChef -- simulate metatranscriptomic datasets from real or simulated reads", usage = msg())
    parser.add_argument("--meta", "-m", help = "metadata file -- MANDATORY  ")
    parser.add_argument("--out", "-o", default = os.getcwd() + "/" ,help = "output folder. default --> current directory")
    parser.add_argument("--run", "-r", default = "all", help = "all, depth. default --> all")
    #argcomplete.autocomplete(parser)
    args = parser.parse_args()
    main()    
