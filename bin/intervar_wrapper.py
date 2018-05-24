#!/Users/sjung/anaconda3/bin/python

import optparse, os, sys, vcf, tempfile, subprocess, shutil
from subprocess import *

CHUNK_SIZE = 2**20 #1mb

#def gnomad2vcf(gene):
#    return gene


def vcf2bed(inputF):
    vcf_to_bed=[]
    vcf_reader = vcf.Reader(open(inputF, 'r'))
    for record in vcf_reader:
        if len(record.ALT) == 1:
        #print(','.join(str(v) for v in record.ALT))
            #line = "%d\t%d\t%d\t%s\t%s" % (int(record.CHROM),int(record.POS),int(record.POS),record.REF, ','.join(str(v) for v in record.ALT))
            line = "%s\t%s\t%s\t%s\t%s" % (record.CHROM,record.POS,record.POS,record.REF, record.ALT[0])
            vcf_to_bed.append(line)
        else:
            # if multiple alternative alleles, print each in a separate line
            for i in range(0,len(record.ALT)):
                line = "%s\t%s\t%s\t%s\t%s" % (record.CHROM,record.POS,record.POS,record.REF, record.ALT[i])
                vcf_to_bed.append(line)
    return(vcf_to_bed)

def job_proc(cmd, tmp_dir, intervar_path):
    stdout = tempfile.NamedTemporaryFile( prefix="intervar-stdout-", dir=tmp_dir )

    stderr = tempfile.NamedTemporaryFile( prefix="intervar-stderr-", dir=tmp_dir )
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=intervar_path)
    return_code = proc.wait()

    if return_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
        stderr.flush()
        stderr.seek(0)
    #while True:
    #    chunk = stderr.read( CHUNK_SIZE )
    #    print(chunk)
    #    if chunk:
    #        stderr_target.write( chunk )
    #    else:
    #        break
    stderr.close()
    stdout.close()

def post_proc(intervarF, outname):
    #intervarF="/var/folders/gq/y3z1z7dd6vjdf2kpryyp6h_h0000gn/T/intervar-jg4_f9n7/cep290.hg19_multianno.txt.intervar"
    intervar_reader = open(intervarF, 'r')
    header = intervar_reader.readline()
    intervar_only=[]
    intervar_only.append(header.rstrip("\n"))
    intervar_clinvar=[]
    intervar_clinvar.append(header.rstrip("\n"))
    clinvar_only=[]
    clinvar_only.append(header.rstrip("\n"))
    for line in intervar_reader:
        values = line.rstrip("\n").split("\t")
        clinvar, intervar = values[12].lower(), values[13].lower()

        if "pathogenic" in intervar:
            intervar_only.append(line.rstrip("\n"))

        if "pathogenic" in clinvar and "pathogenic" in intervar:
            intervar_clinvar.append(line.rstrip("\n"))
        elif "pathogenic" in clinvar and "pathogenic" not in intervar:
            clinvar_only.append(line.rstrip("\n"))
        else:
            continue

    with open("%s.pathogenic.intervar.txt" % outname, "w") as f1, open("%s.pathogenic.clinvar_only.txt" % outname, "w") as f2, open("%s.pathogenic.both_in_intervar_clinvar.txt" % outname,"w") as f3:
        f1.writelines(["%s\n" % item  for item in intervar_only])
        f2.writelines(["%s\n" % item  for item in clinvar_only])
        f3.writelines(["%s\n" % item  for item in intervar_clinvar])
    shutil.copy(intervarF,"%s.intervar" % outname)

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('-i','--inputF', dest='inputF', action='store', type='string', help='Input vcf file')
    parser.add_option('-b','--bedF', dest='bedF', action='store', type='string', help='Output bed file')
    parser.add_option('-o','--outputF', dest='outputF', action='store', type='string', help='Output intervar file')
    (options, args) = parser.parse_args()

    intervar_bin="/Users/sjung/Documents/GitHub/pavar/InterVar/Intervar.py"
    intervar_path="/Users/sjung/Documents/GitHub/pavar/InterVar/"
    tmp_dir = tempfile.mkdtemp(prefix="intervar-")
    print(tmp_dir)
    if not os.path.exists(options.inputF):
        print("input file does not exist!")
        exit()
    else:
        inputF=options.inputF
    bedF = vcf2bed(inputF)

    bedOutput = "%s/%s" % (tmp_dir, options.bedF); print(bedOutput)
    intervarOutput = "%s/%s.tmp" % (tmp_dir, options.bedF); print(intervarOutput)

    with open(bedOutput,"w") as fout:
            fout.writelines(["%s\n" % item  for item in bedF])

    cmd = "%s -b hg19 -i %s -o %s" % (intervar_bin, bedOutput, intervarOutput)
    print(cmd)
    job_proc(cmd, tmp_dir,intervar_path)
    intervarF=intervarOutput+".hg19_multianno.txt.intervar"
    post_proc(intervarF,options.outputF)
    cwd = os.getcwd()
    #shutil.copy(intervarF,cwd)
    shutil.rmtree(tmp_dir)

if __name__=="__main__": __main__()
