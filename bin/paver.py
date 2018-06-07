#!/Users/sjung/anaconda3/bin/python

"""
This program executes a program called InterVar, implementing the ACMG/AMP 2015 guideline,
to process a list of genes to identify (likely) pathogenic rare variants

command: ./bin/paver.py --inputF input/genelist.txt --ref gnomad/gnomad.genomes.r2.0.2.sites.chr1.vcf.gz --out_dir output_gene
         --inputF : input file that lists gene(s) symbols. Note that one gene per line is required.
         --ref    : vcf file that includes the variants call. Requird filed is CHROM\tPOS\tID\tREF\tALT
                    your gene symbol must present in the vcf file.
         --out_dir: this folder will be created and the output files will be generated.
"""

import optparse, os, sys, vcf, tempfile, subprocess, shutil, glob, gzip
from subprocess import *

CHUNK_SIZE = 2**20 #1mb

#def gnomad2vcf(gene):
#    return gene


def extract_gene(inputF, input_vcf, ref_file):
    try:
        with open(inputF, 'r') as geneList:
            for g in geneList:
                gene =g.rstrip("\n").upper()
                #print(gene)
                with open("%s/%s.vcf" % (input_vcf, gene), 'w') as gout, gzip.open(ref_file, 'rt') as ref:
                    gout.writelines(["%s" % line for line in ref if gene in line])
        ref.close()
    except:
        print("vcf generation is failed")
        raise


def vcf2bed(inputF, input_bed):
    for i in inputF:
        with open("%s/%s.bed" % (input_bed, os.path.splitext(os.path.basename(i))[0]),'w') as bedOut:
            vcf_reader = vcf.Reader(open(i, 'r'))
            try:
                if os.stat(i).st_size != 0:
                    for record in vcf_reader:
                        if len(record.ALT) == 1:
                            line = "%s\t%s\t%s\t%s\t%s\n" % (record.CHROM,record.POS,record.POS,record.REF, record.ALT[0])
                            bedOut.write(line)
                        else:
                            # if multiple alternative alleles, print each in a separate line
                            for i in range(0,len(record.ALT)):
                                line = "%s\t%s\t%s\t%s\t%s\n" % (record.CHROM,record.POS,record.POS,record.REF, record.ALT[i])
                                bedOut.write(line)
            except:
                print("vcf file is empty. Please check your genes and the reference vcf file")
                raise

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
    stderr.close()
    stdout.close()

def hwe_proc(cmd, tmp_dir):
    stdout = tempfile.NamedTemporaryFile( prefix="hwe-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="hwe-stderr-", dir=tmp_dir )
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir)
    return_code = proc.wait()

    if return_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
        stderr.flush()
        stderr.seek(0)
    stderr.close()
    stdout.close()

def post_proc(intervarF, outname, out_dir):
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

    with open("%s/%s.pathogenic.intervar.txt" % (out_dir, outname), "w") as f1, open("%s/%s.pathogenic.clinvar_only.txt" % (out_dir, outname), "w") as f2, open("%s/%s.pathogenic.both_in_intervar_clinvar.txt" % (out_dir, outname),"w") as f3:
        f1.writelines(["%s\n" % item  for item in intervar_only])
        f2.writelines(["%s\n" % item  for item in clinvar_only])
        f3.writelines(["%s\n" % item  for item in intervar_clinvar])
    shutil.copy(intervarF,"%s/%s.intervar" % (out_dir, outname))

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('','--inputF', dest='inputF', action='store', type='string', help='input file with a list of gene symbol')
    parser.add_option('','--ref', dest='ref', action='store', type='string', help='reference file')
    parser.add_option('','--out_dir', dest='out_dir', action='store', type='string', help='Output directory')

    (options, args) = parser.parse_args()

    intervar_bin="/Users/sjung/Documents/GitHub/pavar/InterVar/Intervar.py"
    hwe_bin="/Users/sjung/Documents/GitHub/pavar/bin/hwe.R"
    intervar_path="/Users/sjung/Documents/GitHub/pavar/InterVar/"
    tmp_dir = tempfile.mkdtemp(prefix="intervar-")
    #print(tmp_dir)
    if not os.path.exists(options.inputF):
        print("input file does not exist!")
        exit()
    else:
        inputF=options.inputF

    if os.path.exists(tmp_dir):
        input_vcf = "%s/input_vcf" % tmp_dir
        os.mkdir(input_vcf)
        input_bed = "%s/input_bed" % tmp_dir
        os.mkdir(input_bed)

    out_dir=options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    extract_gene(inputF, input_vcf, options.ref)
    VCFs = sorted(glob.glob("%s/*.vcf" % input_vcf))
    vcf2bed(VCFs, input_bed)

    #bedOutput = "%s/%s" % (tmp_dir, options.bedF); print(bedOutput)
    BEDs = sorted(glob.glob("%s/*.bed" % input_bed))

    for i in BEDs:
        intervarOutput = "%s.tmp" % (i); #print(intervarOutput)
        genename=os.path.splitext(os.path.basename(i))[0]

        cmd = "%s -b hg19 -i %s -o %s" % (intervar_bin, i, intervarOutput)
        #print(cmd)
        job_proc(cmd, tmp_dir,intervar_path)
        intervarF=intervarOutput+".hg19_multianno.txt.intervar"
        post_proc(intervarF,genename, out_dir)
    try:
        shutil.move("%s/input_vcf" %tmp_dir,out_dir)
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

    for j in VCFs:
        genesymbol=os.path.splitext(os.path.basename(j))[0]
        paverfile= "%s/%s.pathogenic.intervar.txt" % (out_dir,genesymbol)
        vcffile= "%s/input_vcf/%s.vcf" % (out_dir,genesymbol)
        hwe= "%s/%s_hwe.txt" % (out_dir,genesymbol)
        cmd = "Rscript %s -v %s -i %s -o %s" % (hwe_bin, vcffile, paverfile, hwe)
        #print(cmd); print(out_dir)
        hwe_proc(cmd, out_dir)

    #shutil.rmtree(tmp_dir)






if __name__=="__main__": __main__()
