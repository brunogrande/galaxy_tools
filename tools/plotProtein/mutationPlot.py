import json
import time
from subprocess import Popen, PIPE
import httplib2, sys
http = httplib2.Http(".cache")
server = "http://rest.ensembl.org/"

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--ensembl_gene','-g', dest='ensembl_gene',help='supply an ensembl gene ID to plot')
parser.add_argument('--maf_file','-m', dest='maf_file',help='supply a maf file as input')
parser.add_argument('--tick_spacing','-t',dest='tick_spacing',help="supply a distance to space ticks in image",default=100)

parser.add_argument('--output_pdf','-o',dest='out_file',help="supply a name for the output file")

#ENSG00000166888 STAT6

args = parser.parse_args()

print args
try:
    ensembl_gene = args.ensembl_gene
except AttributeError:
    print "error, no gene specified"
    exit()
try:
    maf = args.maf_file
except AttributeError:
    pass
tick_spacing = args.tick_spacing
out_file = args.out_file

mut_file = "%s_mutations.txt" % ensembl_gene
splice_file = "%s_splice.txt" % ensembl_gene
out_mut = open(mut_file,"w")
features_file = "%s_features.txt" % ensembl_gene
out_features = open(features_file,"w")
exon_file = "%s_exons.txt" % ensembl_gene
exon_out = open(exon_file,"w")
out_splice = open(splice_file,"w")

def parse_maf(maf_file,ensg):
    mafh = open(maf_file,"r")
    maf_data = []
    index_cols = {}
    for line in mafh:
        maf_dict = {}
        if line.startswith("#"):
            continue
        elif line.startswith("Hugo_Symbol"):
            #this is a header row
            header_vals = line.rstrip("\n").split("\t")
            i = 0
            for val in header_vals:
                index_cols[i] = val
                i=i+1
        else:
            vals = line.rstrip("\n").split("\t")
            i = 0
            for val in vals:
                maf_dict[index_cols[i]]= val
                i = i+1
            if maf_dict["Gene"] == ensg:
                maf_data.append(maf_dict)
    return maf_data
def overlap(s1,e1,s2,e2):
    '''Return true if two regions overlap'''
    if s1 == s2 and e1 == e2:
        return 1
    elif (s1 > s2 and e1 < e2) or (s2 > s1 and e2 < e1):
        return 1
    elif (s1 <= s2 and e1 >= s2) or (s2 <= s1 and e2 >= e1):
        return 1
    elif (s1 >= s2 and e1 <= e2) or (s2 >= s1 and e2 <= e1):
        return 1
    else:
        return 0
def getTranscriptLength(enst):
    time.sleep(0.5)
    ext = "/sequence/id/" + enst + "?type=cds"
    print ext
    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    if not resp.status == 200:
        print "Invalid response: ", resp.status
        print ext
        return 0
    #print json.dumps()
    #exit()
    decoded = json.loads(content)
    sequence = decoded["seq"]
    seq_len = len(sequence)
    return(seq_len)
    
def getProteinFeatures(ensp):
    time.sleep(0.5)
    ext = "/overlap/translation/" + ensp 
    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    if not resp.status == 200:
        print "Invalid response: ", resp.status
        print ext
        return 0
    #print json.dumps()
    #exit()
    decoded = json.loads(content)
    #print json.dumps(decoded, indent=4, sort_keys=True)
    return decoded
   
def getGeneDetails(ensg):
    time.sleep(0.5)
    ext = "/lookup/id/" + ensg + "?expand=1"
    print ext
    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    if not resp.status == 200:
        print "Invalid response: ", resp.status
        print ext
        return 0
    #print json.dumps()
    decoded = json.loads(content)
    #print json.dumps(decoded, indent=4, sort_keys=True)
    
    transcripts = decoded["Transcript"]
    longest_trans = ""
    longest_ensp = ""
    best_len = 0
    strand = ""
    for transcript in transcripts:
        enst = transcript["id"]
        strand = transcript["strand"] #-1 for - strand
        if transcript.has_key("Translation"):
            #print transcript["Translation"]
            ensp = transcript["Translation"]["id"]
            #print "checking %s" % ensp
            length = getTranscriptLength(enst)
            if length > best_len:
                longest_trans = enst
                longest_ensp = ensp
                best_len = length 
        else:
            continue
    #get exon positions in this transcript and cds start/end (relative to genome) from Translation
    exon_start_ends = {} #trimmed down to CDS, leaving out any purely UTR exons. Keyed on start position to allow easy sorting
    for transcript in transcripts:
        enst = transcript["id"]
        if enst == longest_trans:
            print "THIS IS THE LONGEST TRANSCRIPT: %s" % enst
            #exit()
            genomic_cds_start = transcript["Translation"]["start"]
            genomic_cds_end = transcript["Translation"]["end"]
            
            for exon in transcript["Exon"]:
                genomic_start = exon["start"]
                genomic_end = exon["end"] #note, even if on - strand, start is less than end
                if genomic_start > genomic_cds_end or genomic_end < genomic_cds_start:
                    continue
                elif genomic_start < genomic_cds_start:
                    genomic_start = genomic_cds_start
                elif genomic_end > genomic_cds_end:
                    genomic_end = genomic_cds_end
                #trimmed to CDS portion of exon
                exon_start_ends[genomic_start] = genomic_end
    exon_starts = exon_start_ends.keys()
    if strand >1:
        exon_starts.sort()
    else:
        exon_starts.sort(reverse=True) #puts them in the order of the protein
    exon_amino_acid_lengths = []
    tot =0
    prot_position = 0
    protein_len = float(best_len)/3
    exons_trimmed = []
    print "0\t%s\tprotein\tprotein" % protein_len
    for start in exon_starts:
        end = exon_start_ends[start]
        cds_len = end - start + 1
        prot_len = float(cds_len) / 3
        prot_start = prot_position
        prot_end = prot_len + prot_position
        #exons_trimmed.append(prot_start,prot_end)  #NEED TO TRACK WHERE IN THE GENOME EACH EXON IS TO ALLOW SPLICE SITE POSITION MAPPING>>>
        print "%s\t%s\texon\texon" % (prot_start,prot_end)
        line = "%s\t%s\n" % (prot_start,prot_end)
        exon_out.write(line)
        prot_position = prot_end 
        tot = tot+prot_len
    #note: although I have calculated the positions of exons, the current visualization script cannot indicate their positions. Needs to be implemented later
    features = getProteinFeatures(longest_ensp)
    merged_features = {}
    skip = {}
    for feature in features:
        if feature["description"] == "":
            continue
        if merged_features.has_key(feature["description"]):
            regions = merged_features[feature["description"]]
            i =  0
            no_overlap = 1
            for region in regions:
                start,end = region
                key="%s-%s" % (start,end)
                if overlap(start,end,feature["start"],feature["end"]):
                    all_four = [start,end,feature["start"],feature["end"]]
                    all_four.sort()
                    regions[i] = (all_four[0],all_four[3]) #min and max for region
                    no_overlap = 0
                i=i+1
            if no_overlap:
                regions.append((feature["start"],feature["end"]))
            
        else:
            merged_features[feature["description"]] = [(feature["start"],feature["end"])]
    colours = ('firebrick','gold','darkslategray3','darkgreen','deepskyblue3','maroon3','lightsalmon2','lightslategray','thistle4','violetred4')
    i = 0
    out_features.write("architecture_name\tstart_site\tend_site\tcolour\n")
    for feature_type in merged_features:
        colour = colours[i]
        i=i+1
        regions = merged_features[feature_type]
        for region in regions:
            print "%s\t%s\t%s\t%s" % (feature_type,region[0],region[1],colour)
            out_line = "%s\t%s\t%s\t%s\n" % (feature_type,region[0],region[1],colour)
            out_features.write(out_line)
    return (protein_len,exons_trimmed)
def getGeneMutations(ensg,exon_positions,maf=None):
    if maf:
        #parse mutations from MAF file instead of loading from db
        gene_symbol = ""
        maf_data = parse_maf(maf,ensg)
        print maf_data
        out_splice.write("ensg\tgene_symbol\tposition\tcolour\tmutation_type\n")
        for maf_dict in maf_data:
            mut_class = maf_dict['Variant_Classification']
            if mut_class == "Missense_Mutation":
                mut_class = "non-synonymous"
            elif mut_class == "Silent":
                mut_class = "synonymous"
            elif mut_class == "Nonsense_Mutation":
                mut_class = "nonsense"
            elif mut_class == "Splice_Site":
                mut_class = "splice"
            elif mut_class == "Frame_Shift_Del":
                mut_class = "indel"
            elif mut_class == "Frame_Shift_Ins":
                mut_class = "indel"
            else:
                continue #not handled yet
            #vcf2maf format
            if mut_class == "splice":
                positions = maf_dict["Protein_position"]
                (pos,length) = positions.split("/")
                try:
                    (pos,other) = pos.split("-")
                except ValueError:
                    pass
                print "%s\t%s\t%s\t%s\t%s" % (ensg,maf_dict["Hugo_Symbol"],pos,'blue',mut_class)
                line = "%s\t%s\t%s\t%s\t%s\n" % (ensg,maf_dict["Hugo_Symbol"],pos,'blue',mut_class)
                out_splice.write(line)
            elif mut_class == "indel":
                #print "THIS IS AN INDEL<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                positions = maf_dict["Protein_position"]
                (pos,length) = positions.split("/")
                ref = "-"
                nref = "fs"
                try:
                    (pos,other) = pos.split("-")
                except ValueError:
                    pass
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (ensg,maf_dict["Hugo_Symbol"],pos,ref,nref,'orange',mut_class)
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ensg,maf_dict["Hugo_Symbol"],pos,ref,nref,'orange',mut_class)
                gene_symbol = maf_dict["Hugo_Symbol"]
                out_mut.write(line)
            else:
                
                aa = maf_dict["Amino_acids"]
                if mut_class == 'synonymous':
                    ref = aa
                    nref = aa
                else:
                    (ref,nref) = aa.split("/")
                positions = maf_dict["Protein_position"]
                (pos,length) = positions.split("/")
                try:
                    (pos,other) = pos.split("-")
                except ValueError:
                    pass
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (ensg,maf_dict["Hugo_Symbol"],pos,ref,nref,'red',mut_class)
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ensg,maf_dict["Hugo_Symbol"],pos,ref,nref,'red',mut_class)
                gene_symbol = maf_dict["Hugo_Symbol"]
                out_mut.write(line)
        return gene_symbol
    else:
        import cancerGenome
        db = cancerGenome.cancerGenomeDB(database_name='lymphoma_full',database_host='jango.bcgsc.ca',database_user='rmorin',database_password='rmorin')
        db_relapse = cancerGenome.cancerGenomeDB(database_name='DLBCL_relapse',database_host='jango.bcgsc.ca',database_user='rmorin',database_password='rmorin')
        gene_ob = cancerGenome.Gene(db_relapse.db,ensembl_id=ensembl_gene)
        print "gene is %s,%s" % (gene_ob.gene_symbol,gene_ob.id)
        limits = {"ensembl_id":"= '%s'" % ensembl_gene}
        limits["disease"] = "='DLBCL'"
        muts = db.getMutations(limits=limits)
        for mut in muts:
            anno = mut.annotation
            annos = anno.split(";")
            anno = annos[-1]
            ref = anno[0:1]
            nref = anno[-1:]
            pos = anno[1:-1]
            mut_class = mut.type
            if nref == "*":
                mut_class = "nonsense"
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (ensg,gene_ob.gene_symbol,pos,ref,nref,'red',mut_class)
            line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ensg,gene_ob.gene_symbol,pos,ref,nref,'red',mut_class)
            out_mut.write(line)
        muts = db_relapse.getMutations(limits=limits)
        for mut in muts:
            anno = mut.annotation
            ref = anno[0:1]
            nref = anno[-1:]
            pos = anno[1:-1]
            mut_class = mut.type
            
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (ensg,gene_ob.gene_symbol,pos,ref,nref,'orange',mut_class)
            line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ensg,gene_ob.gene_symbol,pos,ref,nref,'orange',mut_class)
            out_mut.write(line)
        #splice = db_relapse.getSpliceSiteSNVs(gene=gene_ob)
        #nsplice = len(splice)
        
        #for snv in splice:
        #    positions = maf_dict["Protein_position"]
        #    (pos,length) = positions.split("/")
        #    try:
        #        (pos,other) = pos.split("-")
        #    except ValueError:
        #        pass
        #    print "%s\t%s\t%s\t%s\t%s" % (ensg,maf_dict["Hugo_Symbol"],pos,'blue',mut_class)
        #    line = "%s\t%s\t%s\t%s\t%s\n" % (ensg,maf_dict["Hugo_Symbol"],pos,'blue',mut_class)
        #    out_splice.write(line)
        #indels = db_relapse.getIndels(limits=limits)
        #for indel in indels:
        #    print "%s\t%s" % (indel,indel.annotation)
        #exit()
        
        
        return gene_ob.gene_symbol

(protein_length,exon_positions) = getGeneDetails(ensembl_gene)
gene_symbol = getGeneMutations(ensembl_gene,exon_positions,maf)
out_features.close()
out_mut.close()
exon_out.close()
out_splice.close()
#try actually running the script
command = "/home/bgrande/galaxy-dist/tools/gtools/plotProtein/plotProtein.R %s %s %s %s %s %s %s yes %s no" % (mut_file,features_file,exon_file,splice_file,protein_length,gene_symbol,tick_spacing,out_file)
args = command.split(" ")
print command
#run the command
process = Popen(args,stdout=PIPE, stderr=PIPE)
process.wait()
stdout, stderr = process.communicate()
print stdout
print stderr



