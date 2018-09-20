#!/bin/sh
''''exec python -u -- "$0" ${1+"$@"} # '''
# vi: syntax=python

from __future__ import division
"""This script takes candidate polymorphisms from LAST alignments.
From the raw counts of nucleotides, the script validates SNPs by relative frequencies. These are then weighted by the probabilities
for base transitions (transitions vs. transversions). For human data, SNP positions are compared to dbSNP, extracting also reported
alleles and rs-IDs."""


import sys
import os
import re
import json
from subprocess import Popen, PIPE, STDOUT
import time

"""Functions"""

def calculateBchange(nuc, target, tt_ratio):
    # At equal mut freq: transv more likely than transitions!!!  
    """Returns score for mutation of target to query nucleotide.
    Based on p(transition) = tt_ratio * p(transversion)."""
    score = 1
    transition_list = [("a", "g"), ("g", "a"), ("c", "t"), ("t", "c")]
    transversion_list = [("a", "c"), ("c", "a"), ("a", "t"), ("t", "a"), ("c", "g"), ("g", "c"), ("g", "t"), ("t", "g")]
    if (nuc, target) in transition_list:
        score = score * tt_ratio
        return score
    elif (nuc, target) in transversion_list:
        return score
    # remove it from polymorphism table
    elif nuc == target:  
        score = 0
        return score
    else:
        return "Error: Non-DNA character encountered"


def getSNPwww(queries):

    """A function to look up base multiple positions in the new API of dbSNP."""

    import urllib2
    import urllib
    
    out_dict = {}
    isDBerror = False
    
    try:
        
        query_str =""
        assembly = "GCF_000001405.38"
       
        # Create post request string for API. Example: query_str = "10 557773 iD G N\n 10 557775 iD C C"
        # 50,000 positions can be queried at once. I use 40,000 for security.

        for i in range (0, len(queries), 40000):
            l=queries[i:i+40000]
            query_str= "".join(l)

            api_link = "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/file/set_rsids?assembly=" + assembly
            
            # Try if server responds
            try:
                site = urllib2.Request(api_link, query_str) 
                resp = urllib2.urlopen(site).read().strip("\n")
                
                # Each list entry is one position
                out_list = resp.split("\n") 
            except urllib2.HTTPError:
                isDBerror = True
                print "HTTPError on website  %s." % api_link
                break
    
            except urllib2.URLError:
                isDBerror = True
                print "URLError on website %s." % api_link
                break
            
            else:
                # Retrieve and process results
                for pos in out_list:    
                    pos_data = pos.split("\t")
                    
                    if len(pos_data) >= 5:
                        chr = pos_data[0]
                        base_pos = pos_data[1]
                        rs_id = pos_data[2]
                        allele = pos_data[3]+"/"+pos_data[4]
                        
                        # Get rid of NORSID and second part of error message
                        if re.match(r"rs", rs_id):
                            isID = False
                            
                            if chr in out_dict:
                                
                                if base_pos in out_dict[chr]:
                                    
                                    #Get to the tupels
                                    for list_element in out_dict[chr][base_pos]:# outdict[base position]=[(rs_id, [allels])]
                                        if rs_id in list_element[0]:
                                            isID = True
                                            
                                            # Avoid multiple entries of the same alleles for the same rs_ID
                                            if allele in list_element[1]:  
                                                break
                                            else:
                                                list_element[1].append(allele)
                                                
                                    # The rs_ID was not at all in the list of values
                                    if isID == False:
                                        out_dict[chr][base_pos].append((rs_id, [allele]))
                                        
                                else:
                                    out_dict[chr][base_pos] = [(rs_id, [allele])]
                                    
                            else:
                                out_dict[chr] ={base_pos: [(rs_id, [allele])]}
                        
                    # Write error message from API
                    else:
                        print "".join(pos)
                        
            # At max two queries per second
            time.sleep(0.5)
        
    finally:
        return out_dict, isDBerror


def getSNPplas(chr_numb, print_dict):
    
    """Check candidate SNP positions in print_dict for reports in the local PlasmoDB database. The database is available for Plasmodium falciparum
    and the version 3 (_v3). All IDs in a database file have to be sorted according to the position!
    The function will return a dictionary: outdict[snp]=[(dbID, [dbMaj: dbMajF, dbMin: dbMinF])]"""
    
    outdict = {}
    snp_list = print_dict.keys()
    dbPath = "/bioinf/projects/SNPdbPlasf/"+chr_numb+".txt"
    with open(dbPath, "r") as plasmDB:
        for line in plasmDB:
            
            # Skip header
            if not "[" in line[0]:
                
                # Get ID
                dbID = line.split("\t")[0]
                dbPos = dbID.split(".")[-1]
                
                # Check if SNPs match
                for snp in snp_list:
                    if snp == dbPos:
                        
                        # Major allele and its frequency
                        dbMaj = line.split("\t")[1]
                        dbMajF = line.split("\t")[2]
                        
                        # Minor allele and its frequency
                        dbMin = line.split("\t")[3]
                        dbMinF = line.split("\t")[4]
                        
                        # Write to outdict
                        if snp in outdict:
                            outdict[snp].append((dbID, [dbMaj + ":" + dbMajF, dbMin + ":" + dbMinF]))
                        else:
                            outdict[snp]=[(dbID, [dbMaj + ":" + dbMajF, dbMin + ":" + dbMinF])]
                            
                    # Remove SNPs from the list that can't be found (anymore)      
                    elif int(snp) < int(dbPos):
                        snp_list.remove(snp)       
    return outdict



res_dict={}


# Variables
isQualityAnaly = False
isRep = False
isDBerror = False
target_threshold = 0.8
poly_threshold = 0.2
# Transitions are twice as likely as transversions
tt_ratio = 2  
cover_threshold = 0.3
# Threshold for printing raw coverage
rawCovThresh = 1 

#Get species argument from command line for dbSNP check
## script.py -s human ...

# take arguments, remove script name
arg_list = sys.argv[1:] 
if "-s" in arg_list:
    arg_index = arg_list.index("-s")
    org_index = arg_index + 1
    try:
        if arg_list[org_index] == "human":
            organism = "Homo sapiens"
        elif arg_list[org_index] == "plasf":
            organism = "Plasmodium falciparum"
        else:
            # Might be useful for other databases. Should then be more precise.
            if arg_list[org_index] != "-q":
                organism = "Non-human organism (" + arg_list[org_index] + ")" 
            else:
                organism = "Undetected organism"
                
    except:
        organism = "Undetected organism"
else:
    organism = "Undetected organism"

if "-q" in arg_list:
    isQualityAnaly = True

    
# Get encodings for chromosomes and scaffolds
chr_enc_dict = {}
chr_to_enc_dict = {}
with open("calc.tidmap", "r") as chr_file:
    for line in chr_file:
        line = line.strip("\n")
        line_list = line.split("\t")
        chr_ = line_list[0]
        chr_ = chr_.replace("chr", "")
        chr_enc = line_list[1]
        chr_enc_dict[chr_enc] = chr_
        chr_to_enc_dict[chr_] = chr_enc

# The tidmap file was empty, because the data was insufficient. Exit without raising error.
if not chr_enc_dict:
    sys.exit(0)        

# Get nuccount files from current directory, ignore ".help" files
curr_dir = os.getcwd()
for root, dirs, files in os.walk(curr_dir, topdown = True):
    file_list = files
    break

last_file = file_list[-1]

for file_name in file_list:
    
    # Repetitions of erroneous queries after last element from first file_list was processed
    if file_name == last_file:
        isRep = True
    match = re.search(r"^calc.nuccounts.\d+$", file_name)
    if match:
        print_dict = {}
        cover_dict = {} # for discarding low coverage data
        raw_cov_dict = {} #for printing raw coverage
        
        #sys.stdout.write("Processing: %s: " % file_name)
        
        # Read file....
        with open(file_name, "r") as alignment_file:
            for line in alignment_file.readlines():
                print_list = []
                nuc_dict = {}
                isPoly = False
                isSPECIALchar = False
                poly_dict = {}
                weight_dict = {}
                result_list = []
                line = line.replace("\n", "")
                line_data = line.split("\t")
                
                # .... and extract data
                if ">" in line_data[0]:
                    chr_numb = line_data[0][-1]
                    continue
                pos_start = line_data[0]
                nuc_dict["a"] = int(line_data[1])
                nuc_dict["c"] = int(line_data[2])
                nuc_dict["g"] = int(line_data[3])
                nuc_dict["t"] = int(line_data[4])
                consensus = line_data[5].lower()
                target = line_data[6].lower()
                
                # Total amount of nucleotides per position
                total_nuc = 0
                for nuc_count in nuc_dict.values():
                    total_nuc = total_nuc + nuc_count
                    
                # Save coverage for all positions
                if total_nuc in cover_dict:
                    cover_dict[total_nuc].append(pos_start)
                else:
                    cover_dict[total_nuc] = [pos_start]
                
                if consensus not in ["-", "n"]:
        
                    # SNP is accepted if: (1) target nuc < 0.8 [isPoly = True]; (2) one other nuc > 0.2 [poly_threshold]
                    ## Special consensus symbols according to IUPAC and "X" pass (1) automatically, (2) also fulfilled
                    if consensus in ["m", "r", "w", "s", "y", "k", "v", "h", "d", "b", "x"]:
                        isPoly = True
                    else:
                        if target in ["a", "c", "g", "t"]:
                            if nuc_dict[target] / total_nuc <= target_threshold:
                                isPoly = True
                        else:
                            isPoly = False
        
                    ## Check for condition (2)       
                    if isPoly == True:
                        for nuc in nuc_dict:   
                            rel_nuc = nuc_dict[nuc] / total_nuc
                            rel_nuc = round(rel_nuc, 3)
                            if rel_nuc >= poly_threshold:
                                poly_dict[nuc] = [rel_nuc, calculateBchange(nuc, target, tt_ratio)]  # obtain score for transitions vs transversions
        
                    # Give weight to the relative occurence of a nuc
                    if poly_dict:
                        
                        for nuc in poly_dict:
                            poly_dict[nuc] = poly_dict[nuc][0] * poly_dict[nuc][1]
        
                        # Rescale weighted rel. occurences to 1
                        weight_factor = sum(poly_dict.values()) / 1
                        if weight_factor == 0:
                            continue
                        for nuc in ["a", "c", "g", "t"]:
                            if nuc in poly_dict:
                                resc_prob = round(poly_dict[nuc] / weight_factor, 3)
                                if resc_prob > 0:
                                    print_list.append(str(resc_prob)+"\t")
                                else:
                                    print_list.append("-\t") 
                            else:
                                print_list.append("-\t")
                        print_list.append(target+"\t")
                        print_list.append("\n")
                        print_dict[pos_start] = print_list
                        raw_cov_dict[pos_start] = nuc_dict
                
                else:
                    continue
               
        # Following calculations only for polymorphic nuccount files
        if cover_dict and print_dict: 
            
            # Get the highest coverage of an SNP within a file, calculate a minimum coverage for every SNP
            cover_list = cover_dict.keys()
            cover_list.sort()
            highest_cover = cover_list[-1]
            min_cover = highest_cover * cover_threshold
            
            for coverage in cover_dict:
                if coverage < min_cover:
                    for position in cover_dict[coverage]:
                        if position in print_dict:
                            del print_dict[position]
                            
            if print_dict:  
                     
                # Get reference data for suspect SNPs from dbSNP
                chr_enc = str(file_name.split(".")[-1])  # file name: calc.nuccounts.chr_enc
                chr_numb = chr_enc_dict[chr_enc] #1
                
                
                
                if organism != "Plasmodium falciparum" and organism != "Homo sapiens":
                        print file_name + ": %s, %s not in SNP databases." % (organism, chr_numb)
                        res_dict["**header**"] = "Position\tA\tC\tG\tT\tTarget\traw A\traw C\traw G\traw T\n"
               
                else:
                     for pos in print_dict:
                            print_dict[pos].pop(-1)
                            print_dict[pos].append("\t")
                            
                #Plasmodium falciparum            
                if organism == "Plasmodium falciparum":
                    res_dict["**header**"] = "Position\tA\tC\tG\tT\tTarget\tMatches in PlasmoDB\traw A\traw C\traw G\traw T\n"
                    try:
                        chr_numb = chr_numb.split(":")[0]
                    except:
                        pass
                    
                    # Local database features Plasmodium falciparum v3
                    if re.search(r"^Pf3D7_\d\d_v3", chr_numb):
                        plas_dict = getSNPplas(chr_numb, print_dict)
                    
                        if plas_dict:
                            suc = 0
                            fail = 0
                            for plas_pos in plas_dict:
                                plas_pos = str(plas_pos)
                                if plas_pos in print_dict:
                                    print_dict[plas_pos].pop(-1)
                                    rs_counter = 0
                                    
                                    # List_element is a tupel: (plasID, [dbMaj: dbMajF, dbMin: dbMinF])
                                    plas_str = ""
                                    for list_element in plas_dict[plas_pos]:
                                        rs_counter += 1
                                        plas_ID = list_element[0]
                                        plas_str = plas_str + plas_ID + ": "
                                        allele_counter = 0
                                        for plas_allele in list_element[1]:
                                            allele_counter += 1
                                            plas_str = plas_str + plas_allele                                            
                                            # The last allele should not be followed by a separator
                                            if allele_counter < len(list_element[1]):
                                                plas_str = plas_str + " + "
                                        if rs_counter < len(plas_dict[plas_pos]):
                                            plas_str = plas_str + "; "
                                            
                                    print_dict[plas_pos].append(plas_str+"\t")
                                    suc += 1
                                else:  # Translation of dbSNP rs-ID to base position not successful = request artifact
                                    fail += 1
                            print file_name + ": PlasmoDB local: %s match(es); %s artifact(s)." % (str(suc), str(fail))
                        else:
                           print file_name + "PlasmoDB local: 0 matches."
                           
                    else:
                        print file_name + ": %s, %s not in local PlasmoDB." % (organism, chr_numb)
                        for pos in print_dict:
                            print_dict[pos].pop(-1)
                            print_dict[pos].append("N/A\t")
                        
                
                    
                        
                # Write raw coverage to output
                for pos in print_dict:
                    if print_dict[pos][-1] == "\n":
                            print_dict[pos].pop(-1)

                    for nuc in ["a", "c", "g", "t"]:
                        if raw_cov_dict[pos][nuc] >= rawCovThresh:
                            if nuc == "t":
                                print_dict[pos].append(str(raw_cov_dict[pos][nuc])+"\n")
                            else:
                                print_dict[pos].append(str(raw_cov_dict[pos][nuc])+"\t")
                        else:
                            if nuc == "t":
                                print_dict[pos].append("-\n")
                            else:
                                print_dict[pos].append("-\t")

                
                # Save all files in dictionary
                res_dict[chr_numb]=print_dict

            else:
                
                print file_name + ": No polymorphisms for coverage cutoff."
        else:
            print file_name + ": No polymorphisms."
            
            
# Use the dictionary over all files to query dbSNP            
if res_dict:
    
    # Only human data in dbSNP, besides only chromosomes/transcriptome ids can be looked up
    if organism == "Homo sapiens":
        res_dict["**header**"] = "Position\tA\tC\tG\tT\tTarget\tMatches in dbSNP\traw A\traw C\traw G\traw T\n"
        
        suc = 0
        fail = 0
        queries=[]
        
        for chr in res_dict.keys():
            
            # Variable header for outfiles dependent on database used for SNP lookup
            if chr == "**header**":
                pass
            else:
                if re.search(r"^\d{1,2}$", chr) or re.search(r"^[X,Y]$", chr) or re.match(r"^[A,N,X][C,G,T,W,Z,M,R][_]", chr):
                    
                    # Create search string for dbSNP. For each chr and pos, the target allele and possible mutations have to be queried
                    for pos in res_dict[chr].keys():
                        target = res_dict[chr][pos][4]
                        target = target.strip("\t").upper()
                        
                        for mut in ["A", "C", "G", "T"]: 
                            if mut != target:
                                queries.append(chr + " " + pos + " "  + "iD" + " " + target + " " + mut + "\n")
                                
                else:
                    processed_f = "calc.nuccounts." + chr_to_enc_dict[chr]
                    print processed_f + ": %s, %s not in dbSNP." %(organism, chr)
                    for pos in res_dict[chr]:
                        res_dict[chr][pos][-5] = "N/A\t"
                    
        if queries:                
            db_dict, isDBerror = getSNPwww(queries)
        
        # The query did not reach the db
        if isDBerror == True:
            
            # If an error occurred, try a second time
            print "Retrying..."
            time.sleep(30)
            db_dict, isDBerror = getSNPwww(queries)
            
            if isDBerror == True:
                for chr_ in res_dict:
                    if chr_ != "**header**":
                        for pos in res_dict[chr_]:
                            res_dict[chr_][pos][-5] ="db:error\t"
                    
        # No database error
        if isDBerror == False: 
                       
            if db_dict:
                
                for chr in db_dict:
                    fail = 0
                    suc = 0
                    
                    for db_pos in db_dict[chr]:
                        db_pos = str(db_pos)
                        
                        if chr in res_dict:
                            
                            # Should be the case, because I queried pos from res_dict
                            if db_pos in res_dict[chr]:
                                rs_counter = 0
                                # List_element is a tupel: (db_rs, db_nuc)
                                db_str = ""
                                
                                for list_element in db_dict[chr][db_pos]:
                                    rs_counter += 1
                                    db_rs = list_element[0]
                                    db_str = db_str + db_rs + ": "
                                    nuc_counter = 0
                                    
                                    for db_nuc in list_element[1]:
                                        nuc_counter += 1
                                        db_str = db_str + db_nuc
                                        
                                        # The last allele should not be followed by a separator
                                        if nuc_counter < len(list_element[1]):
                                            db_str = db_str + "+"
                                            
                                    if rs_counter < len(db_dict[chr][db_pos]):
                                        db_str = db_str + ";"
                                        
                                res_dict[chr][db_pos][-5]= db_str+"\t" #res_dict[chr][db_pos][-10]= db_str+"\t"
                                suc += 1
                            
                            # Request artifact
                            else:  
                                fail += 1
                        
                        # Request artifact
                        else:
                            fail += 1
                            
                    processed_f = "calc.nuccounts." + chr_to_enc_dict[chr]
                    print processed_f + ": dbSNP: %s match(es); %s artifact(s)." % (str(suc), str(fail))
                    
            else:
                sys.stdout.write("dbSNP: 0 matches.\n")



    # Launch subsequent Perl script to analyze the alignment quality 
    print "Analyzing alignment quality: %s" %str(isQualityAnaly)  
    if isQualityAnaly == True:
        
        # Add p-value column to header
        header = res_dict["**header**"]
        header_list = header.split("\t")
        header_list.insert(-4, "P-error (local alignment quality)")
        header = "\t".join(header_list)
        res_dict["**header**"] = header
        
        # Send res_dict as json to perl script
        json_res = json.dumps(res_dict)

        p = Popen(["/bioinf/projects/nanopipe2/calculate/nanopipe_qualreg.pl"], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        result, error = p.communicate(input=json_res)
        print result
        
        if error:
            print "Error: " + error
       

    else:
        # Print Data to separate output files
        # Positions are sorted
        for chr in res_dict.keys():
            if chr != "**header**":
                print_string = res_dict["**header**"]
                sorted_key_list = []
                
                for pos in res_dict[chr].keys():
                    sorted_key_list.append(int(pos))
                sorted_key_list.sort()
                
                for pos in sorted_key_list:
                    print_string = print_string +str(pos) + "\t" + "".join(res_dict[chr][str(pos)])
                    
                print_string = print_string.strip("\n")
                encode = chr_to_enc_dict[chr]
                outfile = curr_dir + "/" + "calc.nuccounts." + encode + ".poly"
                 
                with open(outfile, "w") as output_file:
                    output_file.write(print_string)