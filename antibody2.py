import os
import pyperclip
import glob
import sys
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import NoSuchElementException
import unittest, time, re

#IMGT accepts a max of 50 sequences per run, so I divided the sequences into batches of 48 since we usually submit half or full 96 well plates for sequencing.

def FASTA_convert(fn, pathway):
    read_files = glob.glob(pathway + "*.seq")
    x = len(read_files)
    num_runs = x/48
    if x%48 != 0:
        num_runs += 1
    print num_runs
    d = {}
    for i in range(1, num_runs+1):
        current_files = read_files[48*(i-1):48+48*(i-1)]
        outfile = open(pathway + fn + str(i) + ".txt", "wb")
        for f in current_files:
            outfile.write(">" + f.split("/")[-1] + "\n")
            for line in open(f, "r"):
                 outfile.write(line)
        outfile.close()

def run_VQUEST(fasta_file, organism):
    g = open(fasta_file, "r")
    pyperclip.copy(g.read())
    g.close()

    fp = webdriver.FirefoxProfile()
    fp.set_preference("dom.max_chrome_script_run_time", 0)
    fp.set_preference("dom.max_script_run_time", 0)
    driver = webdriver.Firefox(firefox_profile=fp)
    driver.get("http://www.imgt.org/IMGT_vquest/vquest?livret=0&Option=" + organism + "Ig")
    driver.find_element_by_xpath("//input[@name='l01p01c05'][2]").click()
    driver.find_element_by_name("l01p01c10").clear()
    driver.find_element_by_name("l01p01c10").send_keys(Keys.COMMAND, 'v')
    driver.find_element_by_name("l01p01c12").click()
    driver.find_element_by_name("l01p01c13").click()
    driver.find_element_by_name("l01p01c06").click()
    driver.find_element_by_name("l01p01c14").click()
    driver.find_element_by_name("l01p01c15").click()
    driver.find_element_by_name("l01p01c16").click()
    driver.find_element_by_name("l01p01c41").click()
    driver.find_element_by_name("l01p01c22").click()
    driver.find_element_by_name("l01p01c23").click()
    driver.find_element_by_name("l01p01c19").click()
    driver.find_element_by_name("l01p01c18").click()
    driver.find_element_by_css_selector("input[type=\"SUBMIT\"]").click()

    output_file = fasta_file[:-4] + "_vquest.txt"
    print "Storing VQuest data to " + output_file + "\n"
    #make an output file and write the html output to that txt file                                     
    a= open(output_file, "w")
    a.write(driver.page_source)
    a.close()
    driver.quit()

def parse_VQUEST(vquest_file, pathway):
    open_file  = open(vquest_file, 'r')
    z  = open_file.read()
    open_file.close()
    parsed_output_name = vquest_file[:-4] + "_parsed.csv"
    print parsed_output_name
    #write headers for csv output file
    d = open(parsed_output_name, 'w')
    d.write('sample, productive, low ID, V gene, V identity, J gene, J identity, D gene, D identity, V nt mut, V s\
ilent mut, V nonsilent mut, FR1 nt mut, FR1 silent mut, FR1 nonsilent mut, CDR1 nt mut, CDR1 silent mut, CDR1 nons\
ilent mut, FR2 nt mut, FR2 silent mut, FR2 nonsilent mut, CDR2 nt mut, CDR2 silent mut, CDR2 nonsilent mut, FR3 nt\
 mut, FR3 silent mut, FR3 nonsilent mut, CDR3 nt mut, CDR3 silent mut, CDR3 nonsilent mut, V aa mut\n')

    #split text file between sequences                                      
    y = z.split('------------------------------')

    all_seq = len(re.findall(r'Sequence\snumber\s\d*\s', z))
    bad_seq = len(re.findall(r'non\sresults', z))
    a = all_seq + bad_seq
    b = 1
    while b <= a:
        find_seqs = re.findall(r'Result\ssummary', y[b])
        non_results = re.findall(r'non\sresults', y[b])

        # identifies successful sequence analysis                                             
        if find_seqs == ['Result summary']:
            #splits off Result summary section                                                
            x = y[b].split('\nResult')

            sample = re.findall(r'\nSequence\snumber\s\d*.*:(.*.\w)', y[b])
            d.write(str(sample) + ',')

            productive = re.findall(r'Productive', y[b]) + re.findall(r'Unproductive', y[b]) + re.findall(r'No rearrangement', y[b])
            d.write(str(productive) + ',')

            low_identity = re.findall(r'Low\sV-REGION', x[1])
            if low_identity == ['Low V-REGION']:
                d.write('low ID,')
            else:
                d.write(' ,')
            V_gene = re.findall(r'V-GENE\sand\sallele;\w*\s(IG.V.*?)\sF', x[1])
            d.write(str(V_gene) + ',')
            V_identity = re.findall(r'V-GENE\w*.*?identity\s\=\s([\d\.]*\%)', x[1])
            d.write(str(V_identity) + ',')
            J_gene = re.findall(r'J-GENE\sand\sallele;\w*\s(IG.J.*?)\sF', x[1])
            d.write(str(J_gene) + ',')
            J_identity = re.findall(r'J-GENE\w*.*?identity\s\=\s([\d\.]*\%)', x[1])
            d.write(str(J_identity) + ',')
            D_gene = re.findall(r'\w*\s(IG.D.*?)\sF', x[1])
            d.write(str(D_gene) + ',')
            D_identity = re.findall(r'D-GENE\w*.*?identity\s\=\s([\d\.]*\%)', x[1])
            d.write(str(D_identity) + ',')

            u = x[1].split('\nAmino')
            V_nt = re.findall(r'V-REGION;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(V_nt) + ',')
            V_nt_silent = re.findall(r'V-REGION;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(V_nt_silent) + ',')
            V_nt_nonsilent = re.findall(r'V-REGION;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(V_nt_nonsilent) + ',')
            FR1_nt = re.findall(r'FR1-IMGT;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR1_nt) + ',')
            FR1_nt_silent = re.findall(r'FR1-IMGT;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR1_nt_silent) + ',')
            FR1_nt_nonsilent = re.findall(r'FR1-IMGT;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR1_nt_nonsilent) + ',')
            CDR1_nt = re.findall(r'CDR1-IMGT;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR1_nt) + ',')
            CDR1_nt_silent = re.findall(r'CDR1-IMGT;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR1_nt_silent) + ',')
            CDR1_nt_nonsilent = re.findall(r'CDR1-IMGT;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR1_nt_nonsilent) + ',')
            FR2_nt = re.findall(r'FR2-IMGT;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR2_nt) + ',')
            FR2_nt_silent = re.findall(r'FR2-IMGT;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR2_nt_silent) + ',')
            FR2_nt_nonsilent = re.findall(r'FR2-IMGT;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR2_nt_nonsilent) + ',')
            CDR2_nt = re.findall(r'CDR2-IMGT;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR2_nt) + ',')
            CDR2_nt_silent = re.findall(r'CDR2-IMGT;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR2_nt_silent) + ',')
            CDR2_nt_nonsilent = re.findall(r'CDR2-IMGT;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR2_nt_nonsilent) + ',')
            FR3_nt = re.findall(r'FR3-IMGT;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR3_nt) + ',')
            FR3_nt_silent = re.findall(r'FR3-IMGT;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR3_nt_silent) + ',')
            FR3_nt_nonsilent = re.findall(r'FR3-IMGT;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(FR3_nt_nonsilent) + ',')
            CDR3_nt = re.findall(r'CDR3-IMGT;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR3_nt) + ',')
            CDR3_nt_silent = re.findall(r'CDR3-IMGT;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR3_nt_silent) + ',')
            CDR3_nt_nonsilent = re.findall(r'CDR3-IMGT;.*?;.*?;.*?;.*?;.*?;(\d*).*?;', u[0])
            d.write(str(CDR3_nt_nonsilent) + ',')
            V_aa = re.findall(r'V-REGION;.*?;.*?;.*?;(\d*).*?;', u[1])
            d.write(str(V_aa) + ',')
            d.write('\n')

        elif non_results == ['non results']:
            sample2 = re.findall(r'Sequence\snumber\s\d*\s(.*.\w);', y[b])
            d.write(str(sample2) + ',')
            d.write('non result' + '\n')

        # deals with random blank spaces                                                     
        else:
            print ""
        b = b + 1
    d.close()

    input = open(parsed_output_name, "rb")
    lines = input.readlines()
    conversion = '[]\''
    newtext = ''
    outputlines = []
    for line in lines:
        temp = line[:]
        for c in conversion:
            temp = temp.replace(c, newtext)
        outputlines.append(temp)

    output = open(parsed_output_name, 'w')
    for line in outputlines:
        output.write(line)
    output.close()

if (__name__ == "__main__"):
    fn = sys.argv[1]
    pathway = sys.argv[2]
    organism = sys.argv[3]    

    FASTA_convert(fn, pathway)
    fasta_files = glob.glob(pathway + "*.txt")
    for i in fasta_files:
        run_VQUEST(i, organism)
    vquest_files = glob.glob(pathway + "*vquest.txt")
    print vquest_files
    for j in vquest_files:
        parse_VQUEST(j, pathway)

