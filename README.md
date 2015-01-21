Antibody2.0
===========

Antibody is a tool we developed for efficient analysis of Sanger sequencing of antibody genes using IMGT V-QUEST. It takes individual .seq files, combines them into a single text file in FASTA format and submits them to V-QUEST. Antibody then parses out the results in a .csv file that can be used for analysis in R.

Antibody requires the Selenium webdriver, Firefox, and Pyperclip. 

I've included a folder of test sequences. You can run it using python antibody2.py output_file mouse_test_seqs/ mouse