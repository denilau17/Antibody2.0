Antibody2.0
===========

Antibody is a tool we developed for efficient analysis of Sanger sequencing of antibody genes using IMGT V-QUEST. It takes individual .seq files and returns their sequence alignment and mutation analysis results in a .csv file suitable for analysis in R. 

Antibody requires the Selenium webdriver, Firefox, and Pyperclip. 

I've included a folder of test sequences. You can run it using python antibody2.py output_file mouse_test_seqs/ mouse