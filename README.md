Antibody2.0
===========

Antibody is a tool we developed for efficient analysis of Sanger sequencing of antibody genes using IMGT V-QUEST. It takes individual .seq files, combines them into a single text file in FASTA format and submits them to V-QUEST. Antibody then parses out the results in a .csv file that can be used for analysis in R.

Antibody requires the use of the Selenium webdriver, Firefox, and Pyperclip. 

To run the program, type python antibody2.py output_filename location_of_files organism

I've included a folder of test files. You can run them using the command:
python antibody2.py test mouse_test_seqs/ mouse