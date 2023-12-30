import tkinter as tk
import time
import numpy as np
from tkinter import ttk
from tkinter import scrolledtext
from tkinter import filedialog
from tkinter import messagebox

st = time.time()

def process_data():
    # Extract data from the text input
    fasta_data = fasta_text.get("1.0",tk.END).split('\n')
    # Extract selected option from the dropdown menu
    selected_option = option_var.get()
    # FASTA text reader
    new_dict = {}
    for item in fasta_data:
        if item.startswith('>'):
            key = item
            new_dict[key] = []
        else:
            new_dict[key].append(item)
    for key,value in new_dict.items():
        value = ''.join(value)
        new_dict[key] = value
    value = []
    for values in new_dict.values():
        value.append(values)

    return value

def check_seq():
    seq = process_data()[0]
    bases = ['A','T','G','C','N']
    count_ukw = 0
    count_upper = 0
    count_len = 0
    
    for i in range(0,len(seq)):
     
        if seq[i] not in bases:
            count_ukw += 1
            
        elif not seq[i].isupper():
            count_upper += 1
            
        elif len(seq) < 11:
            count_len += 1
           
    if count_ukw > 0:
        messagebox.showinfo("InValid Sequnece","InValid sequence: Unknown Base")
        action_button.config(state=tk.DISABLED)
    elif count_upper > 0:
        messagebox.showinfo("Invalid sequence","Please make the sequence all capital")
        action_button.config(state=tk.DISABLED)
    elif count_len > 0:
        messagebox.showinfo("Invalid sequence","Please enter a seqeuence length greater than 11")
        action_button.config(state=tk.DISABLED)
    else:
        messagebox.showinfo("Valid sequence","Valid sequence")
        action_button.config(state=tk.NORMAL)
        
        
def upload_file():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta")])
    if file_path:
        with open(file_path,'r') as file:
            fasta_text.delete("1.0", tk.END)
            fasta_text.insert(tk.END, file.read())


    
#Blast algorithm  
def fasta_reader(file_name):
    with open(file_name) as f:
        new_dict = {}
        
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                
                key = line.lstrip('>')
                new_dict[key] = []
            else:
                
                new_dict[key].append(line)
    for key,value in new_dict.items():
        value = ''.join(value)
        new_dict[key] = value
    value = []
    for values in new_dict.values():
        value.append(values)
    return value
     
def reverse_complement(seq):
    rev_comp_seq = ''
    seq = seq[::-1]
    for i in range(0,len(seq)):
        if seq[i] == 'A':
            rev_comp_seq += 'T'
        elif seq[i] == 'T':
            rev_comp_seq += 'A'
        elif seq[i] == 'G':
            rev_comp_seq += 'C'
        else:
            rev_comp_seq += 'G'
    return rev_comp_seq
            
def edit_dist(query,seq):
    matrix = np.zeros((len(query)+1,len(seq)+1))
    for i in range(0,len(seq)+1):
        matrix[i][0] = i 
        matrix[0][i] = i 
    for i in range(1,len(seq)+1):
        for j in range(1,len(query)+1):
            left = matrix[i][j-1]
            right =  matrix[i-1][j]
            diag = matrix[i-1][j-1]
            minim  =  min(left,right,diag)
            if seq[i-1] == query[j-1]:
                matrix[i][j] = minim
            elif seq[i-1] != query[j-1]:
                matrix[i][j] = minim + 1
    return matrix[len(seq)][len(query)]

def kmer(query,w=11):
    hash_map = {}
    for i in range(0,len(query)-w+1):
        sub_seq = query[i:i+w]
        if sub_seq in hash_map:
            hash_map[sub_seq].append(i)
        else:
            hash_map[sub_seq] = [i]
    return hash_map
def hits(seq,query_map,w=11):
    k_index = []
    for i in range(0,len(seq)- w + 1):
        seq_mer = seq[i:i+w]
        if seq_mer in query_map:
            l = query_map[seq_mer]
            #print(l)
            for ind in l:
                k_index.append((ind,i))
    return k_index
    
def hamming_dist(seq,query):
    count = 0
    length = len(query)
    for i in range(0,len(seq)):
        if seq[i] == query[i]:
            pass
        elif seq[i] != query[i]:
            count += 1
    percent = (count/length) * 100
    return percent
    
def extend_hits(query,seq,indexes,w=11):
    index_q = []
    index_s = []
    answers = []
    threshold_score= option_var_threshold.get()
    sensitivity= option_var_sensitivity.get()
    print("Function called")
    for i in range(0,len(indexes)):
        index_q.append(indexes[i][0])
        index_s.append(indexes[i][1])
    
    for i in range(0,len(index_q)):
        
        index_f_q = index_q[i]
        index_b_q = index_q[i]
        
        index_f_s = index_s[i]
        index_b_s = index_s[i]
        score = 0
        count = 0
        
        
        
        if index_q[i] >-1:
            while (index_f_q < (len(query)-1) and index_f_s < (len(seq)-1)) and (score < threshold_score):
                score = hamming_dist(seq[index_b_s:index_f_s+1],query[index_b_q:index_f_q+1])
                
                index_f_q += 1
                index_f_s += 1
                if index_b_q > 0 and index_b_s > 0:
                    index_b_q -= 1
                    index_b_s -= 1
                elif index_b_q <= 0:
                    index_b_q = 0
                if index_b_s <= 0:
                    index_b_s = 0
                
                count += 1
                #print(count)
        if (score < threshold_score) and (abs(index_q[i] - index_f_q) > len(query) - sensitivity):
            answers.append([index_q[i],index_s[i],index_f_q,index_f_s])
        
    with open("answers.txt",'w') as f:
        for i in answers:
            f.write("Index of the query sequence that match: \n")
            f.write(str(i[0]))
            f.write("\t")
            f.write(str(i[2]))
            f.write("\n")
            f.write("Index of the database sequence that matches the query: \n")
            f.write(str(i[1]))
            f.write("\t")
            f.write(str(i[3]))
            f.write("\n")
            f.write("\n")
            f.write("\n")
            
            
        
    return answers        

def process():
    query = process_data()[0]
    selected_option = option_var.get()
     
    if selected_option == "Staphylococcous":
        # Action for staph
        db_seq = fasta_reader('sequence.fasta')[0]
        rev_seq = reverse_complement(db_seq)
        hash_map = kmer(query)
        hitu = hits(rev_seq,hash_map)
        print(extend_hits(query,rev_seq,hitu))
        et = time.time()
        elapsed_time = et-st
        print("Execution time: ",elapsed_time,"seconds")
        messagebox.showinfo("Job Done","Job Done Subsequence found, file written")

    elif selected_option == "E.coli":
        # Action for Ecoli
        db_seq = fasta_reader('sequence.fasta')[0]
        rev_seq = reverse_complement(db_seq)
        hash_map = kmer(query)
        hitu = hits(rev_seq,hash_map)
        print(extend_hits(query,rev_seq,hitu))
        et = time.time()
        elapsed_time = et-st
        print("Execution time: ",elapsed_time,"seconds")
        messagebox.showinfo("Job Done","Job Done Subsequence found, file written")

    elif selected_option == "Actinobacter":
        # Action for Actinobacter
        db_seq = fasta_reader('sequence.fasta')[0]
        rev_seq = reverse_complement(db_seq)
        hash_map = kmer(query)
        hitu = hits(rev_seq,hash_map)
        print(extend_hits(query,rev_seq,hitu))
        et = time.time()
        elapsed_time = et-st
        print("Execution time: ",elapsed_time,"seconds")
        messagebox.showinfo("Job Done","Job Done Subsequence found,file written")

            
       

root = tk.Tk()
root.title("Tkinter Application")

root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)

#FASTA format input
fasta_text = scrolledtext.ScrolledText(root, width=50, height=10, wrap=tk.WORD)
fasta_text.insert(tk.END, ">Example 1 \nACGT")
fasta_text.grid(column=0, row=0,sticky="nsew")

#Upload file button
upload_button = tk.Button(root, text="Upload File", command=upload_file)
upload_button.grid(column = 0, row = 1)

#dropdown menu
options = ["Staphylococcous", "E.coli", "Actinobacter"]
option_var = tk.StringVar(root)
option_var.set(options[0])  # Set default option
option_menu = tk.OptionMenu(root, option_var, *options)
option_menu.grid(column=0,row=2)

#Dropdown menu
options_threshold = [95, 90, 85]
option_var_threshold = tk.IntVar(root)
option_var_threshold.set(options_threshold[0])  # Set default option
option_menu_threshold = tk.OptionMenu(root, option_var_threshold, *options_threshold)
option_menu_threshold.grid(column=1,row=2)

#Dropdown menu
options_sensitivity = [50, 40, 30]
option_var_sensitivity = tk.IntVar(root)
option_var_sensitivity.set(options_sensitivity[0])  # Set default option
option_menu_sensitivity = tk.OptionMenu(root, option_var_sensitivity, *options_sensitivity)
option_menu_sensitivity.grid(column=2,row=2)

#check data button
check_button = tk.Button(root, text="Check Data", command=check_seq)
check_button.grid(column=0,row=4)


#Button for chossing database sequence
action_button = tk.Button(root, text="Process", command=process)
action_button.config(state=tk.DISABLED)
action_button.grid(column=0,row=6)


# Run the Tkinter event loop

root.mainloop()
