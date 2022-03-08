# hw1_110761504.py
# ©©© Mequanent Argaw Muluneh ©©©

# **** My Final local paths when writing and testing my code *******
# python_path C:\Users\USER\Desktop\BioInfoHW1\hw1_110761504.py
# Input_path: C:\Users\USER\Desktop\BioInfoHW1\mut.txt
# Output_path: C:\Users\USER\Desktop\BioInfoHW1\output.txt

# Copy the following command to run the code, with the respective path changes
# python C:\Users\USER\Desktop\BioInfoHW1\hw1_110761504.py --input C:\Users\USER\Desktop\BioInfoHW1\mut.txt --pam 250 --output C:\Users\USER\Desktop\BioInfoHW1\output.txt
# The arguments (--input input_path_mut.txt --pam x --output output_path_pamx.txt) can be shuffled

import sys # to forward the outputs into file
import argparse # For parsing arguments when they are shuffled

import numpy as np
np.set_printoptions(suppress=True) # To suppress displaying the scientific notation

# ****************************************************************************************
# ******** Arguments' parsing section to handle their shuffling **************************
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input_path_mut.txt") # default type: string
parser.add_argument("--output", help = "output_path_pamx.txt")
parser.add_argument("--pam", help = "x for PAMx", type = int)

args = parser.parse_args()

input_path = args.input
output_path = args.output
x = args.pam
# ***************** End of parsing *************************************************************
# **********************************************************************************************

sys.stdout=open(output_path,"w") # Open an output file to write all printed statements and values  

# **********************************************************************************************
# ************ Preprocess the input data, mut.txt, to get the PAM1 matrix ***********************
input_file = open(input_path, "r")
lines = input_file.readlines()[2:] # remove the first two lines
input_file.close()

mut = []
for line in lines:
  line = line[1:].strip().split() # remove the letters from left and the white spaces
  mut.append(line)

mut = np.array(mut).astype(int) # convert the string array into numpy array of integers
print('\nChecking input mut size: ', mut.shape) # confirm if we get the right matrix size from the text
PAM1 = np.true_divide(mut, 10000) # Change the four digit values into values between 0 & 1
# Normalized frequencies of acids keeping the order in the given matrix
fi = [0.087, 0.041, 0.040, 0.047, 0.033, 0.038, 0.050, 0.089, 0.034, 0.037, 0.085, 0.081,
      0.015, 0.040, 0.051, 0.070, 0.058, 0.010, 0.030, 0.065]

print("\nPrinted values include: PAM1, PAMx(round to nearest int) and Log-odds of PAMx as nearest int")

print("\n*******************************************************")
print("****************** PAM1 *******************************")
print("*******************************************************")
print(PAM1)
# ************** End of pre-processing the required inputs  *******************************
# *****************************************************************************************

# ******************  The PAMx Function  **************************************************

def pam(x):
  PAM_str = "PAM" + str(x) # To make x in the PAMx flexible when displaying
  if x < 1:     # Check if the value of x is valid
    print("You entered {}, please enter a number >= 1".format(x))
    sys.exit()
  elif x == 1:
    return PAM1
  else:
    PAMx = PAM1
    for i in range(2, x+1):
      PAMx = np.matmul(PAMx, PAM1) # np.matmul(mut, PAMx) also works the same
    return PAMx, PAM_str

# ************** End of the PAMx Function ******************************************************
# **********************************************************************************************

PAMx, PAM_str = pam(x) # get the PAMx value and 'PAMx' for displaying purpose


print("\n*******************************************************")
print("****************** {} *****************************".format(PAM_str))
print("*******************************************************\n")
print(np.rint(PAMx * 100).astype(int)) # Part of post-processing to display PAMx values as nearest -
                              # integers resolving formats like 1.0, 0.0 changed to 1, 0 respctively

# **********************************************************************************************
# ************************  The Log-odds matrix function   *************************************
def log_odds():
  S = []    
  for i in range(20):
    S.append(PAMx[i] / fi[i]) # for equation Rij = Mij / fi
  S = np.array(S) # Here S holds Rij values for decreasing complexity for later use

  for i in range(20):       # for loop is used rather than numpy since this is
    for j in range(20):     # more suitable to handle log10(0) exceptions. 
      if S[i][j] > 0:
        S[i][j] = 10 * np.log10(S[i][j])
      else:
        S[i][j] = 0 # To set log10(0) = 0

  # The log-odds matrix should be made symmetric
  for i in range(20):
    for j in range(20):
      S[i][j] = (S[i][j] + S[j][i]) / 2
      S[j][i] = S[i][j]

  return np.rint(S).astype(int) # Part of the post-process to round into nearest integer and
                                # change into int to avoid decimal formats like 2.0, 0.0 ...

# ************************** End of the Log-odds Function  *************************************
# **********************************************************************************************

print("\n*********************************************************")
print("****** The Log-odds matrix for {} *******************".format(PAM_str))
print("*********************************************************")
print(log_odds()) # The post-processed values streamed to the output file

sys.stdout.close() # Close the output file

