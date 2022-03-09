## INTRODUCTION:
- Code supplement for the paper "Optimizing Freight Procurement for Transportation-Inventory Systems Under Supply and Demand Uncertainty" by Lingxiao Wu, Yossiri Adulyasak, and Jean-François Cordeau.
- If you need help using the code, please send an email to lingxiaowu513[at]gmail[dot]com.
- The code and data sets are also available from https://github.com/LingxiaoWu2021/SFPTMP.

## DATASETS:
- The instance data are provided in the folder "data". 
- For a detailed explanation of the instance data, find the "README" file in the folder "data". 

## RESULTS:
- The detailed results obtained by each solution method for all instances are provided in the folder "data". 

## CODE:
- Code of all algorithms used in our computational experiments is provided in the folder "sourcecode". 
- List of .cpp files in the subfolder "src":
  1. S0：   code for S0
  2. S1：   code for S1
  3. S2：   code for S2
  4. S3：   code for S3
  5. CPLEX: code for running CPLEX on model P 
  6. BM1：  code for BM1
  7. BM2：  code for BM2

- List of .h files in the subfolder "inc":
  1. Avgminmax02.h:         user-defined c++ library header file
  2. Seqinsertion.h:         user-defined c++ library header file

## USAGE:
- To run an algorithm for solving an instance:
  1. Copy the instance data from the "data" folder;
  2. Load the data into the code for the algorithm between lines "//input data starts here" and "//input data ends here";
  3. Build and run the code;
  4. Obtain the results from the console window.
