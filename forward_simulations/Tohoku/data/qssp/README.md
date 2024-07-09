To reproduce results from Zhang et al 2020, do the following steps:


### (1) Compiling QSSP_PEGS
1. Get the 2020 version of QSSP_PEGS from:
ftp://ftp.gfz-potsdam.de/home/turk/wang/

2. Correct one of the source files.

   +  Navigate to `SourceCode`
   + Open `qpshkern.f`
   + Change **line 29** from
      ```fortran
      if(ldeg.le.1.or.lyr.lt.lyup.or.lyr.gt.lylw.or.lylw.lt.lylwa)return
      ```
      to 
      ```fortran
      if(ldeg.le.0.or.lyr.lt.lyup.or.lyr.gt.lylw.or.lylw.lt.lylwa)return
      ```
         
   + Note that this bug probably makes no difference. It was communicated to me by Rongjiang Wang for a different application but is worth doing anyway. 
      
3. Make the executable. Still in the `SourceCode` directory  by just typing `make`. 
   - If it throws an error ensure you do not already have a `../bin` directory. 
   - You may also get an error if you are in your `(base)` conda environment. Try `conda deactivate` and then `make`. 

### (2) Running QSSP_PEGS

1. Ensure you have `tohoku_PEGS_input_Eaton.inp` in your `qssp_PEGS` directory. 
2. Run the code with 
```bash 
     $ ./bin/qssp tohoku_PEGS_input_Eaton.inp
```
3. If it prompts you to enter the input file just re-enter `tohoku_PEGS_input_Eaton.inp` 

It should take ~ 1 hour to run on a laptop if it is the first time you run it since it need
to build the Green's function database.  

The output files will have the format `ps_vapf_0d__******_.dat`

### (3) Post-processing
1. Once you have your outputs from QSSP we need them in a SAC format. First, copy them to
the `prem_elastic/raw' directory. 
2. Next, run the script `postprocess.sh` 

Note that you will need to change the paths in your post process shell file for your own 
SAC binary path and also finding a way to run the python script - e.g. with a conda 
environment. 


