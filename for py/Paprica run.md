# For paprica on linux. 


```bash

#Full PAPRICA workflow

#### Step 1: Paprica installation (Python)

git clone https://github.com/bowmanjeffs/paprica.git
cd paprica
chmod a+x *py
chmod a+x *sh


### the python package
conda create -n paprica_env python=3.10.12
pip install pandas
pip install biopython 
pip install joblib
pip install termcolor
pip install ete3
pip install seqmagick
conda install -c bioconda epa-ng -y
conda install -c bioconda gappa -y
conda install -c bioconda infernal -y
conda install -c bioconda hmmer -y 
conda install -c easel


## test
./paprica-run.sh test bacteria


#### Step 2: individual excel files of seqtab for each sample  (R)
library(dplyr)
df<- readRDS(#your ps object) # input seqtab file, please change accordingly

a <- df %>% t %>% data.frame # convert to dataframe, transposed

for (x in 1:ncol(df)) { 
  b <- a[x]
  b
  c <- cbind(rownames(b), b)
  colnames(c) <- c('sequence', 'abundance')
  write.csv(c, paste0(colnames(a)[x], '.csv'), row.names = F)
}

# you should have individual csv files of all your samples from your seqtab file

#### Step 3: Looping all CSV into Fasta files (R)
all_files <- list.files()
csv_files <- all_files[grepl("\\.csv$", all_files)]

for (i in csv_files){
  #print(i)
  
  df = read.csv(i, header=T, row.names=1)
  #head(df,2)
  
  m <- as.matrix(df)
  csum <- colSums(m)
  idx <- unlist(lapply(csum, seq_len), use.names=FALSE)
  res <- matrix(c(sprintf(">%s_%d", rep("seq", csum), idx), # id
                  rep(rownames(m)[row(m)], m)),                   # sequence
                nrow=2, byrow=TRUE)
  
  out <- paste0(substr(i, 1, nchar(i) - 4), ".fasta")
  
  writeLines(res, out)
  
}


# you should have individual fasta files after this step. 


#### Step 4: PAPRICA

cd paprica
conda activate paprica_env
# I would suggest putting all your fasta files into a file and enter it

ls -1 | sed -e 's/\.fasta$//' > samples.txt #this would create a txt file of all your sample names

cd paprica #return to your paprica folder
nano my_analysis.sh #create a paprica loop command 
#!/bin/bash
while read f;do
   ./paprica-run.sh $f bacteria
done < samples.txt  #runs and loops all the samples in sample.txt into paprica


#### Step 5: Combining all paprica results into csv files

paprica-combine_results.py -domain bacteria -o 20240508out # change the name after -o to your preferred name for your files. 

```