# AntonymDetection
Implementation of Word Embedding-based Antonym Detection using Thesauri and Distributional Information in NAACL2015

This implementation is based on GloVe and ivLBL, and uses the SpaceSaving algorithm.  Please refer to the following papers for detail.
* Masataka Ono, Makoto Miwa, and Yutaka Sasaki. "Word Embedding-based Antonym Detection using Thesauri and Distributional Information", NAACL/HLT-2015, 2015.
* Metwally, Ahmed, Divyakant Agrawal, and Amr El Abbadi. "Efficient computation of frequent and top-k elements in data streams." Database Theory-ICDT 2005. Springer Berlin Heidelberg, 2005. 398-412.
* Jeffrey Pennington, Richard Socher, and Christopher D. Manning. "GloVe: Global Vectors for Word Representation.", EMNLP 2014, 2014.
* Andriy Mnih and Koray Kavukcuoglu. Learning word embeddings efficiently with noise-contrastive estimation. In NIPS 26, 2013.

## Prerequisites
* [Eigen](http://eigen.tuxfamily.org/)
* [gzstream](http://www.cs.unc.edu/Research/compgeom/gzstream/)
* [boost](http://www.boost.org/)
* g++ 4.8 or more
* 2GB or more RAM 
* 10GB or more disk (temporal files are stored in tmp/)
* capable of running 10 or more threads (Please change the code if required.)

## Compilation:

`make`

## Usage:

### Use a pretrained model

```
wget http://tti-coin.jp/data/antonym_model.tar.gz
tar xzf antonym_model.tar.gz
./TestAntonym model/100000.wc.bin model/we-td.model data/antonym/devset.txt
```

### Training a new model

1. Prepare Data
   * put gzipped text files (each sentence per line) into a directory (e.g., texts)
   * put (replace) antonym dictionary into a data/dict directory
2. Count top-N frequence words (with MINIMUM_COUNT) 
   
   ```
   ./WordCount TEXT N FREQFILE [MIN_COUNT]
   (e.g., ./WordCount texts/ 100000 100000.wc.bin 0)
   ```
   
3. Count co-occurrences in a sentence 
   
   ```
   ./SequenceCount DATA FREQFILE COOCCFILE WINDOW_SIZE WEIGHT
   (WEIGHT  2: ivlbl, 1: glove, 0: no weighting)
   (e.g., ./SequenceCount data 100000.wc.bin 100000.seq.bin 5 2) 
   ```
   
4. Train Antonym model on co-occurrences and dictionary 
   
   ```
   ./TrainAntonym FREQFILE COOCCFILE MODEL DIM ITER
   (e.g., ./TrainAntonym 100000.wc.bin 100000.seq.bin we-td.model 100 20)
   ```
   
5. Test Antonym model 
   
   ```
   ./TestAntonym FREQFILE MODEL TEST_FILE
   (e.g., ./TestAntonym 100000.wc.bin we-td.model data/antonym/devset.txt)
   ```
   

## Notes
Please cite our paper when using this tool.
* Masataka Ono, Makoto Miwa and Yutaka Sasaki. Word Embedding-based Antonym Detection using Thesauri and Distributional Information. In the Proceedings of the 2015 Conference of the North American Chapter of the Association for Computational Linguistics - Human Language Technologies (NAACL-HLT 2015). pp. 984--989, June 2015. 

