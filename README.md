# Masking NewHope

This project implements Masked NewHope KEM, a cryptographic scheme designed for secure key encapsulation mechanisms (KEM). This repository includes the source code and tools to test the performance of the implementation.
It also contains the source code of masked Kyber copied from https://github.com/fragerar/HOTableConv/tree/main/Masked_KEMs for performance comparision.

The masking idea we use is based on [CGMZ21](https://eprint.iacr.org/2021/1615), please refer to the paper for more details.

## Features
- Implementation of Masked NewHope IND-CPA PKE
- Implementation of Masked NewHope IND-CCA KEM decapsulation
- Performance evaluation tools

## Getting Started

### Clone the Repository
```bash
git clone https://github.com/skad0d0/Masking-NewHope.git
cd Masking-NewHope/code/masked-newhope/src
```

## Running the Speed Test
To evaluate the performance of the Masked NewHope implementation at different masking order, run the following command:
```bash
python3 run_speedtest.py
```
The results are in bench_res.txt

For performance comparision with the Masked Kyber, run:
```bash
python3 run_benchmarks.py
```
The results are in bench_res.txt


